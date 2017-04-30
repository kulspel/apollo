clear all
close all
clc
%%
% number of degrees of freedom per node
ndf = 1;  % temperature is a scalar!

% obtain coordinates and element properties
[ndof, nel, ex, ey, edof]=input_Apollo_CALFEM(ndf);

% plot mesh, element numbers and node number
eldraw2(ex,ey,[1,2,1])

%%
% introduce the needed material parameters 

%thermal conductivity of titanium 
%22 W/(m ·K)
k=22;

%mass density of titanium
%4540 kg/m3
rho=4540;

%mass speci?c heat capacity of titanium 
%520 J/(kg·K)
c=520;



%%
%Other usefule parameters

%convection coefficients
%alpha_1 = 120 W/(m2 ·K), alpha_2 = 20 W/(m2 ·K) 
alpha1 = 120;
alpha2 = 20;

%Far field temperature
%T1? = 500?C ,T2? = ?50?C 
T1_inf=500;
T2_inf=-50;

%Thickness
%1m
t=1;
ep=t;

%Volume heat source
eq=0;

%%
% introduce the constitutive matrix D
D=-k*eye(2);

% initialise global tangent stiffness matrix K and global force vector
C=zeros(ndof);
K=zeros(ndof);
f=zeros(ndof,1);

%Boundaries
[NodesAB,NodesBC,NodesCD,NodesDA,NodesTubes] = nodes_2tria3();
%[NodesAB,NodesBC,NodesCD,NodesDA,NodesTubes] = BCNodes()

upper_boundary = horzcat(NodesBC,NodesCD,NodesDA);

keyboard

for e = 1:nel
    %C-matrix for shape functions
   % C=[1 xi yi;
   %    1 xj yj;
   %    1 xk yk];
    % call plantml to calculate the element heat capacity matrix for element e
    Ce=plantml(ex(e,:),ey(e,:),c*rho);
    
    % call flw2te to calculate the element stiffness matrix for element e
    [Ke,fe] =flw2te(ex(e,:),ey(e,:),ep,D,eq);
    % call assem to assamble the element stiffness matrix to the global one
    [K,f]=assem(edof,K,Ke,f,fe);
end


% introduce the boundary vector fb
fb=zeros(ndof,1);

%%
% call solveq
%[a,Q]=solve(K,fb,bc);
keyboard

%%
%===== start of the main FEM program =====================================
% separate the degrees of freedoms into
% 1. freeDofs, the unkown node temperature we need to compute
% 2. drltDofs, prescribed temperature values at the boundary

% a vector contains all dofs
allDofs = (1:1:nnp*ndf)';

% dofs with Dirichlet boundary conditions, prescribed temperature values
numDrltDofs = length(NodesTubes);
drltDofs = zeros(numDrltDofs,1);
drltValues = zeros(numDrltDofs, 1);

%Dirichelet BC. Prescribed temperature in cooling pipes
Tg = 100;

for i=1:numDrltDofs
    node = NodesTubes(i);
    drltDofs(i) = (node-1)*ndf + 1;
    drltValues(i) = Tg;
end

% freeDofs = allDofs "-" drltDofs
freeDofs = setdiff(allDofs, drltDofs);

% dofs with Robin boundary conditions
numRobinDofs = length(NodesAB);

lower_robinDofs = zeros(numRobinDofs,1);
lower_robinValues = zeros(numRobinDofs, 1);

%Convecion coefficient for AB boundary (W/(m2 ·K))
alpha1=120;
%Celsius
T1_inf=500;

for i=1:numRobinDofs
    node = NodesAB(i);
    lower_robinDofs(i) = (node-1)*ndf + 1;
    lower_robinValues(i) = -alpha1*T1_inf;
end

numRobinDofs = numRobinDofs+length(upper_boundary);

upper_robinDofs = zeros(numRobinDofs-length(NodesAB),1);
upper_robinValues = zeros(numRobinDofs-length(NodesAB), 1);

%Convecion coefficient for upper boundary (W/(m2 ·K))
alpha2=20;
%Celsius
T2_inf=-50;

for i=1:numRobinDofs-length(NodesAB)
    node = upper_boundary(i);
    upper_robinDofs(i) = (node-1)*ndf + 1;
    upper_robinValues(i) = -alpha2*T2_inf;
end

robinDofs = vertcat(lower_robinDofs,upper_robinDofs);
robinValues = vertcat(lower_robinValues,upper_robinValues);

clear NodesAB NodesBC NodesCD NodesDA NodesTubes upper_boundary upper_robinDofs upper_robinValues lower_boundary lower_robinDofs lower_robinValues node


%=========================================================================
% fe-analysis
%=========================================================================

% unknown node temperatures
a = zeros(ndf*nnp,1);

% initialise global volume force vector
fvol = zeros(ndf*nnp,1);