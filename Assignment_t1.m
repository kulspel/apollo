clear all
close all
clc
% obtain coordinates and element properties
[x, elem] = input_2tria3();
%[x, elem, nel] = input_Apollo();

% number of spatial dimensions
ndm = 2; % 1d heat equation
% number of degrees of freedom per node
ndf = 1; % temperature is a scalar
% number of nodes points in the mesh
nnp = length(x);

% plot mesh, element numbers and node number
%Assignment_plot(elem, x);

[Ex,Ey]=coordxtr(elem.cn,x,Dof,nen)
eldraw2(x(:,1)',x(:,2)')

[NodesAB,NodesBC,NodesCD,NodesDA,NodesTubes] = nodes_2tria3();
%[NodesAB,NodesBC,NodesCD,NodesDA,NodesTubes] = BCNodes()

upper_boundary = horzcat(NodesBC,NodesCD,NodesDA);


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