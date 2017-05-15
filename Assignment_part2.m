clear all
close all
clc

%%
% number of degrees of freedom per node
ndf = 2;  % 2-dimensional

% obtain coordinates and element properties (needs ndf!)
input_Apollo;

keyboard

%Boundaries
BCNodes;

nen = size(elem(1).cn,2);

fig = figure;
% plot mesh, element numbers and node number
eldraw2(ex,ey,[1,2,1])

%%
% introduce the needed material parameters

%Poisson's ratio nu= 0.36
nu=0.36;

%Youngs modulus E=110 GPa
E=110*10^9;

%mass density of titanium
%4540 kg/m3
rho=4540;

%Other useful parameters

%Maximum pressure pmax=0.1 MPa
pmax = 0.1*10^6;

%Gravity acceleration constant g=9.81 m/s^2
g=9.81;

%Thickness
%1m
t=1;

%%

%Polynomial approx for the force
%y values goes from 0 to pmax to 0
Y=[0 pmax 0];

%x-values goes from the minimum possible x value to he maximum possible
X=linspace(min([min(ex(:,1)) min(ex(:,2)) min(ex(:,3))]), max([max(ex(:,1)) max(ex(:,2)) max(ex(:,3))]),3);
P=polyfit(X,Y,2);

clear Y X

%Body force, gravity acting ne;gative in the y-dir
eq=[0; -g*rho];

% Analysis type:Plane stress, ptype=2
ptype=2;

ep=[ptype t];

% introduce the constitutive matrix D
D = hooke(ptype,E,nu);

% initialise tangent stiffness matrix K and global force vector
K=zeros(ndof);
f=zeros(ndof,1);

%keyboard

%%
%Element loop
for e = 1:nel
    
    % call plante to calculate the element stiffness matrix for element e
    [Ke,fe] =plante(ex(e,:),ey(e,:),ep,D,eq);
    %keyboard
    %%
    
    %Counters to see how many nodes are on each boundary
    nAB = zeros(1,nen);
    %nTubes = 0;
    
    
    %Loop over the nodes in the element and see if at least 2 of them are
    %on a boundary
    for n=1:nen
        
        %Picking the node we are currently looking at
        node = elem(e).cn(n);
        %keyboard
        
        %Check if this node is on the AB boundary
        %The following steps are repeated for every boundary except the
        %tubes
        if (any(node==NodesAB)==1)
            
            %If it is on the boundary increase the counter, since we only
            %care if 2 nodes are on the boundary
            nAB(1)=nAB(1)+1;
            
            %THis is to be able to pick out which of the 2 nodes are on the
            %boundary since the calculations look slightly different
            %depending on which node is not on the boundary
            if(nAB(2)==0)
                nAB(2)=node;
            else
                nAB(3)=node;
            end
        end
    end
    %end for the boundary loop
    
    
    %%
    %keyboard
    fe_additionAB=0;
    
    %Now that we know how many nodes are on each boundary we can send them
    %into the boundary conditions calculator
    if(nAB(1)==2)
        %keyboard
        fe_additionAB=neumann(ex(e,:),ey(e,:),t,P,find(elem(e).cn==setdiff(elem(e).cn,nAB(2:end))));
        %keyboard
    end
    
    clear n nAB
    fe_addition=fe_additionAB;
    
    clear fe_additionAB
    %keyboard
    
    fe=fe+fe_addition;
    keyboard
    
    % call assem to assamble the element stiffness matrix to the global one
    [K,f]=assem(edof(e,:),K,Ke,f,fe);
    keyboard
end
%End of element loop
clear Ke Ke_addition fe fe_addition nu E t ep eq

%%
%Set Dirichelet BC
bc = vertcat([nodedof(NodesCD(:),3),zeros(size(NodesCD,2),1)],[nodedof(2,3),0]);

clear NodesAB NodesBC NodesCD NodesDA NodesTubes

% Call solveq
a=solveq(K,f,bc);

%%
%Post processing

%Fix post processing, needs to calculate von mises stress
% call extract
ed=extract(edof,a);
%keyboard

Seff_el=zeros(nel,1);

%Element loop
for e = 1:nel
    [es,~]=plants(ex(e,:),ey(e,:),ep,D,ed);
    Seff_el(e)=sqrt(es(1)^2+es(2)^2+3*es(3)^2);
end

Seff_nod=zeros(nnp,1);

for i=1:nnp
    [c0,c1]=find(conn(:,2:4)==i);
    Seff_nod(i,1)=sum(Seff_el(c0))/size(c0,1);
end

keyboard

%fill(ex',ey',ed')
fill(ex',ey',Seff_nod')

%Settings for the graphical presentation
colormap('jet')
h=colorbar;
h.Label.String='Von mises stress';
%caxis([-50 500]);
title(strcat('Temperature distribution after ',{' '}, int2str(time), ' seconds'));

%Saving of the image
% saveas(fig,strcat('Time_',int2str(time),'_timestep',int2str(dt),'_step',int2str(step)),'png')