clear all
close all
clc

%%
% number of degrees of freedom per node
ndf = 2;  % temperature is a scalar!

% obtain coordinates and element properties (needs ndf!)
input_Apollo;

%Boundaries
BCNodes;

nen = size(elem(1).cn,2);

fig= figure;
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

%Gravity acceleration constant g=9.81 m/s^2
g=9.81;

%Thickness
%1m
t=1;

%%

%Body force, gravity acting negative in the y-dir 
eq=[0; -g*rho};

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
    Ke_additionAB=0;
    fe_additionAB=0;
    
    %Now that we know how many nodes are on each boundary we can send them
    %into the boundary conditions calculator
    if(nAB(1)==2)
        %keyboard
        [Ke_additionAB,fe_additionAB]=elem_cond(ex(e,:),ey(e,:),ep,alpha1,T1_inf,find(elem(e).cn==setdiff(elem(e).cn,nAB(2:end))));
        %keyboard
    end
    
    clear n nAB 
    
    Ke_addition=Ke_additionAB;
    fe_addition=fe_additionAB;
    
    clear Ke_additionAB fe_additionAB 
    %keyboard
    
    Ke=Ke+Ke_addition;
    fe=fe+fe_addition;
    %keyboard
    
    % call assem to assamble the element stiffness matrix to the global one
    [K,f]=assem(edof(e,:),K,Ke,f,fe);
    %keyboard
end
%End of element loop
clear Ke Ke_addition fe fe_addition nu E t ep eq

%%
%Set Dirichelet BC and initial conditions
bc = [nodedof(NodesTubes(:),2),ones(size(NodesTubes,2),1)*Tg];

clear NodesAB NodesBC NodesCD NodesDA NodesTubes


% loop over all time/load steps
while (time < stopTime)
    time = time + dt;
    step = step + 1;
    
    %keyboard
    %use [a,Q]=solve(K,f,bc); for steady state
    % call solveq
    %[a,Q]=solve(K,f,bc);
    
    
    ftot=f+(C/dt)*a;
    
    [a,Q]=solve(K,ftot,bc);
    %a_pre=a;
    %keyboard
    
    fprintf(1,'step= %2d time= %8.4e dt= %8.4e\n', step, time, dt);
    
    %Saving of a picture every 100 seconds
    %Activating this increases runtime alot!
    %     if(mod(time,100)==0)
    %         %Post processing
    %         % call extract
    %         ed=extract(edof,a);
    %         %keyboard
    %
    %         %fill(ex',ey',ed')
    %         fill(ex',ey',ed')
    %
    %         %Settings for the graphical presentation
    %         colormap('jet')
    %         h=colorbar;
    %         h.Label.String='Temperature in centigrade';
    %         caxis([-50 500]);
    %
    %         title(strcat('Temperature distribution after ',{' '}, int2str(time), ' seconds'));
    %         saveas(fig,strcat('Time_',int2str(time),'_timestep',int2str(dt),'_step',int2str(step)),'png')
    %     end
    
end
%End of time loop

clear f ftot K Ktot

%%
%Post processing
% call extract
ed=extract(edof,a);
%keyboard

%fill(ex',ey',ed')
fill(ex',ey',ed')

%Settings for the graphical presentation
colormap('jet')
h=colorbar;
h.Label.String='Temperature in centigrade';
caxis([-50 500]);
title(strcat('Temperature distribution after ',{' '}, int2str(time), ' seconds'));

%Saving of the image
% saveas(fig,strcat('Time_',int2str(time),'_timestep',int2str(dt),'_step',int2str(step)),'png')