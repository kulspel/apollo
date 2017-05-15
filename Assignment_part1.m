clear all
close all
clc

%%
% number of degrees of freedom per node
ndf = 1;  % temperature is a scalar!

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

%thermal conductivity of titanium
%22 W/(m ·K)
k=22;

%mass density of titanium
%4540 kg/m3
rho=4540;

%mass speci?c heat capacity of titanium
%520 J/(kg·K)
c=520;

%Other usefule parameters

%convection coefficients
%alpha_1 = 120 W/(m2 ·K), alpha_2 = 20 W/(m2 ·K)
alpha1 = 120;
alpha2 = 20;

%Far field temperature
%T1? = 500?C ,T2? = ?50?C
T1_inf=500;
T2_inf=-50;

%Tubetemp deg celsius
Tg=100;

%Initial temperature T_0=-50 deg celsius
T_0=-50;

%Thickness
%1m
t=1;
ep=t;

%Volume heat source
eq=0;


%%
% introduce the constitutive matrix D
D=k*eye(2);

% initialise global heat capacity matrix, tangent stiffness matrix K and global force vector
C=zeros(ndof);
K=zeros(ndof);
f=zeros(ndof,1);

%keyboard

%%
%Element loop
for e = 1:nel
    % call plantml to calculate the element heat capacity matrix for element e
    Ce=plantml(ex(e,:),ey(e,:),c*rho);
    
    % call flw2te to calculate the element stiffness matrix for element e
    [Ke,fe] =flw2te(ex(e,:),ey(e,:),ep,D,eq);
    %keyboard
    %%
    
    %Counters to see how many nodes are on each boundary
    nAB = zeros(1,nen);
    nBC = zeros(1,nen);
    nCD = zeros(1,nen);
    nDA = zeros(1,nen);
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
        
        
        if (any(node==NodesBC)==1)
            nBC(1)=nBC(1)+1;
            if(nBC(2)==0)
                nBC(2)=node;
            else
                nBC(3)=node;
            end
        end
        
        if (any(node==NodesCD)==1)
            nCD(1)=nCD(1)+1;
            if(nCD(2)==0)
                nCD(2)=node;
            else
                nCD(3)=node;
            end
        end
        
        if (any(node==NodesDA)==1)
            nDA(1)=nDA(1)+1;
            if(nDA(2)==0)
                nDA(2)=node;
            else
                nDA(3)=node;
            end
        end
        %keyboard
        
    end
    %end for the boundary loop
    
    
    %%
    %keyboard
    Ke_additionAB=0;
    Ke_additionBC=0;
    Ke_additionCD=0;
    Ke_additionDA=0;
    fe_additionAB=0;
    fe_additionBC=0;
    fe_additionCD=0;
    fe_additionDA=0;
    
    %Now that we know how many nodes are on each boundary we can send them
    %into the boundary conditions calculator
    if(nAB(1)==2)
        %keyboard
        [Ke_additionAB,fe_additionAB]=elem_cond(ex(e,:),ey(e,:),ep,alpha1,T1_inf,find(elem(e).cn==setdiff(elem(e).cn,nAB(2:end))));
        %keyboard
    end
    
    if(nBC(1)==2)
        [Ke_additionBC,fe_additionBC]=elem_cond(ex(e,:),ey(e,:),ep,alpha2,T2_inf,find(elem(e).cn==setdiff(elem(e).cn,nBC(2:end))));
        %keyboard
    end
    
    if(nCD(1)==2)
        [Ke_additionCD,fe_additionCD]=elem_cond(ex(e,:),ey(e,:),ep,alpha2,T2_inf,find(elem(e).cn==setdiff(elem(e).cn,nCD(2:end))));
        %keyboard
    end
    
    if(nDA(1)==2)
        [Ke_additionDA,fe_additionDA]=elem_cond(ex(e,:),ey(e,:),ep,alpha2,T2_inf,find(elem(e).cn==setdiff(elem(e).cn,nDA(2:end))));
        %keyboard
    end
    
    clear n nAB nBC nCD nDA
    
    Ke_addition=0;
    fe_addition=0;
    
    Ke_addition=Ke_additionAB+Ke_additionBC+Ke_additionCD+Ke_additionDA;
    fe_addition=fe_additionAB+fe_additionBC+fe_additionCD+fe_additionDA;
    
    clear Ke_additionAB Ke_additionBC Ke_additionCD Ke_additionDA fe_additionAB fe_additionBC fe_additionCD fe_additionDA
    %keyboard
    
    Ke=Ke+Ke_addition;
    fe=fe+fe_addition;
    %keyboard
    
    % call assem to assamble the element stiffness matrix to the global one
    [K,f]=assem(edof(e,:),K,Ke,f,fe);
    C=assem(edof(e,:),C,Ce);
    %keyboard
end
%End of element loop
clear Ce Ke Ke_addition fe fe_addition k c t ep eq

%%
%Set Dirichelet BC and initial conditions
bc = [nodedof(NodesTubes(:),2),ones(size(NodesTubes,2),1)*Tg];
a=ones(ndof,1)*T_0;

clear NodesAB NodesBC NodesCD NodesDA NodesTubes

% simulation time and time/load step
stopTime= 8*60*60;

%dt is timestep in seconds
dt = 100;
time = 0;
step = 0;

% plot the temperatures of the elements of the steady state

% [a,Q]=solve(K,f,bc);
% ed=extract(edof,a);
% fill(ex',ey',ed')
% colormap('jet')
% h=colorbar;
% h.Label.String='Temperature in centigrade';
% title('Temperature distribution for the steady state solution')
% caxis([-50 500]);
% keyboard

K=K+C/dt;

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