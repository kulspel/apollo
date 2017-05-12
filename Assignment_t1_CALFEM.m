clear all
close all
clc

%%
% number of degrees of freedom per node
ndf = 1;  % temperature is a scalar!

% obtain coordinates and element properties
[ndof, nel, ex, ey, edof, elem, nodedof]=input_Apollo_CALFEM(ndf);
%[ndof, nel, ex, ey, edof, elem, nodedof]=input_2tria3_CALFEM(ndf);

nen = size(elem(1).cn,2);
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

%mass specific heat capacity of titanium 
%520 J/(kg·K)
c=520;

%%
%Other useful parameters

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

%Boundaries
%[NodesAB,NodesBC,NodesCD,NodesDA,NodesTubes] = nodes_2tria3();
%NodesAB=[3 4]; NodesBC=[2 3]; NodesCD=[]; NodesDA=[3 4]; NodesTubes=[1 2];
[NodesAB,NodesBC,NodesCD,NodesDA,NodesTubes] = BCNodes();

%keyboard

% simulation time and time/load step
stopTime = 1;
dt = 1;

% loop over all time/load steps
time = 0;
step = 0;
while (time < stopTime)
    time = time + dt;
    step = step + 1;

    % initialise global heat capacity matrix, tangent stiffness matrix K and global force vector
    C=zeros(ndof);
    K=zeros(ndof);
    f=zeros(ndof,1);
    
    %If it's the first step then set initial conditions
    if(step==1)
        a=ones(ndof,1)*T_0;
    end
    
    %%
    %Element loop
    for e = 1:nel
        % call plantml to calculate the element heat capacity matrix for element e
        Ce=plantml(ex(e,:),ey(e,:),c*rho);
        
        %element temperatures
        ae=a(edof(e,2:end));
        
        % call flw2te to calculate the element stiffness matrix for element e
        [Ke,fe] =flw2te(ex(e,:),ey(e,:),ep,D,eq);
        %keyboard
    
        %Counters to see how many nodes are on each boundary
        nAB = zeros(1,nen);
        nBC = zeros(1,nen);
        nCD = zeros(1,nen);
        nDA = zeros(1,nen);
        %nTubes = 0;
    
        %%
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
            [Ke_additionAB,fe_additionAB]=robin_heat_tria2(ex(e,:),ey(e,:),ep,alpha1,T1_inf,find(elem(e).cn==setdiff(elem(e).cn,nAB(2:end))));
            %keyboard
        end
    
        if(nBC(1)==2)
            [Ke_additionBC,fe_additionBC]=robin_heat_tria2(ex(e,:),ey(e,:),ep,alpha2,T2_inf,find(elem(e).cn==setdiff(elem(e).cn,nBC(2:end))));
            %keyboard
        end
    
        if(nCD(1)==2)
            [Ke_additionCD,fe_additionCD]=robin_heat_tria2(ex(e,:),ey(e,:),ep,alpha2,T2_inf,find(elem(e).cn==setdiff(elem(e).cn,nCD(2:end))));
            %keyboard
        end
    
        if(nDA(1)==2)
            [Ke_additionDA,fe_additionDA]=robin_heat_tria2(ex(e,:),ey(e,:),ep,alpha2,T2_inf,find(elem(e).cn==setdiff(elem(e).cn,nDA(2:end))));
            %keyboard
        end
    
        Ke_addition=0;
        fe_addition=0;
    
        Ke_addition=Ke_additionAB+Ke_additionBC+Ke_additionCD+Ke_additionDA;
        fe_addition=fe_additionAB+fe_additionBC+fe_additionCD+fe_additionDA;
    
        clear Ke_additionAB Ke_additionBC Ke_additionCD Ke_additionDA fe_additionAB fe_additionBC fe_additionCD fe_additionDA
        %keyboard
    
        Ke=Ke-Ke_addition;
        fe=fe-fe_addition;
        %keyboard
        
        % call assem to assamble the element stiffness matrix to the global one
        [K,f]=assem(edof(e,:),K,Ke,f,fe);
        C=assem(edof(e,:),C,Ce);
        %keyboard
    end
    %End of element loop
    
    %%
    bc = [nodedof(NodesTubes(:),2),ones(size(NodesTubes,2),1)*Tg]; 

    %keyboard

    clear nAB nBC nCD nDA Ce Ke Ke_addition fe fe_addition k
    %clear NodesAB NodesBC NodesCD NodesDA NodesTubes
    
    % call solveq
    [a,Q]=solve(K+C/dt,f+(C/dt)*a,bc);
    %[a,Q]=solve(K,f,bc);
    %keyboard

    % call extract 
    ed=extract(edof,a);
    %keyboard

    %fill(ex',ey',ed')
    fill(ex',ey',ed')

    % plot the temperatures of the elements
    colormap('jet')
    fprintf(1,'step= %2d time= %8.4e dt= %8.4e\n', step, time, dt);
    %keyboard
end
%End of time loop
%d0=ones(ndof,1)*T_0;
%T=5;
%dt=0.1;
%ti=[0:dt:T];

%ip=[dt T 1 [length(ti),ndof,ti,edof']]
%[Tsnap,D,V]=step1(K,C,d0,ip,f,bc)
