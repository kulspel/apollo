function [Ke_addition, fe_addition]=robin_heat_tria2(ex, ey,ep, alpha, T_inf, node_not_on_boundary)
% 
% [Ke_addition,fe_addition]=robin_heat_tria2(ex,ey,ep,alpha,T_inf,node_not_on_boundary)
%-------------------------------------------------------------
% PURPOSE
%  Compute the addition to the element stiffness (conductivity) matrix for a 
%  triangular field element due to Robin BC's.
%
% INPUT:  ex = [x1 x2 x3]
%         ey = [y1 y2 y3]               element coordinates
%
%         ep = [t]                      element thickness  	 
%                             
%         alpha = [alpha]               convection coefficient
%
%         T_inf = [T_inf]               Far field temperature
%
%         node_not_on_boundary = [n]    number of the node not on the
%                                       boundary
%
% OUTPUT: Ke_addition :  addition to the element 'stiffness' matrix (3 x3) 
%                        due to robin BC
%
%         fe_addition : addition to the element load vector (3 x 1) due to
%                       the Robin BC
%-------------------------------------------------------------

if(node_not_on_boundary>3 || node_not_on_boundary<1)
'error:node_not_on_boundary invalid value, must be between 1 and 3 was'
node_not_on_boundary
Ke_addition=0;
fe_addition=0;
return
end
%keyboard
%Remove the node thats not on the boundary to easier get the coordinates
%for the boundary nodes
ex_temp= ex;
ex_temp(node_not_on_boundary)=[];

ey_temp= ey;
ey_temp(node_not_on_boundary)=[];

x1=ex_temp(1);
x2=ex_temp(2);

y1=ey_temp(1);
y2=ey_temp(2);

xdiff=x2-x1;
ydiff=y2-y1;

L=sqrt(xdiff^2+ydiff^2);

switch node_not_on_boundary
    case 1
        Ke_addition=alpha*L*ep/6*[0 0 0;0 2 1;0 1 2];
        %check sign on fe_addition
        fe_addition=alpha*T_inf*L/2*[0;1;1];
        
    case 2
        Ke_addition=alpha*L*ep/6*[2 0 1;0 0 0;1 0 2];
        fe_addition=alpha*T_inf*L/2*[1;0;1];
    case 3
        Ke_addition=alpha*L*ep/6*[2 1 0;1 2 0;0 0 0];
        fe_addition=alpha*T_inf*L/2*[1;1;0];
    otherwise
        %This should never happen
        'error'
        Ke_addition=0;
        fe_addition=0;
end



end