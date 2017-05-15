
function [fe_addition ]=elem_cond(ex, ey, t, P,node_not_on_boundary)
%
% [fe_addition]=elem_cond(ex,ey,ep,alpha,T_inf,node_not_on_boundary)
%-------------------------------------------------------------
% PURPOSE
%  Compute the addition to the element load vector matrix for a
%  triangular field element due to secon order polynomial Neumann BC's.
%
% INPUT:  ex = [x1 x2 x3]
%         ey = [y1 y2 y3]               element coordinates
%
%         t                             element thickness
%
%         P                             Polynomial coefficients
%
%         node_not_on_boundary = [n]    number of the node not on the
%                                       boundary
%
% OUTPUT: fe_addition : addition to the element load vector (3 x 1) due to
%         the Neumann BC
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

xmid = L/2;
xmin=min(ex_temp(1),ex_temp(2));

%Esitmate the polynomial by the vvalue att he midpoint
p=polyval(P, xmin+xmid);

switch node_not_on_boundary
    case 1
        fe_addition=t*p*L/2*[0;0;0;1;0;1];
    case 2  
        fe_addition=t*p*L/2*[0;1;0;0;0;1];
    case 3
        fe_addition=t*p*L/2*[0;1;0;1;0;0];
    otherwise
        %This should never happen
        'error'
        fe_addition=0;
end
%end of switch

end


