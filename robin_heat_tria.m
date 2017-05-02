function [Ke_addition, fe_addition]=robin_heat_tria(ex, ey, ep, alpha,T_inf, node_not_on_boundary)

if(node_not_on_boundary>3 || node_not_on_boundary<1)
'error:node_not_on_boundary invalid value, must be between 1 and 3 was'
node_not_on_boundary
Ke_addition=0;
fe_addition=0;
return
end
%format long

C=[1 ex(1) ey(1);1 ex(2) ey(2);1 ex(3) ey(3)];

C_inv=inv(C);

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

if xdiff>0
    k=ydiff/xdiff;
    
    m=y1-k*x1;
    
    %keyboard
    
    NT_int=([ x2;
             x2^2/2;
             k*x2^2/2+m*x2]-...
           [ x1;
             x1^2/2;
             k*x1^2/2+m*x1])*sqrt(1+k^2);
                
    N_kvadrat_int=([x2               x2^2/2              k*x2^2/2+m*x2;
                   x2^2/2           x2^3/3              k*x2^3/3+m*x2^2/2;
                   k*x2^2/2+m*x2    k*x2^3/3+m*x2^2/2   k^2*x2^3/3+m*k*x2^2+m^2*x2]-...
                  [x1               x1^2/2              k*x1^2/2+m*x1;
                   x1^2/2           x1^3/3              k*x1^3/3+m*x1^2/2;
                   k*x1^2/2+m*x1    k*x1^3/3+m*x1^2/2   k^2*x1^3/3+m*k*x1^2+m^2*x1])*sqrt(1+k^2);
             
    
else
    k=xdiff/ydiff;
    
    m=x1-k*y1;
    
    NT_int=([ y2;
             k*y2^2/2+m*y2;
             y2^2/2]-...
           [ y1;
             k*y1^2/2+m*y1;
             y1^2/2])*sqrt(1+k^2);
    
    N_kvadrat_int=([y2               k*y2^2/2+m*y2                       y2^2/2;
                   k*y2^2/2+m*y2    k^2*y2^3/3+2*m*k*y2^2/2+m^2*y2      k*y2^3/3+m*y2^2/2;
                   y2^2/2           k*y2^3/3+m*y2^2/2                   y2^3/3]-...
                  [y1               k*y1^2/2+m*y1                       y1^2/2;
                   k*y1^2/2+m*y1    k^2*y1^3/3+2*m*k*y1^2/2+m^2*y1      k*y1^3/3+m*y1^2/2;
                   y1^2/2           k*y1^3/3+m*y1^2/2                   y1^3/3])*sqrt(1+k^2);
end


%C_inv'*N_kvadrat_int*C_inv*6/L
%C_inv'*NT_int/L

Ke_addition=C_inv'*N_kvadrat_int*C_inv*alpha*ep;
%check sign on fe_addition
fe_addition=C_inv'*NT_int*alpha*T_inf*ep;

%integrating between two nodes
end