function [k, m] = coeff_2d( x1, y1, x2, y2)
%returns k and m for straight line
%   out: k - line coefficient
%        m - constant   

%   in:  x1,y1,x2,y2 - point coordinates   

k=(y2-y1)/(x2-x1);

m=y1-k*x1;

end

