function [ex, ey, elem] = input_2tria3_CALFEM()

% global node coordinates
%       x      y
X = [ 0.0    0.0;
      2.0    0.0;
      2.0    1.0;
      0.0    1.0    ];

% number of node points and number of spatial dimensions
[nnp, ndm] = size(X);
 
% connectivity list
connectivity = [ 1       2       3; 
                 1       3       4];

[nel,nen] = size(connectivity);
             
ndf = 1; % temperature is a scalar quantity

% initialize and build up element struct elem
% containing the element connectivity list, cn,
% the element degree of freedom list, edof,
% and the element state variable list, stateVar.
elem = repmat(struct('cn',zeros(nen,1), ...
                     'edof',zeros(nen*ndf,1) ...
              ),nel,1);
% for CALFEM all element dofs in one matrix
edof = zeros(nel, nen+1);
% for CALFEM all element x- and y- coord in ex and ey matrix
ex = zeros(nel, nen);
ey = zeros(nel, nen);

for e = 1:nel
    
    elem(e).cn = connectivity(e,:);
    % dofs belonging to the nodes of element e
    elem(e).edof = elem(e).cn; 
    
    % CALFEM 
    edof(e,:) =  [e elem(e).cn];
    % CALFEM
    ex(e,:) = X(elem(e).cn, 1);
    ey(e,:) = X(elem(e).cn, 2);
end


% a similar JAVA class would be 
% public class Elem
% {
%     Vector<int> cn = new Vector<int>(nen);
%     Vector<int> edof = new Vector<int>(nen*ndf);
%     Vector<double> stateVar0 = new Vector<double>(nsv);
%     Vector<double> stateVar = new Vector<double>(nsv);
% }
%
% Vector<Elem> elem = new Vector<Elem>(nel);