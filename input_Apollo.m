%input_Apollo_CALFEM(ndf)

x = dlmread('Apollo.nodes'); % 1st colum containes node number
x = x(:,2:end);              % remove 1st column

x(:,1) = x(:,1)+0.75;
x(:,2) = x(:,2)-3;


% number of node points and number of spatial dimensions
[nnp, ndm] = size(x);

% total number of dofs
ndof = nnp*ndf;

conn = dlmread('Apollo.conn'); % 1st colum containes element number
conn = conn(:,2:end);          % remove 1st column

[nel,nen] = size(conn);

% allocate memory for the element structure 
% 1. organize connectivity list as a struct
% 2. determine edof, element degrees of freedom
% 3. reserve memory for state variables
%    gradT  [gradT_x    gradT_y] and 
%    q      [    q_x        q_y]
% number of state variables
nsv = 1;

elem = repmat(struct('cn',zeros(nen,1), ...
                     'edof', zeros(nen*ndf,1), ...
                     'stateVar0', zeros(1,nsv), ...
                     'stateVar', zeros(1,nsv) ...
                    ),  nel, 1 );

% for CALFEM all element dofs in one matrix
edof = zeros(nel, nen*ndf+1);
% for CALFEM all element x- and y- coord in ex and ey matrix
ex = zeros(nel, nen);
ey = zeros(nel, nen);

%keyboard 
nodedof = zeros(nnp,ndf+1);

for n=1:nnp
    nodedof(n,1) = n;
    nodedof(n,2:end) = [n*ndf-(ndf-1):n*ndf];
end

clear n
%keyboard

for e = 1:nel
    
    elem(e).cn = conn(e,:);
    % dofs belonging to the nodes of element e
    elem(e).edof = elem(e).cn; 
    
    % CALFEM
    edof(e,:) =  [e nodedof(elem(e).cn(1),2:end) nodedof(elem(e).cn(2),2:end) nodedof(elem(e).cn(3),2:end)];
    % CALFEM
    ex(e,:) = x(elem(e).cn, 1);
    ey(e,:) = x(elem(e).cn, 2);
end

clear e ndf ndm nsv

