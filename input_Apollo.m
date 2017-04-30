function [x,elem,nel]=input_Apollo()

x = dlmread('Apollo.nodes'); % 1st colum containes node number
x = x(:,2:end);              % remove 1st column

x(:,1) = x(:,1)+0.75;
x(:,2) = x(:,2)-3;

conn = dlmread('Apollo.conn'); % 1st colum containes element number
conn = conn(:,2:end);          % remove 1st column
[nel,nen] = size(conn);

elem=repmat(struct('cn',zeros(1,nen)),nel,1);


for e = 1:nel
    elem(e).cn = conn(e,:);
end

end