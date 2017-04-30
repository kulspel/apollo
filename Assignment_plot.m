function  Assignment_plot(elem, x)

% nel number of elements in the mesh
nel = length(elem);

% nnp number of node points
% ndm number of spatial dimensions
[nnp, ndm] = size(x);

figure(1)
hold on
grid on
xlabel('coordinate x in cm')
ylabel('coordinate y in cm')

for e=1:nel
    % obtain element coordinates
    x_e = [x(elem(e).cn,1),x(elem(e).cn,2)];
    
    % calculate element center
    x_c = [(1/3*(x_e(1,1)+x_e(2,1)+x_e(3,1))),(1/3*(x_e(1,2)+x_e(2,2)+x_e(3,2)))];
    
    figure(1)
    % plot triangular mesh
    patch(x_e(:,1),x_e(:,2),'white','FaceAlpha',0.5)
    
    %Finds number of element nodes
    nen = size(elem(e).cn,2);
     
    % print node numbers in blue at node positions
     for i = 1:nen
     plot(x_e(i,1),x_e(i,2),'ob','MarkerSize',10,'MarkerFaceColor','w')
     text(x_e(i,1),x_e(i,2),num2str(elem(e).cn(i)),'Color','b','FontSize',6)
     end
    
    % print element numbers into element center in red
     plot(x_c(1),x_c(2),'sr','MarkerSize',10,'MarkerFaceColor','w')
     text(x_c(1),x_c(2),num2str(e),'Color','r','FontSize',6)
    
end
