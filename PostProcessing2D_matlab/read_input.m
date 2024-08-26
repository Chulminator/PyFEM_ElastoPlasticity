clear all; close all;
% title = "./mesh/SENT_0"
title = "../Result/PlateWithHole_113_93"
flag = 'disp'
ReadMesh
%%
figure
set(gcf, 'renderer', 'painters');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
%%
patch('Faces',A(:,2:end),'Vertices',Node(:,1:2),'FaceVertexCData',Node(:,5),'FaceColor','interp'); 
% p.EdgeColor='none';   % no edge
hold on

%%
for i = 1:length(Cohesive)

%     A(Cohesive(i,1))
%     A(Cohesive(i,2))

% 
%     tmp = Cohesive(i,3)+2+1;
%     if isnan(A(Cohesive(i,1),tmp))
%         tmp = 2;
%     end
%     x = [Node(A(Cohesive(i,1),Cohesive(i,3)+1+1),1), Node(A(Cohesive(i,1),tmp),1)];
%     y = [Node(A(Cohesive(i,1),Cohesive(i,3)+1+1),2), Node(A(Cohesive(i,1),tmp),2)];
%     plot(x,y,'r',LineWidth=1.5 )
% 
    tmp = Cohesive(i,4)+2+1;
    if isnan(A(Cohesive(i,2),tmp))
        tmp = 2;
    end
    x = [Node(A(Cohesive(i,2),Cohesive(i,4)+1+1),1), Node(A(Cohesive(i,2),tmp),1)];
    y = [Node(A(Cohesive(i,2),Cohesive(i,4)+1+1),2), Node(A(Cohesive(i,2),tmp),2)];
    plot(x,y,'g',LineWidth=3 )
end
%%
% 
%     x = [Node(A(Cohesive(i,2),3),1), Node(A(Cohesive(i,4),tmp),1)];
%     y = [Node(A(Cohesive(i,2),3),2), Node(A(Cohesive(i,4),tmp),2)];
%     plot(x,y,'k',LineWidth=5 )
% 
%     x = [Node(A(Cohesive(i,1),4),1), Node(A(Cohesive(i,1),2),1)];
%     y = [Node(A(Cohesive(i,1),4),2), Node(A(Cohesive(i,1),2),2)];
%     plot(x,y,'b',LineWidth=5 )

% pbaspect([max(Node(:,1))-min(Node(:,1)),max(Node(:,2))-min(Node(:,2)),1])
% NNode = size(Node,1);
% % for i = 1: NNode
% for j = 1: 3
%     j = A(k,j+1)
% %     plot(Node(i,1), Node(i,2), 'bo')
%     text(Node(j,1), Node(j,2),num2str(j),FontSize=15)
% end
% tmp = [ 704 738 8651 8617];

% p1=patch('Faces', coh_Elem, 'Vertices', Node(:,1:2),'FaceColor','w') ;
% p1=patch('Faces', tmp, 'Vertices', Node(:,1:2),'FaceColor','w') ;
% p1.EdgeColor='black';
% p.FaceColor='flat';

% axis off

colormap(jet)
colorbar
% caxis([-6e+08, 6e+08])

% plot(Node(:,1),Node(:,2),'ro','markerfacecolor','r');

% hold on
% contourf(Node(:,1), Node(:,2), Node(:,3));
% contourf(Node(:,1), Node(:,2), Node(:,1).* Node(:,2)), 'interp',... 
        

