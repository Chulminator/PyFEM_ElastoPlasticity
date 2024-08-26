clear all; close all;
% % % % % % % % % % % Input % % % % % % % % % % % 
% title = "./mesh/SENT_0"
% title = "../Result/PlateWithHole_113_93/PlateWithHole_113_93";
% title = "../Result/PlateWithHole_113_93_plastic/PlateWithHole_113_93_Stress";
title = "../Result/PlateWithHole_113_93_plastic/PlateWithHole_113_93_CustomizedAttribute";
% title = "../Result/PlateWithHole_113_93/PlateWithHole_113_93_CustomizedAttribute";

% option - None;
% option = 0;
% 'disp';
% option = 1;
% 'stress';
% option = 2;
% 'dispstress';
% option = 3;
% 'custom';
option = 4;
% % % % % % % % % % % Input % % % % % % % % % % % 

str = sprintf('title\t: %s\n',title);
disp(str)
%%
mesh = ReadMesh(title);
PlotMesh(mesh, option);
colorbar

% %  Node Num
% for i =1:size(mesh.Node,1)
%     text(mesh.Node(i,2),mesh.Node(i,3),num2str(i))
% end

% % Element Num
% for i =1:size(mesh.Node,1)
%     text(sum(mesh.Node(mesh.Element(i,2:5),2))/4,sum(mesh.Node(mesh.Element(i,2:5),3))/4,num2str(i))
% end

plot_size