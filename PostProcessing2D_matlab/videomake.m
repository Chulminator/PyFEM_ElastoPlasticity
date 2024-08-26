% videomake
% clear all; close all;
% title = "./mesh/SENT_456000_stress"
for i = 1:80+10
    title =sprintf('./mesh/SENT_%s_stress',num2str((i-1)*10000));
    read_input
end