function PlotMesh(mesh,option);
    figure
    set(gcf, 'renderer', 'painters');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [4 2]);
    hold on
    %%
    if option == 0
        patch('Faces', mesh.Element(:,2:end),'Vertices',mesh.Node(:,2:3),'FaceColor','none'); 
    elseif option == 1
        patch('Faces', mesh.Element(:,2:end),'Vertices',mesh.Node(:,2:3)+mesh.Node(:,4:5),'FaceColor','none'); 
    elseif option == 2
        patch('Faces', mesh.Element(:,2:end),'Vertices',mesh.Node(:,2:3),'FaceVertexCData',mesh.Node(:,4),'FaceColor','interp'); 
%         patch('Faces', mesh.Element(:,2:end),'Vertices',mesh.Node(:,2:3),'FaceVertexCData',mesh.Node(:,6),'FaceColor','interp'); 
    elseif option == 4
        patch('Faces', mesh.Element(:,2:end),'Vertices',mesh.Node(:,2:3),'FaceVertexCData',mesh.Node(:,6),'FaceColor','interp'); 
%         patch('Faces', mesh.Element(:,2:end),'Vertices',mesh.Node(:,2:3),'FaceVertexCData',mesh.Node(:,6),'FaceColor','interp'); 
    else
        assert(false,'Invalid option keys')
    end

%     Edgecolor
%     p.EdgeColor='none';   % no edge

    hold off
end