function rows_used = Plot_Yeo_reference(L,cortex_dir,scortex)
%
% Render the 7 resting-state networks of Yeo et al. (2011) in the top right
% corner of a pyramid figure (empty part of the pyramid), in the same view
% as the FC states, with the name of each network below, as a reference
% for the colors of the FC states.
%
% INPUT:
% L           layout from Pyramid_layout
% cortex_dir  view of the brains ('TopView' or 'SideView')
% scortex     smoothed brain mask used for the transparent cortex
%
% OUTPUT:
% rows_used   height used by the reference, in rows of the pyramid (to
%             place the legend below it)
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

[YeoColor,Yeo_names] = Yeo_colors;
S = load('ParcelsMNI2mm.mat','V_Yeo7');
N_Yeo = 7;

% Area of the reference: right half of the figure, top rows
x0 = 0.55;
width = (1 - L.right - x0)/N_Yeo;
brain_height = 1.4*L.h;
name_height = 0.4*L.h;
y_brain = 1 - L.top - brain_height;

for net = 1:N_Yeo
    axes('Position',[x0+(net-1)*width y_brain width brain_height]);
    hold on
    % Transparent cortex, as in the pyramid
    cortexpatch = patch(isosurface(scortex,0.1),'FaceColor',[0.9 0.9 0.9],'EdgeColor','none','FaceAlpha',0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    % Network
    patch(isosurface(smooth3(S.V_Yeo7 == net),0.3),'FaceColor',YeoColor(net,:),'EdgeColor','none');
    material dull
    lighting gouraud
    switch cortex_dir
        case 'TopView'
            view(-90,90)
        case 'SideView'
            view(0,0)
    end
    daspect([1 1 1])
    camlight;
    xlim([5 105])
    ylim([5 85])
    zlim([10 80])
    axis off

    % Name of the network below
    annotation('textbox',[x0+(net-1)*width y_brain-name_height width name_height],'String',Yeo_names{net},...
               'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','top',...
               'FontSize',9,'FontWeight','bold','Color',YeoColor(net,:)*0.8,'FitBoxToText','off');
end

rows_used = (brain_height + name_height)/L.h;
