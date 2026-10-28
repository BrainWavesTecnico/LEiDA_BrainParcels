function Plot_K_repertoire(data_dir,stats_file,save_dir,selectedK,parcellation)
%
% Plot the centroids in 3D glass brain, in vector format and the boxplots
% of the fractional occupancy values
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% stats_file    full path of LEiDA_Stats_FracOccup.mat of the study design
%               ('' to skip the fractional occupancy)
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% parcellation  parcellation template used to segment the brain
%
% OUTPUT:
% .fig/.png     Plot of summary information for the repertoire of K
%               centroids
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';
% Load required data:
load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
plot_FO = ~isempty(stats_file) && isfile(stats_file);

% Matrix of dimension selectedK*90 where each row represents one FC state
V = Kmeans_results{rangeK == selectedK}.C;
% Scale each cluster centroid by its maximum value and transpose the matrix
V = V'./max(abs(V'));

% Number of areas considered for analysis
n_areas = size(V,1);

% Colors of the 7 Yeo et al. (2011) networks (defined in Yeo_colors.m)
[YeoColor,Yeo_names] = Yeo_colors;

% Compute the overlap between cluster centroids obtained for each value of
% K and the resting-state networks defined by Yeo et al. 2011
[dice_V_yeo7, overlap_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);
clear Kmeans_results

% LEiDA networks colored according to closest RSN
Volume = struct2array(load('ParcelsMNI2mm',['V_' parcellation]));
% Only the areas considered
Volume(Volume > n_areas) = 0;
Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

% Order of the areas (left then right hemisphere in AAL)
Order = Areas_order(parcellation,n_areas);

disp(' ');
disp(['Plotting summary information for the ' num2str(selectedK) ' PL states:'])
Fig = figure('Position', get(0, 'Screensize')); 
colormap(jet)
for c = 1:selectedK
    
    [~, net] = max(dice_V_yeo7(rangeK == selectedK,c,:) + 2*overlap_yeo7(rangeK == selectedK,c,:));
    
    % First Plot view from top
    subplot_tight(8,selectedK,c+selectedK,0.02)
    hold on
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    
    % Color all areas with positive V with transparency proportional to contribution
    n_pos = find(V(:,c) > 0);
    if numel(n_pos) > 0
        for region = 1:length(n_pos)
            region_Vol = Volume == n_pos(region); % 3D volume with 1 only in the voxels of that area
            sregion = smooth3(region_Vol > 0);
            
            if overlap_yeo7(rangeK == selectedK,c,net)
                patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            else
                patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            end
        end
    else % Plot of the global mode with "gray/transparent" color
        region_Vol = Volume > 0; % 3D volume with 1 only in the voxels of that area
        sregion = smooth3(region_Vol > 0);
        patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
    end
    
    if overlap_yeo7(rangeK == selectedK,c,net)
        title({['FC state ' num2str(c), newline, Yeo_names{net}]} ) %,Yeo_names{c}})
    elseif numel(n_pos) == 0
        title({['FC state ' num2str(c), newline, 'Global Mode']} ) %,Yeo_names{c}})
    else
        title({['FC state ' num2str(c), newline, '']} ) %,Yeo_names{c}})
    end
    
    material dull; lighting gouraud
    daspect([1 1 1])
    view(-90,90)    % Top view    Side view:   view(0,0)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    %  Same but view from the side
    subplot_tight(8,selectedK,c+2*selectedK,0.025)
    hold on
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    
    n_pos = find(V(:,c) > 0);
    if n_pos
        for region = 1:length(n_pos)
            region_Vol = Volume == n_pos(region); % 3D volume with 1 only in the voxels of that area
            sregion = smooth3(region_Vol > 0);
            
            if overlap_yeo7(rangeK == selectedK,c,net)
                patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            else
                patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
            end
        end
    else
        region_Vol = Volume > 0; % 3D volume with 1 only in the voxels of that area
        sregion = smooth3(region_Vol > 0);
        patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
    end
    
    material dull; lighting gouraud
    daspect([1 1 1])
    view(0,0)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    % PLOT THE VECTOR OF EACH FC STATE
    subplot_tight(8,selectedK,[c+3*selectedK c+4*selectedK],0.045)
    Vo = V(Order,c);
    hold on
    barh(Vo.*(Vo<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(Vo.*(Vo>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 n_areas+1])
    xlim([-1 1])
    set(gca,'XTick',[-1,0,1],'Fontsize',8)
    set(gca,'YTick',[10:20:n_areas],'Fontsize',8)
    ax = gca;
    ax.XAxis.FontSize = 8;
    grid on
    
    % FRACTIONAL OCCUPANCY IN EACH CONDITION (one point per subject or
    % scan, as in the tests; bars for the pairs of conditions with p<0.05,
    % colored by significance)
    if plot_FO
        subplot_tight(8,selectedK,[c+5*selectedK c+6*selectedK],0.035)
        if overlap_yeo7(rangeK == selectedK,c,net)
            Plot_FO_boxplot(stats_file,find(rangeK == selectedK),c,rangeK,YeoColor(net,:),true);
        else
            Plot_FO_boxplot(stats_file,find(rangeK == selectedK),c,rangeK,[0 0 0],true);
        end
        xtickangle(30)
        if c == 1
            ylabel('Fractional Occupancy');
        end
    end
end

Save_png(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_RepertoireSummary.png']));
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_RepertoireSummary.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_RepertoireSummary']);

