function Plot_K_overlap_yeo_nets(data_dir,save_dir,selectedK,parcellation)
%
% Compute the overlap between each PL state (1 to selectedK) and the RSNs
% defined in Yeo et al., (2011). Plot the overlap between each PL state and
% each of the Yeo RSNs.
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% parcellation  parcellation template used to segment the brain
%
% OUTPUT:
% .fig/.png     Plot of the overlap of each PL state with each of the RSNs
%               defined in Yeo et al., (2011)
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';

% Load required data:
if isfile([data_dir file_cluster])
    load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
end

% Matrix of dimension selectedK*90 where each row represents one FC state
V = Kmeans_results{rangeK == selectedK}.C;

% Number of areas considered for analysis
n_areas = size(V,2);

% Colors of the 7 Yeo et al. (2011) networks (defined in Yeo_colors.m)
YeoColor = Yeo_colors;

% Compute the overlap (weighted Dice coefficient) between cluster centroids
% obtained for each value of K and the resting-state networks defined by
% Yeo et al. 2011, and the largest Dice reachable with this parcellation
[dice_V_yeo7, overlap_yeo7, dice_max] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);
clear Kmeans_results

disp(' ');
disp(['Plotting overlap with Yeo RSNs for K = ' num2str(selectedK) ' clusters.'])
disp(' ')
Fig = figure('Position', get(0, 'Screensize'));

b = bar(1:selectedK,squeeze(dice_V_yeo7(rangeK == selectedK,1:selectedK,1:7)));
for i = 1:7
    set(b(i),'FaceColor',YeoColor(i,:));
end
xlabel(['K = ' num2str(selectedK) ' clusters'], 'FontSize',14);
ylabel({'Weighted Dice coefficient with 7 RSNs';'(Yeo et al. 2011)'}, 'FontSize', 14);
% Legend with the largest Dice reachable for each network with this parcellation
[~,Yeo_names] = Yeo_colors;
legend(strcat(Yeo_names,arrayfun(@(d) sprintf(' (max %.2f)',d),dice_max,'UniformOutput',false)), ...
       'Location','northwest','FontSize',14);
legend boxoff;
box off;
% Mark with * the network assigned to each FC state (largest Dice
% coefficient, large enough relative to the max, see Overlap_LEiDA_Yeo)
for centroidK = 1:selectedK
    NetYeo = find(overlap_yeo7(rangeK == selectedK,centroidK,:));
    if ~isempty(NetYeo)
        text(b(NetYeo).XEndPoints(centroidK),b(NetYeo).YEndPoints(centroidK),'*','fontsize',25,...
             'Color',[0 0 0],'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end

Save_png(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_OverlapYeoNets.png']));
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_OverlapYeoNets.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_OverlapYeoNets']);
