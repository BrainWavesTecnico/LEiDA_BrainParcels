function Plot_K_boxplot_FO(data_dir,stats_file,save_dir,selectedK,parcellation)
%
% Plot the fractional occupancy (FO) of each condition for each of the FC
% states of the specified K: boxplots with one point per subject (or scan,
% as in the statistical tests), colored with the Yeo network of the FC
% state, and bars for the pairs of conditions colored by significance
% (dark red p<0.05, dark green p<0.05/K, dark blue p<0.05/sum(rangeK)),
% labelled with the p-value and effect size if K <= 6 (see Plot_FO_boxplot).
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% stats_file    full path of LEiDA_Stats_FracOccup.mat of the study design
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% parcellation  parcellation atlas used to segment the brain
%
% OUTPUT:
% .fig/.png     Boxplots of the fractional occupancy values for each
%               condition across FC states
%
% Authors: Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt
%          Joana Cabral, University of Minho, joanacabral@med.uminho.pt

load([data_dir 'LEiDA_Clusters.mat'], 'Kmeans_results', 'rangeK');
k = find(rangeK == selectedK);
n_areas = size(Kmeans_results{k}.C,2);

% Yeo network assigned to each FC state
YeoColor = Yeo_colors;
[dice_V_yeo7, overlap_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);
clear Kmeans_results

% Show the p-values and effect sizes only if there is room for them
compact = selectedK > 6;

disp(' ');
disp(['Plotting the fractional occupancy of each of the ' num2str(selectedK) ' FC states by condition:'])
Fig = figure('Position', get(0, 'Screensize'));
for c = 1:selectedK
    if selectedK <= 10
        subplot(1,selectedK,c)
    else
        subplot(2,10,c)
    end
    [~, net] = max(dice_V_yeo7(k,c,:) + 2*overlap_yeo7(k,c,:));
    if overlap_yeo7(k,c,net)
        color = YeoColor(net,:);
    else
        color = [0 0 0];
    end
    Plot_FO_boxplot(stats_file,k,c,rangeK,color,compact);
    xtickangle(30)
    title(['FC state ' num2str(c)],'FontSize',11);
    set(gca,'color','none')
end
Save_png(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_BoxplotFracOccup.png']));
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_BoxplotFracOccup.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_BoxplotFracOccup']);
