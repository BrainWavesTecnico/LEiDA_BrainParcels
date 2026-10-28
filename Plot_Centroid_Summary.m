function Plot_Centroid_Summary(results_dir,stats_file,save_dir,selectedK,selectedC,parcellation,n_areas)
%
% Summary figure of one FC state (centroid):
%   - the FC state rendered in a transparent brain (side and top views),
%     colored with the Yeo network assigned to it (black if none)
%   - the strongest links between the areas with positive values in the
%     centroid (top view)
%   - the overlap with each of the 7 Yeo networks: weighted Dice, largest
%     Dice reachable with the parcellation, and fraction of the FC state
%     within each network (see Overlap_LEiDA_Yeo)
%   - the areas with positive values in the centroid, ordered by value
%   - the fractional occupancy of each condition (one point per subject or
%     scan, as in the tests), with the p-value, effect size (Hedges' g) and
%     its 95% bootstrap confidence interval for each pair of conditions
%
% INPUT:
% results_dir   directory with LEiDA_Clusters.mat
% stats_file    full path of LEiDA_Stats_FracOccup.mat ('' to skip the
%               fractional occupancy panel)
% save_dir      directory to save the figure
% selectedK     number of FC states K
% selectedC     FC state to plot (1 <= selectedC <= selectedK)
% parcellation  parcellation (atlas in ParcelsMNI2mm.mat)
% n_areas       number of areas ([] to use the size of the centroids)
%
% OUTPUT:
% .fig/.png     K<selectedK>C<selectedC>_Summary in save_dir
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

load(fullfile(results_dir,'LEiDA_Clusters.mat'),'Kmeans_results','rangeK');
k = find(rangeK == selectedK);
if isempty(k)
    error('LEiDA:centroid','K = %d is not in the range of K computed (%d to %d)',selectedK,rangeK(1),rangeK(end));
end
if selectedC < 1 || selectedC > selectedK || selectedC ~= round(selectedC)
    error('LEiDA:centroid','The FC state must be an integer between 1 and %d',selectedK);
end
if isempty(n_areas)
    n_areas = size(Kmeans_results{k}.C,2);
end
if isempty(Atlas_Nparcels(parcellation))
    error('LEiDA:centroid','Parcellation ''%s'' is not available in ParcelsMNI2mm.mat',parcellation);
end

% Centroid vector (sum of entries < 0, as in the clustering), scaled by
% its largest absolute value
Vc = Kmeans_results{k}.C(selectedC,1:n_areas);
if sum(Vc) > 0
    Vc = -Vc;
end
Vc = Vc/max(abs(Vc));

% Overlap with the Yeo networks and network assigned to the FC state
[dice_V_yeo7,overlap_yeo7,dice_max,frac_V_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);
dice_c = squeeze(dice_V_yeo7(k,selectedC,:))';
frac_c = squeeze(frac_V_yeo7(k,selectedC,:))';
net = find(overlap_yeo7(k,selectedC,:));
clear Kmeans_results
[YeoColor,Yeo_names] = Yeo_colors;
if isempty(net)
    state_color = [0 0 0];
    net_name = 'no network assigned';
else
    state_color = YeoColor(net,:);
    net_name = Yeo_names{net};
end

% Parcels, their centers and names, and brain mask
S = load('ParcelsMNI2mm.mat',['V_' parcellation]);
Volume = S.(['V_' parcellation]);
Volume(Volume > n_areas) = 0;
MNI_coord = zeros(n_areas,3);
for n = 1:n_areas
    [X,Y,Z] = ind2sub(size(Volume),find(Volume == n));
    MNI_coord(n,:) = mean([Y,X,Z],1);
end
labels = Parcellation_labels(parcellation,n_areas);
Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

% Voxels of the FC state (parcels with positive values)
sregion = smooth3(ismember(Volume,find(Vc > 0)));

disp(' ')
disp(['Plotting summary of FC state ' num2str(selectedC) ' (K = ' num2str(selectedK) '):'])
Fig = figure('Position',[50 50 1700 950],'Name',['K' num2str(selectedK) ' C' num2str(selectedC)]);
t = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
title(t,['FC state ' num2str(selectedC) ' of K = ' num2str(selectedK) ' (' parcellation ', ' ...
         num2str(n_areas) ' areas): ' net_name],'FontSize',16,'FontWeight','bold')

% FC state in the brain, side and top views
views = {[0 0],[-90 90]};
for v = 1:2
    nexttile(v)
    hold on
    cortexpatch = patch(isosurface(scortex,0.1),'FaceColor',[0.9 0.9 0.9],'EdgeColor','none','FaceAlpha',0.2);
    reducepatch(cortexpatch,0.05);
    isonormals(scortex,cortexpatch);
    if any(Vc > 0)
        patch(isosurface(sregion,0.3),'FaceColor',state_color,'EdgeColor','none');
    end
    material dull
    lighting gouraud
    view(views{v})
    daspect([1 1 1])
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([10 80])
    axis off
end

% Strongest links between areas with positive values (2% of all links)
nexttile(3)
Draw_links_in_cortex(Vc,MNI_coord,scortex,round(0.02*n_areas*(n_areas-1)/2));
view(-90,90)
camlight;

% Areas with positive values, ordered by value (at most 40)
nexttile(4,[2 1])
[Vo,order] = sort(Vc,'descend');
n_show = min(40,nnz(Vo > 0));
if n_show > 0
    barh(1:n_show,Vo(n_show:-1:1),'FaceColor',state_color,'EdgeColor','none','BarWidth',0.7)
    set(gca,'YTick',1:n_show,'YTickLabel',labels(order(n_show:-1:1)),'FontSize',8,'TickLabelInterpreter','none')
    ylim([0.4 n_show+0.6])
    xlim([0 1])
    xlabel('Centroid value (scaled)')
    title(['Areas in the FC state (' num2str(nnz(Vc > 0)) ' of ' num2str(n_areas) ')'],'FontSize',11)
else
    axis off
    text(0.5,0.5,{'All areas have the same sign:','global FC state (all areas in phase)'},...
         'HorizontalAlignment','center','FontSize',11)
end
box off

% Overlap with the 7 Yeo networks
nexttile(5)
b = bar(1:7,[dice_c; frac_c]','grouped','EdgeColor','none');
b(1).FaceColor = 'flat'; b(1).CData = YeoColor;
b(2).FaceColor = 'flat'; b(2).CData = YeoColor; b(2).FaceAlpha = 0.35;
hold on
plot(b(1).XEndPoints,dice_max,'_k','MarkerSize',10,'LineWidth',1.5)
hold off
set(gca,'XTick',1:7,'XTickLabel',Yeo_names,'FontSize',9)
xtickangle(45)
ylim([0 1])
ylabel('Overlap')
legend({'Weighted Dice','Fraction of the state','Max Dice (parcellation)'},'Location','northwest','FontSize',8)
legend boxoff
title(['Overlap with Yeo networks: ' net_name],'FontSize',11)
box off

% Fractional occupancy in each condition
nexttile(6,[1 2])
if ~isempty(stats_file) && isfile(stats_file)
    Plot_FO_boxplot(stats_file,k,selectedC,rangeK,state_color,false);
else
    axis off
    text(0.5,0.5,'No fractional occupancy statistics (run LEiDA\_Start with conditions)','HorizontalAlignment','center')
end

if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
file_name = ['K' num2str(selectedK) 'C' num2str(selectedC) '_Summary'];
Save_png(Fig,fullfile(save_dir,[file_name '.png']));
saveas(Fig,fullfile(save_dir,[file_name '.fig']),'fig');
disp(['- Plot successfully saved as ' fullfile(save_dir,file_name)]);
end

