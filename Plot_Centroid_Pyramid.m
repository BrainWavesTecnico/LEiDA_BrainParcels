function Plot_Centroid_Pyramid(data_dir,cond,parcellation,n_areas,cortex_dir,stat_of_interest,design_dir)
%
% Plot of the rendered centroids on cortical space and overlap with Yeo
% functional brain networks (Yeo et al. 2011). If the fractional occupancy
% statistics are available (LEiDA_Stats_FracOccup.mat), the mean +- SEM
% fractional occupancy of each condition is plotted next to each brain, and
% the error bars are colored according to the significance of the
% difference between the pair of conditions stat_of_interest, or between
% any pair of conditions: dark red p<0.05, dark green p<0.05/K, dark blue
% p<0.05/sum(rangeK)
%
% INPUT:
% data_dir      directory where the centroids obtained for each clustering
%               solution are saved
% cond          tags of each condition considered in the experiment
% parcellation  parcellation used to segment the brain (atlas available
%               in ParcelsMNI2mm.mat); if empty or not available, the
%               plot is skipped
% n_areas       number of areas of the parcellation to be considered
%               ([] to use the number of areas of the centroids)
% cortex_dir    direction of observation of the rendering of the FC states
% stat_of_interest  index of the pair of conditions whose p-values are
%               shown (pairs ordered 1-2, 1-3, ..., 2-3, ...), or 'any'
%               to color the states with a significant difference between
%               any pair of conditions: the smallest p-value over all
%               pairs, Bonferroni-corrected for the number of pairs.
%               Default: 1 with 2 conditions, 'any' with more
% design_dir    directory with the statistics of the study design
%               (LEiDA_Stats_FracOccup.mat), where the figure is saved
%               (default data_dir)
%
% OUTPUT:
% Fig           FC states rendered on a transparent cortex with their
%               corresponding overlap with the resting-state networks
%               defined by Yeo et al., 2011
%
% Authors: Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

% File with the clustering results (output from LEiDA_cluster.m)
file_clusters = 'LEiDA_Clusters.mat';

% The FC states can only be rendered with an atlas in ParcelsMNI2mm.mat
N_parcels = Atlas_Nparcels(parcellation);
if isempty(N_parcels)
    disp('Rendering of the centroids skipped (parcellation not available in ParcelsMNI2mm.mat)')
    return
end

% Load required data:
load([data_dir file_clusters],'Kmeans_results', 'rangeK');
if nargin < 7 || isempty(design_dir)
    design_dir = data_dir;
end

if isempty(n_areas)
    n_areas = size(Kmeans_results{1}.C,2);
end
if n_areas > N_parcels
    disp(['Rendering of the centroids skipped: the data has ' num2str(n_areas) ...
          ' areas but ' parcellation ' has only ' num2str(N_parcels) ' parcels'])
    return
end

% Number of conditions of the experiment
n_Cond = numel(cond);

if nargin < 6 || isempty(stat_of_interest)
    if n_Cond > 2
        stat_of_interest = 'any';
    else
        stat_of_interest = 1;
    end
end

% Fractional occupancy statistics (output from LEiDA_stats_FracOccup.m)
file_FracOccup = 'LEiDA_Stats_FracOccup.mat';
plot_FO = n_Cond >= 2 && isfile([design_dir file_FracOccup]);
if plot_FO
    load([design_dir file_FracOccup], 'P', 'P_pval2sided', 'Index_Conditions', 'Subjects', 'stats_level');
    if ~exist('stats_level','var')
        stats_level = 'scan';
    end

    % Pairs of conditions (ordered as in LEiDA_compare_conditions: 1-2,
    % 1-3, ..., 2-3, ...) and p-values shown
    [cond2,cond1] = find(triu(ones(n_Cond),1)');
    n_pairs = numel(cond1);
    if ischar(stat_of_interest) || isstring(stat_of_interest)
        % Difference between any pair of conditions: smallest p-value over
        % all pairs, Bonferroni-corrected for the number of pairs
        stat_tag = strjoin(cond,'-');
        sig_desc = ['any pair of ' strjoin(cond,', ') ', smallest p over ' num2str(n_pairs) ...
                    ' pairs x ' num2str(n_pairs)];
        P_pval2sided = reshape(min(min(P_pval2sided,[],1)*n_pairs,1),size(P,2),size(P,3));
    else
        stat_tag = [cond{cond1(stat_of_interest)} '-' cond{cond2(stat_of_interest)}];
        sig_desc = [cond{cond1(stat_of_interest)} ' vs ' cond{cond2(stat_of_interest)}];
        P_pval2sided = reshape(P_pval2sided(stat_of_interest,:,:),size(P,2),size(P,3));
    end
    disp(['P-values reported for ' sig_desc])

    % Mean and SEM of the fractional occupancy in each condition (over
    % subjects, averaging the scans of each subject, if stats_level is
    % 'subject', as in the statistical tests)
    [mean_P,sem_P] = FracOccup_mean_sem(P,Index_Conditions,Subjects,stats_level,n_Cond);
    % Two subplots per FC state: brain and fractional occupancy
    n_cols = 2*rangeK(end);
else
    n_cols = rangeK(end);
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT CENTROIDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

% Compute the overlap between cluster centroids obtained for each value of
% K and the resting-state networks defined by Yeo et al. 2011
disp('Computing overlap of centroids (PL states) with Yeo resting state networks')
[dice_V_yeo7, overlap_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);

% Colors of the 7 Yeo et al. (2011) networks (defined in Yeo_colors.m)
YeoColor = Yeo_colors;

% LEiDA networks colored according to closest RSN
Atlas = Load_atlas(parcellation);
Volume = Atlas.V; 

Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

% Orientation of the cortex
disp(' ')
disp(['Rendering the centroids (using ' cortex_dir '):'])
Fig = figure('Position', get(0, 'Screensize'));
% One row per K (plus space at the bottom for labels)
L = Pyramid_layout(length(rangeK),n_cols);
for k = 1:length(rangeK)
    
    % Matrix with the cluster centroids for a given K
    VLeida = Kmeans_results{k}.C;
    disp(['- K = ' num2str(rangeK(k))])
    
    for Centroid = 1:rangeK(k)
        
        % Vector of the coordinates of the cluster centroid
        V = VLeida(Centroid,:);
        if sum(V) > 0 % for eigenvectors with sum of entries > 0
            V = -V;
        end
        %V = round(V*1e2)/1e2;
        
        % Network assigned to the FC state (if any; see Overlap_LEiDA_Yeo)
        [~,net] = max(dice_V_yeo7(k,Centroid,:) + 2*overlap_yeo7(k,Centroid,:));
        
        Centroid_Vol = zeros(size(Volume));
        
        for n = find(V > 0)
            Centroid_Vol(Volume == n) = 1;
        end 
        sregion = smooth3(Centroid_Vol > 0);
        
        % View networks
        if plot_FO
            Pyramid_axes(L,k,2*Centroid-1,1,[0.05 0.05]);
        else
            Pyramid_axes(L,k,Centroid,1,[0.05 0.05]);
        end
             
        hold on
        % First plot a transparent cortex
        cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
        reducepatch(cortexpatch,0.01);
        isonormals(scortex,cortexpatch);
        
        if overlap_yeo7(k,Centroid,net) 
            patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %'FaceAlpha', V(n))
        else
            % Color in Black if no significant overlap with any of the Yeo resting-state networks
            patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %'FaceAlpha', V(n))
        end
        
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
        
        % Plot the fractional occupancy of each condition (mean +- SEM)
        % next to the brain, colored by the significance of the difference:
        % p<0.05, p<0.05/K, p<0.05/sum(rangeK)
        if plot_FO
            sig_color = Significance_color(P_pval2sided(k,Centroid),rangeK(k),rangeK);

            % Fractional occupancy of each condition (error bars without bars)
            Pyramid_axes(L,k,2*Centroid,1,[0.2 0.15]);
            errorbar(mean_P(:,k,Centroid),sem_P(:,k,Centroid),'_','MarkerSize',2,'LineStyle','none',...
                     'Color',sig_color,'CapSize',4,'LineWidth',1);
            set(gca,'color','none','YTick',[],'XTick',1:n_Cond,'XTickLabel','','Fontsize',6,'TickLabelInterpreter','none')
            % Names of the conditions under the bottom row
            if k == length(rangeK)
                set(gca,'XTickLabel',cond)
                xtickangle(90)
            end
            xlim([.4 n_Cond+.6])
            box off
        end
    end
end

% Labels of the rows (K) and of the FC states (c)
Pyramid_labels(L,rangeK,n_cols/rangeK(end));

% Reference of the 7 Yeo networks (colors of the FC states) in the top
% right corner, in the same view
rows_reference = Plot_Yeo_reference(L,cortex_dir,scortex);

if plot_FO
    % Legend of the fractional occupancy plots
    Pyramid_legend(L,{['Fractional occupancy (mean ' char(177) ' SEM over ' stats_level 's)'], ...
                      ['of ' strjoin(cond,', ') ' (left to right)'], ...
                      ['Colors: ' sig_desc]},rangeK,rows_reference+0.3);
    save_name = ['Centroid_Pyramid_' regexprep(stat_tag,'[^a-zA-Z0-9_-]','_') '_' cortex_dir];
else
    save_name = ['Centroid_Pyramid_' cortex_dir];
end
Save_png(Fig, fullfile(design_dir, [save_name '.png']));
saveas(Fig, fullfile(design_dir, [save_name '.fig']),'fig');
disp(['- Plot successfully saved as ' save_name]);

disp(' ')
