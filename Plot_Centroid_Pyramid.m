function Plot_Centroid_Pyramid(data_dir,cond,parcellation,n_areas,cortex_dir,stat_of_interest)
%
% Plot of the rendered centroids on cortical space and overlap with Yeo
% functional brain networks (Yeo et al. 2011). If the fractional occupancy
% statistics are available (LEiDA_Stats_FracOccup.mat), the mean +- SEM
% fractional occupancy of each condition is plotted next to each brain, and
% the FC states with significant differences between the pair of
% conditions stat_of_interest are marked with asterisks and colors:
% * (dark red) p<0.05, ** (dark green) p<0.05/K, *** (dark blue)
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
%               shown (pairs ordered 1-2, 1-3, ..., 2-3, ...; default 1)
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

if isempty(n_areas)
    n_areas = size(Kmeans_results{1}.C,2);
end
if n_areas > N_parcels
    disp(['Rendering of the centroids skipped: the data has ' num2str(n_areas) ...
          ' areas but ' parcellation ' has only ' num2str(N_parcels) ' parcels'])
    return
end

if nargin < 6 || isempty(stat_of_interest)
    stat_of_interest = 1;
end

% Number of conditions of the experiment
n_Cond = numel(cond);

% Fractional occupancy statistics (output from LEiDA_stats_FracOccup.m)
file_FracOccup = 'LEiDA_Stats_FracOccup.mat';
plot_FO = n_Cond >= 2 && isfile([data_dir file_FracOccup]);
if plot_FO
    load([data_dir file_FracOccup], 'P', 'P_pval2sided', 'Index_Conditions', 'Subjects', 'stats_level');
    if ~exist('stats_level','var')
        stats_level = 'scan';
    end

    % Pair of conditions whose p-values are shown
    % (pairs ordered as in LEiDA_compare_conditions: 1-2, 1-3, ..., 2-3, ...)
    [cond2,cond1] = find(triu(ones(n_Cond),1)');
    stat_tag = [cond{cond1(stat_of_interest)} '-' cond{cond2(stat_of_interest)}];
    disp(['P-values reported for ' cond{cond1(stat_of_interest)} ' vs ' cond{cond2(stat_of_interest)}])
    P_pval2sided = reshape(P_pval2sided(stat_of_interest,:,:),size(P,2),size(P,3));

    % Mean and SEM of the fractional occupancy in each condition (over
    % subjects, averaging the scans of each subject, if stats_level is
    % 'subject', as in the statistical tests)
    mean_P = zeros(n_Cond,size(P,2),size(P,3));
    sem_P = zeros(n_Cond,size(P,2),size(P,3));
    for j = 1:n_Cond
        P_j = reshape(P(Index_Conditions == j,:,:),nnz(Index_Conditions == j),[]);
        if strcmpi(stats_level,'subject')
            [~,~,g] = unique(Subjects(Index_Conditions == j));
            P_j = splitapply(@(x) mean(x,1),P_j,g);
        end
        mean_P(j,:,:) = reshape(mean(P_j,1),size(P,2),size(P,3));
        sem_P(j,:,:) = reshape(std(P_j,0,1)/sqrt(size(P_j,1)),size(P,2),size(P,3));
    end
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
[cc_V_yeo7, p_V_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);

% Colors of the 7 Yeo et al. (2011) networks (defined in Yeo_colors.m)
YeoColor = Yeo_colors;

% Colors of the significance levels: p<0.05, p<0.05/K and p<0.05/sum(K)
darkRed = [0.6, 0, 0];
darkGreen = [0.3, 0.5, 0];
darkBlue = [0.07, 0.62, 1.00];
% Vertical position of the significance asterisks above each brain
% (normalized axes units: 1 is the top of the axes)
star_height = 0.95;

% LEiDA networks colored according to closest RSN
Volume = struct2array(load('ParcelsMNI2mm',['V_' parcellation])); 

Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

% Orientation of the cortex
disp(' ')
disp(['Rendering the centroids (using ' cortex_dir '):'])
Fig = figure('Position', get(0, 'Screensize'));
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
        
        [~,net] = max(cc_V_yeo7(k,Centroid,:));
        
        Centroid_Vol = zeros(size(Volume));
        
        for n = find(V > 0)
            Centroid_Vol(Volume == n) = 1;
        end 
        sregion = smooth3(Centroid_Vol > 0);
        
        % View networks
        if plot_FO
            subplot_tight(length(rangeK),n_cols,2*Centroid-1+(k-1)*n_cols,0.01)
        else
            subplot_tight(length(rangeK),n_cols,Centroid+(k-1)*n_cols,0.008)
        end
             
        hold on
        % First plot a transparent cortex
        cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
        reducepatch(cortexpatch,0.01);
        isonormals(scortex,cortexpatch);
        
        if p_V_yeo7(k,Centroid,net) < 0.05/rangeK(k) 
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
        
        % Mark the FC states with significant differences in fractional
        % occupancy with asterisks and plot the fractional occupancy of each
        % condition (mean +- SEM) next to the brain:
        % * p<0.05, ** p<0.05/K, *** p<0.05/sum(rangeK)
        if plot_FO
            pval = P_pval2sided(k,Centroid);
            if pval <= 0.05/sum(rangeK)
                stars = '***'; star_color = darkBlue;
            elseif pval <= 0.05/rangeK(k)
                stars = '**'; star_color = darkGreen;
            elseif pval <= 0.05
                stars = '*'; star_color = darkRed;
            else
                stars = ''; star_color = 'k';
            end
            if ~isempty(stars)
                % Asterisks just above the brain (vertical position in
                % normalized axes units; decrease to move them down)
                text(0.5,star_height,stars,'Units','normalized','HorizontalAlignment','center',...
                     'VerticalAlignment','middle','Fontsize',16,'FontWeight','bold','Color',star_color)
            end

            % Fractional occupancy of each condition (error bars without bars)
            subplot_tight(length(rangeK),n_cols,2*Centroid+(k-1)*n_cols,0.01)
            errorbar(mean_P(:,k,Centroid),sem_P(:,k,Centroid),'_','MarkerSize',2,'LineStyle','none',...
                     'Color',star_color,'CapSize',4,'LineWidth',1);
            set(gca,'color','none')
            set(gca,'XTick',1:n_Cond,'XTickLabel','','YTick',[])
            xlim([.4 n_Cond+.6])
            box off
        end
    end
end

if plot_FO
    % Legend of the fractional occupancy plots
    annotation('textbox',[0.005 0 0.9 0.025],'EdgeColor','none','Interpreter','none','FontSize',10,...
        'String',['Fractional occupancy (mean ' char(177) ' SEM over ' stats_level 's) of ' strjoin(cond,', ') ...
                  ' (left to right). Significance (' cond{cond1(stat_of_interest)} ' vs ' ...
                  cond{cond2(stat_of_interest)} '): * p<0.05, ** p<0.05/K, *** p<0.05/' num2str(sum(rangeK))])
    save_name = ['Centroid_Pyramid_' regexprep(stat_tag,'[^a-zA-Z0-9_-]','_') '_' cortex_dir];
else
    save_name = ['Centroid_Pyramid_' cortex_dir];
end
saveas(Fig, fullfile(data_dir, [save_name '.png']),'png');
saveas(Fig, fullfile(data_dir, [save_name '.fig']),'fig');
disp(['- Plot successfully saved as ' save_name]);

disp(' ')
