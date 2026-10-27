function Plot_centroid_links_FO(data_dir,selectedK,selectedC,parcellation,Subcortical)
%
% Plot a 1x3 figure for a single centroid (PL state):
%   (1) the strongest links (link_pct % of all possible links) between the
%       areas with positive values in the centroid vector in a glass
%       cortex, top view
%   (2) the same links in a glass cortex, side view
%   (3) boxplot of the fractional occupancy (FO) values by condition
%
% INPUT:
% data_dir      directory where LEiDA results are stored. The figure is
%               saved in subfolder K<selectedK>c<selectedC> (e.g., K18c10)
% selectedK     K defined by the user
% selectedC     centroid (PL state) to plot, between 1 and selectedK
% parcellation  parcellation template used to segment the brain
% Subcortical   1 to add subcortical areas to the parcellation, 0 otherwise
%
% OUTPUT:
% .fig/.png     Figure with the links of the selected centroid in a glass
%               brain (top and side views) and the boxplot of its FO
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';
% File with results for the fractional occupancy (output from LEiDA_stats_FracOccup.m)
file_P = 'LEiDA_Stats_FracOccup.mat';

% Load required data:
if isfile([data_dir file_cluster])
    load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
else
    error(['File ' file_cluster ' not found in ' data_dir]);
end
if isfile([data_dir file_P])
    load([data_dir file_P], 'cond', 'P', 'P_pval2sided', 'Index_Conditions');
else
    error(['File ' file_P ' not found in ' data_dir]);
end

if ~any(rangeK == selectedK)
    error(['K = ' num2str(selectedK) ' is not in the range of K computed']);
end
if selectedC < 1 || selectedC > selectedK || selectedC ~= round(selectedC)
    error(['Centroid must be an integer between 1 and ' num2str(selectedK)]);
end

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Matrix (selectedK*n_areas) where each row represents one PL state
V = Kmeans_results{rangeK == selectedK}.C;

% Select row selectedC of the set of centroids and scale it by its maximum value
Vc = V(selectedC,:);
Vc = Vc./max(abs(Vc));

clear Kmeans_results

% Percentage of all possible links (n_areas*(n_areas-1)/2) to show in the
% glass cortex, keeping only the strongest ones
link_pct = 2;

% Get centers of gravity of each brain area from the parcellation provided
Volume = struct2array(load('ParcelsMNI2mm',['V_' parcellation]));
n_parcels = max(Volume(:));
% Description of the parcellation for the figure title
% (e.g., 'Glasser 360 cortical + 19 subcortical parcels')
parcel_name = regexprep(parcellation,'\d+$','');
if Subcortical
    Vol_sub= struct2array(load('ParcelsMNI2mm','V_Subcortical'));
    parcel_title = [parcel_name ' ' num2str(n_parcels) ' cortical + ' num2str(max(Vol_sub(:))) ' subcortical parcels'];
    Vol_sub(Vol_sub>0)= Vol_sub(Vol_sub>0) + n_parcels;
    Volume(Vol_sub>0) = Vol_sub(Vol_sub>0);
else
    parcel_title = [parcel_name ' ' num2str(n_parcels) ' parcels'];
end
n_areas = max(Volume(:));
n_links = round(link_pct/100 * n_areas*(n_areas-1)/2);
sz = size(Volume);
MNI_coord = zeros(max(Volume(:)),3);
for n = 1:max(Volume(:))
    ind_N = find(Volume == n);
    [X,Y,Z] = ind2sub(sz,ind_N);
    MNI_coord(n,:) = mean([Y,X,Z]);
end
Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);


disp(' ');
disp(['Plotting links in cortex and FO boxplot for K = ' num2str(selectedK) ', centroid ' num2str(selectedC) ':'])
Fig = figure('Position', [100 100 1200 450]);
t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% Panel 1: links in cortex, top view
nexttile
plot_links_in_cortex(Vc, MNI_coord, scortex, n_links)
view(-90,90)
camlight;

% Panel 2: links in cortex, side view
nexttile
plot_links_in_cortex(Vc, MNI_coord, scortex, n_links)
view(0,0)
camlight;
% Shift the side view slightly to the left within its panel
xlim([5 105] + 10)
set(gca,'Clipping','off')

% Panel 3: boxplot of fractional occupancy by condition
nexttile

P_cond = cell(1,n_Cond);
g = cell(1,n_Cond);
for j = 1:n_Cond
    P_cond{j} = P(Index_Conditions == j,rangeK == selectedK,selectedC);
    g{j} = repmat(cond(j), length(P(Index_Conditions == j,rangeK == selectedK,selectedC)),1);
end
P_data = vertcat(P_cond{:});
g_data = vertcat(g{:});

bx = boxplot(P_data, g_data, 'Symbol', 'k.','OutlierSize',7, 'Widths',0.5);
set(bx,{'linew'},{1})
set(bx,'Color','k')   % Boxes, whiskers and medians in black

hold on
set(gca,'XTickLabel',cond)
ax = gca;
ax.XAxis.FontSize = 18;   % Condition labels (e.g., CN, MCI)
ax.YAxis.FontSize = 12;
ax.FontAngle = 'normal';
ylabel('Fractional Occupancy', 'Fontsize', 18);
set(gca,'color','none','LineWidth',1)

X_locations = zeros(n_Cond*(n_Cond-1)/2,2);
cond_pair = 1;
for cond1 = 1:n_Cond-1
    for cond2 = cond1+1:n_Cond
        X_locations(cond_pair,1) = cond1;
        X_locations(cond_pair,2) = cond2;
        cond_pair = cond_pair + 1;
    end
end
X_locations(:,1) = X_locations(:,1) + 0.1;
X_locations(:,2) = X_locations(:,2) - 0.1;
Max_Y =  max(P_data);
Y_LOCATIONS = Max_Y + 0.05*(abs(X_locations(:,1) - X_locations(:,2)));

% Significance bars with the p-value and Hedge's effect size (g) of each
% pair of conditions (shown even when not significant) and asterisks for
% the significance level:
%   *   p <= 0.05 (uncorrected)
%   **  p <= 0.05/K (corrected for the number of states in this K)
%   *** p <= 0.05/sum(rangeK) (corrected for all states across all K)
p_vals = P_pval2sided(:,rangeK == selectedK,selectedC);
for k = 1:numel(p_vals)
    plot(X_locations(k,:),[Y_LOCATIONS(k) Y_LOCATIONS(k)],'-k','LineWidth',1)

    if p_vals(k) < 0.01
        p_str = ['p = ' num2str(p_vals(k),'%.2e')];
    else
        p_str = ['p = ' num2str(p_vals(k),'%.2f')];
    end

    if p_vals(k) <= 0.05/sum(rangeK)
        stars = '***';
    elseif p_vals(k) <= 0.05/selectedK
        stars = '**';
    elseif p_vals(k) <= 0.05
        stars = '*';
    else
        stars = '';
    end

    % Effect size between the two conditions of this pair (difference of
    % means / pooled standard deviation; permutation_test.m also applies
    % the small-sample correction of Hedges' g)
    x1 = P_cond{round(X_locations(k,1))};
    x2 = P_cond{round(X_locations(k,2))};
    n1 = numel(x1);
    n2 = numel(x2);
    s_pooled = sqrt(((n1-1)*var(x1) + (n2-1)*var(x2)) / (n1+n2-2));
    if s_pooled > 0
        eff = abs(mean(x1) - mean(x2)) / s_pooled;
    else
        eff = 0;
    end
    eff_str = ['g = ' num2str(eff,'%.2f')];

    if isempty(stars)
        label = {p_str, eff_str};
    else
        label = {p_str, eff_str, stars};
    end
    text(mean(X_locations(k,:)), Y_LOCATIONS(k), label, ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','Fontsize',14)
end

ylim([0 max(Y_LOCATIONS)*1.05+0.25])

hold off
box off

title(t, parcel_title, 'Fontsize', 20, 'Interpreter', 'none')

% Save the figure in folder K<selectedK>c<selectedC> (e.g., K18c10) inside data_dir
save_dir = fullfile(data_dir, ['K' num2str(selectedK) 'c' num2str(selectedC)]);
if ~isfolder(save_dir)
    mkdir(save_dir);
end
file_name = ['K' num2str(selectedK) '_C' num2str(selectedC) '_LinksCortex_FO'];
saveas(Fig, fullfile(save_dir, [file_name '.png']),'png');
saveas(Fig, fullfile(save_dir, [file_name '.fig']),'fig');
disp(['- Plot successfully saved as ' fullfile(save_dir, file_name)]);


end


function plot_links_in_cortex(Vc, MNI_coord, scortex, n_links)
% Plot the glass cortex, the n_links strongest links between regions with
% positive value in the centroid vector and the nodes colored by their
% contribution, in the current axes. The view and camlight are set by the
% caller.

hold on

% Cortex patch
cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.1);
reducepatch(cortexpatch,0.5);
isonormals(scortex,cortexpatch);

% Plot the strongest links between regions with positive value in centroid
% vector, where the strength of the link between areas n and p is Vc(n)*Vc(p)
n_strong = find(Vc > 0);

if numel(n_strong) > 1
    pairs = nchoosek(n_strong,2);
    link_strength = Vc(pairs(:,1)).*Vc(pairs(:,2));
    [~, order] = sort(link_strength,'descend');
    pairs = pairs(order(1:min(n_links,numel(order))),:);
    for l = 1:size(pairs,1)
        c1 = MNI_coord(pairs(l,1),:);
        c2 = MNI_coord(pairs(l,2),:);
        plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)],'Color',[0.8 0 0],'LineWidth',1);
    end
end

a = 2;
[x,y,z] = sphere;
x = a*x;
y = a*y;
z = a*z;

% Scale positive and negative values separately so that the most positive
% node is pure red and the most negative node is pure blue
max_pos = max([Vc(Vc > 0) eps]);
max_neg = max([-Vc(Vc < 0) eps]);

for n = 1:length(Vc)
    if Vc(n) > 0 % color nodes with positive contribution in a gradient of red colors
        surf(x+MNI_coord(n,1), y+MNI_coord(n,2),z+MNI_coord(n,3),'FaceColor',[Vc(n)/max_pos 0 0],'EdgeColor','none','FaceAlpha',0.7);
    elseif Vc(n) < 0 % color nodes with negative contribution in a gradient of blue colors
        surf(x+MNI_coord(n,1), y+MNI_coord(n,2),z+MNI_coord(n,3),'FaceColor',[0 0 -Vc(n)/max_neg],'EdgeColor','none','FaceAlpha',0.7);
    end
end

axis off;
axis equal;

material dull;
lighting gouraud;
daspect([1 1 1])
xlim([5 105]); ylim([5 85]); zlim([0 80])

hold off

end
