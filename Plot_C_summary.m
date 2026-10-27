function Plot_C_summary(data_dir,save_dir,selectedK,centroid,parcellation,n_areas,Subcortical)
%
% Plot the centroids in 3D glass brain, in vector format and the boxplots
% of the fractional occupancy values
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% centroid      centroid defined by the user
% parcellation  parcellation template used to segment the brain
%
% OUTPUT:
% .fig/.png     Plot of summary information for the selected centroid
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';
% File with results for the fractional occupancy (output from LEiDA_stats_FracOccup.m)
file_P = 'LEiDA_Stats_FracOccup.mat';

% Load required data:
load([data_dir file_cluster], 'Kmeans_results', 'rangeK');

load([data_dir file_P], 'cond', 'P', 'P_pval2sided', 'Index_Conditions');

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Matrix of dimension selectedK*90 where each row represents one FC state
V = Kmeans_results{rangeK == selectedK}.C(centroid,:);
% Scale each cluster centroid by its maximum value and transpose the matrix
V = V'./max(abs(V'));


% Colors of the 7 Yeo et al. (2011) networks (defined in Yeo_colors.m)
[YeoColor,Yeo_names] = Yeo_colors;

% Compute the overlap between cluster centroids obtained for each value of
% K and the resting-state networks defined by Yeo et al. 2011
[cc_V_yeo7, p_V_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,Kmeans_results,rangeK);
clear Kmeans_results

% LEiDA networks colored according to closest RSN
Volume = struct2array(load('ParcelsMNI2mm',['V_' parcellation]));
if Subcortical
    Vol_sub= struct2array(load('ParcelsMNI2mm','V_Subcortical')); 
    Vol_sub(Vol_sub>0)= Vol_sub(Vol_sub>0) + n_areas;
    Volume(Vol_sub>0) = Vol_sub(Vol_sub>0);
end
Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

% Get labels from parcellation used
switch parcellation
    case 'AAL116'
        load('ParcelsMNI2mm','label116');
        labels = label116([1:n_areas],:);
        clear label116
    case 'AAL120'
        load('ParcelsMNI2mm','label120');
        labels = label120([1:n_areas],:);
        clear label120    
    case 'Schaefer100'
        labels = struct2array(load('ParcelsMNI2mm','labelSchaefer100P7Networks'));
    case 'Desikan'
        labels = struct2array(load('ParcelsMNI2mm','labelDesikan70'));
end

if Subcortical
        load('ParcelsMNI2mm','label_Glasser_19subcortical')
        labels=cat(1,labels,label_Glasser_19subcortical);
end



disp(' ');
disp(['Plotting summary information for PL state ' num2str(centroid) ':'])
Fig = figure('Position', get(0, 'Screensize')); 
    
% Plot of vector with parcels ordered by contribution
[Vo, asc_ind] = sort(V);
numMarkers = size(Vo,1);
if centroid == 1
    markerColors = jet(numMarkers*2);
else
    markerColors = jet(numMarkers);
end
miny = min(Vo);
maxy = max(Vo);
colorMapRows = round(rescale((Vo - miny) / (maxy - miny), 1, numMarkers));
subplot(2,4,[1 5])
for k = 1:length(Vo)
    thisMarkerColor = markerColors(colorMapRows(k), :);
    barh(k,Vo(k),'FaceColor',thisMarkerColor,'EdgeColor','none','Barwidth',.5)
    hold on;
    ylim([0 n_areas+1])
    xlim([-1 1])
    set(gca,'XTick',[-1 0 1])
    set(gca,'XTickLabel',[-1 0 1])
    set(gca,'YTick',1:n_areas)
    set(gca, 'YTickLabel', labels(asc_ind,:), 'FontSize', 6, 'TickLabelInterpreter', 'none')
    ax = gca;
    ax.XAxis.FontSize = 10;
    grid on;
    title({'Signal Phase','projection'},'Fontsize',10)
end
    
% Plot of centroid in 3D glass brain
[~, net] = max(cc_V_yeo7(rangeK == selectedK,centroid,:));
subplot(2,4,2)
hold on
cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
reducepatch(cortexpatch,0.01);
isonormals(scortex,cortexpatch);

% Color all areas with positive V with transparency proportional to
% contribution
n_pos = find(V > 0);
if numel(n_pos) > 0
    for region = 1:length(n_pos)
        region_Vol = Volume == n_pos(region); % 3D volume with 1 only in the voxels of that area
        sregion = smooth3(region_Vol > 0);

        if p_V_yeo7(rangeK == selectedK,centroid,net) < 0.05/selectedK
            patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
        else
            patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
        end
    end
else % Plot of the global mode with "gray/transparent" color
    region_Vol = Volume > 0 & Volume <= n_areas; % 3D volume with 1 only in the voxels of that area
    sregion = smooth3(region_Vol > 0);
    patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
end

title(['PL state ' num2str(centroid)], 'Fontsize', 10)    
material dull; lighting gouraud
daspect([1 1 1])
view(-90,90)    % Top view    Side view:   view(0,0)
camlight;
xlim([5 105]); ylim([5 85]); zlim([0 80])
axis off
    
subplot(2,4,3)
hold on
cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
reducepatch(cortexpatch,0.01);
isonormals(scortex,cortexpatch);

% Color all areas with positive V with transparency proportional to
% contribution
n_pos = find(V > 0);
if numel(n_pos) > 0
    for region = 1:length(n_pos)
        region_Vol = Volume == n_pos(region); % 3D volume with 1 only in the voxels of that area
        sregion = smooth3(region_Vol > 0);

        if p_V_yeo7(rangeK == selectedK,centroid,net) < 0.05/selectedK
            patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
        else
            patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
        end
    end
else % Plot of the global mode with "gray/transparent" color
    region_Vol = Volume > 0 & Volume <= n_areas; % 3D volume with 1 only in the voxels of that area
    sregion = smooth3(region_Vol > 0);
    patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none') %,'FaceAlpha', V(n_pos(region),c))
end

title(['PL state ' num2str(centroid)], 'Fontsize', 10)    
material dull; lighting gouraud
daspect([1 1 1])
view(0,0)    % Top view    Side view:   view(0,0)
camlight;
xlim([5 105]); ylim([5 85]); zlim([0 80])
axis off

    
% Plot of matrix format
subplot(2,4,4)
colormap(jet)
imagesc(V*V',[-1 1])
set(gca, 'XTick',10:20:n_areas, 'YTick',10:20:n_areas)
axis square
title(['V_{C_{' num2str(centroid) '}}' '.V_{C_{' num2str(centroid) '}}' '^T'], 'Fontsize', 10)
    
    
% Plot of boxplot of fractional occupancy values by condition
subplot(2,4,8)
P_cond = cell(1,n_Cond);
g = cell(1,n_Cond);
for j = 1:n_Cond
    P_cond{j} = P(Index_Conditions == j,rangeK == selectedK,centroid);
    g{j} = repmat(cond(j), length(P(Index_Conditions == j,rangeK == selectedK,centroid)),1);
end
P_data = vertcat(P_cond{:});
g_data = vertcat(g{:});
    
bx = boxplot(P_data, g_data, 'Symbol', 'k.','OutlierSize',7, 'Widths',0.5);
set(bx,{'linew'},{1})
b = get(get(gca,'children'),'children');   % Get the handles of all the objects
f = get(b,'tag');   % List the names of all the objects
med_col = b(n_Cond+1:n_Cond*2);
set(med_col, 'Color', 'k');
box_col = b(n_Cond*2+1:n_Cond*3);
if p_V_yeo7(rangeK == selectedK,centroid,net) < 0.05/selectedK
    set(box_col, 'Color', YeoColor(net,:));
else
    set(box_col, 'Color', 'k');
end   
    
hold on
set(gca,'XTickLabel',cond,'Fontsize',14)
xtickangle(30)
ylabel('Fractional Occupancy (%)','Fontsize',14);
set(gca,'color','none')
    
X_locations = zeros(n_Cond*(n_Cond-1)/2,2);
cond_pair = 1;
for cond1 = 1:n_Cond-1
    for cond2 = cond1+1:n_Cond
        X_locations(cond_pair,1) = cond1;
        X_locations(cond_pair,2) = cond2;
        cond_pair = cond_pair + 1;
    end
end
% X_locations = sortrows(X_locations,2);
X_locations(:,1) = X_locations(:,1) + 0.1;
X_locations(:,2) = X_locations(:,2) - 0.1;
Max_Y =  max(P_data);
Y_LOCATIONS = Max_Y + 0.05*(abs(X_locations(:,1) - X_locations(:,2)));

% Red asterisks
asterisks = find(P_pval2sided(:,rangeK == selectedK,centroid) <= 0.05);
plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-r', 'LineWidth',1)
plot(mean(X_locations(asterisks,:),2) - 0.05, Y_LOCATIONS(asterisks)*1.05, '*','Color','r','Markersize',10)

if asterisks
    pvals = P_pval2sided(asterisks, rangeK == selectedK, centroid);
    pstr  = arrayfun(@(p) sprintf('%.*f', max(2, ceil(-log10(p))+1), p), pvals, 'UniformOutput', false);
    text(mean(X_locations(asterisks,:), 2) + 0.05, Y_LOCATIONS(asterisks)*1.05, ...
        strcat('p=', pstr), 'FontSize', 16)
end

% Green asterisks
asterisks = find(P_pval2sided(:,rangeK == selectedK,centroid) <= 0.05/selectedK);
plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-g', 'LineWidth',1)
plot(mean(X_locations(asterisks,:),2) - 0.05, Y_LOCATIONS(asterisks)*1.05, '*','Color','g','Markersize',10)


% Blue asterisks
asterisks = find(P_pval2sided(:,rangeK == selectedK,centroid) <= 0.05/sum(rangeK));
plot(X_locations(asterisks,:)',[Y_LOCATIONS(asterisks) Y_LOCATIONS(asterisks)]','-b', 'LineWidth',1)
plot(mean(X_locations(asterisks,:),2) - 0.05, Y_LOCATIONS(asterisks)*1.05, '*b','Markersize',10)



ylim([0 max(Y_LOCATIONS)+0.15])

hold off
box off

%%
subplot(2,4,6)  % Plot Links

Vc=V;

sz = size(Volume);
MNI_coord = zeros(max(Volume(:)),3);
for n = 1:max(Volume(:))
    ind_N = find(Volume == n);
    [X,Y,Z] = ind2sub(sz,ind_N);
    MNI_coord(n,:) = mean([Y,X,Z]);
end
Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

    hold on

    % Cortex patch
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.1);
    reducepatch(cortexpatch,0.5);
    isonormals(scortex,cortexpatch);

    % Plot the links between regions with positive value in centroid vector
    n_strong = find(Vc > 0);
    % Change treshold to show links in the plot  (e.g., Vc > 0.5 for less links)

    if numel(n_strong) > 1
        for a = 1:numel(n_strong)-1
            n = n_strong(a);
            for b = a+1:numel(n_strong)
                p = n_strong(b);
                c1 = [MNI_coord(n,1) MNI_coord(n,2) MNI_coord(n,3)];
                c2 = [MNI_coord(p,1) MNI_coord(p,2) MNI_coord(p,3)];
                plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)],'Color','r','LineWidth',1);
            end
        end
    end

    a = 2;
    [x,y,z] = sphere;
    x = a*x;
    y = a*y;
    z = a*z;

    for n = 1:length(Vc)
        if Vc(n) > 0 % color nodes with positive contribution in a gradient of red colors
            surf(x+MNI_coord(n,1), y+MNI_coord(n,2),z+MNI_coord(n,3),'FaceColor',[Vc(n) 0 0],'EdgeColor','none','FaceAlpha',0.7);
        elseif Vc(n) < 0 % color nodes with negative contribution in a gradient of blue colors
            surf(x+MNI_coord(n,1), y+MNI_coord(n,2),z+MNI_coord(n,3),'FaceColor',[0 0 -Vc(n)],'EdgeColor','none','FaceAlpha',0.7);
        end
    end

    axis off;
    axis equal;

    material dull;
    lighting gouraud;
    daspect([1 1 1])
    view(-90,90)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    
    
%%
subplot(2,4,7)  % Plot Links

Vc=V;

sz = size(Volume);
MNI_coord = zeros(max(Volume(:)),3);
for n = 1:max(Volume(:))
    ind_N = find(Volume == n);
    [X,Y,Z] = ind2sub(sz,ind_N);
    MNI_coord(n,:) = mean([Y,X,Z]);
end
Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

    hold on

    % Cortex patch
    cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.1);
    reducepatch(cortexpatch,0.5);
    isonormals(scortex,cortexpatch);

    % Plot the links between regions with positive value in centroid vector
    n_strong = find(Vc > 0);
    % Change treshold to show links in the plot  (e.g., Vc > 0.5 for less links)

    if numel(n_strong) > 1
        for a = 1:numel(n_strong)-1
            n = n_strong(a);
            for b = a+1:numel(n_strong)
                p = n_strong(b);
                c1 = [MNI_coord(n,1) MNI_coord(n,2) MNI_coord(n,3)];
                c2 = [MNI_coord(p,1) MNI_coord(p,2) MNI_coord(p,3)];
                plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)],'Color','r','LineWidth',1);
            end
        end
    end

    a = 2;
    [x,y,z] = sphere;
    x = a*x;
    y = a*y;
    z = a*z;

    for n = 1:length(Vc)
        if Vc(n) > 0 % color nodes with positive contribution in a gradient of red colors
            surf(x+MNI_coord(n,1), y+MNI_coord(n,2),z+MNI_coord(n,3),'FaceColor',[Vc(n) 0 0],'EdgeColor','none','FaceAlpha',0.7);
        elseif Vc(n) < 0 % color nodes with negative contribution in a gradient of blue colors
            surf(x+MNI_coord(n,1), y+MNI_coord(n,2),z+MNI_coord(n,3),'FaceColor',[0 0 -Vc(n)],'EdgeColor','none','FaceAlpha',0.7);
        end
    end

    axis off;
    axis equal;

    material dull;
    lighting gouraud;
    daspect([1 1 1])
    view(0,0)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    
    


saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) 'C' num2str(centroid) '_Summary.png']),'png');
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) 'C' num2str(centroid) '_Summary.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) 'C' num2str(centroid) '_Summary']);

close all;


