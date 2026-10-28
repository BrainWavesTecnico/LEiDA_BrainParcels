function Plot_K_links_in_cortex(data_dir,save_dir,selectedK,parcellation)
%
% Plot the strongest links between the areas with positive values in the
% centroid of each FC state in a glass cortex (top view), with the nodes
% colored by their value (red positive, blue negative).
%
% INPUT:
% data_dir      directory where LEiDA results are stored
% save_dir      directory to save results for selected optimal K
% selectedK     K defined by the user
% parcellation  parcellation template used to segment the brain
%
% OUTPUT:
% .fig/.png     Plot of the links between positive areas in the centroids
%               in a glass brain for each FC state
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

load([data_dir 'LEiDA_Clusters.mat'], 'Kmeans_results', 'rangeK');

% Matrix (selectedK x n_areas) where each row represents one FC state,
% each scaled by its largest absolute value
V = Kmeans_results{rangeK == selectedK}.C;
V = V./max(abs(V),[],2);
n_areas = size(V,2);
clear Kmeans_results

% Centers of the areas and brain mask
MNI_coord = Parcel_centers(parcellation,n_areas);
Brain_Mask = niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex = smooth3(Brain_Mask > 0);

% Strongest 2% of all possible links (as in Plot_Centroid_Summary)
n_links = round(0.02*n_areas*(n_areas-1)/2);

disp(' ');
disp('Plotting links in cortex:')
Fig = figure('Position', get(0, 'Screensize'));
for c = 1:selectedK
    subplot_tight(1,selectedK,c,0.01)
    Draw_links_in_cortex(V(c,:),MNI_coord,scortex,n_links);
    view(-90,90)
    camlight;
    if c == 1
        title('Global Mode', 'Fontsize', 12);
    else
        title(['FC state ' num2str(c)], 'Fontsize', 12)
    end
end
Save_png(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_LinksCortex.png']));
saveas(Fig, fullfile(save_dir, ['K' num2str(selectedK) '_LinksCortex.fig']),'fig');
disp(['- Plot successfully saved as K' num2str(selectedK) '_LinksCortex']);
