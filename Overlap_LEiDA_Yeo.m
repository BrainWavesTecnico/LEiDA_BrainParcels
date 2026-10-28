function [dice_V_yeo7,overlap_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,centroids,rangeK,dice_threshold)
%
% Overlap between the LEiDA centroids (FC states) and the 7 resting-state
% networks of Yeo et al. (2011), measured with the Dice coefficient in
% voxel space, as in LEiDA_Voxel (Plot_Mode_TransparentBrain.m):
%   - each centroid is mapped to MNI 2mm space: a voxel belongs to the FC
%     state if its parcel has a positive value in the centroid vector
%     (the same voxels rendered in the brain plots);
%   - within the cortical voxels of the Yeo parcellation (networks 1 to 7),
%     Dice = 2|state & network| / (|state| + |network|) for each network.
% Subcortical and cerebellar voxels (labels 8 and 9 of V_Yeo7) are not
% considered, as they are not part of the 7 networks.
%
% INPUT:
% parcellation    parcellation used to segment the brain (atlas in
%                 ParcelsMNI2mm.mat)
% n_areas         number of areas of the parcellation considered
% centroids       centroids obtained from K-means (Kmeans_results)
% rangeK          range of clustering solutions
% dice_threshold  minimum Dice coefficient to assign a FC state to its
%                 best-overlapping network (default 0.3)
%
% OUTPUT:
% dice_V_yeo7     (n_K x max(K) x 7) Dice coefficient between each FC state
%                 and each Yeo network
% overlap_yeo7    (n_K x max(K) x 7) true for the network with the largest
%                 Dice coefficient of each FC state, if it is at least
%                 dice_threshold (FC states colored with that network)
%
% Note: the overlap is limited by how well the parcels match the Yeo
% networks. E.g. with the first 90 areas of AAL116, the largest Dice that
% any set of parcels can reach is about 0.82 (Visual), 0.73 (Somatomotor),
% 0.29 (Dorsal Attention), 0.48 (Ventral Attention), 0.58 (Limbic), 0.56
% (Frontoparietal) and 0.62 (Default Mode).
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

if nargin < 5 || isempty(dice_threshold)
    dice_threshold = 0.3;
end
N_Yeo = 7;

% Parcels of the parcellation and Yeo networks in MNI 2mm space
S = load('ParcelsMNI2mm.mat',['V_' parcellation],'V_Yeo7');
V_Parcels = S.(['V_' parcellation]);
V_Yeo = S.V_Yeo7;

% Cortical voxels of the Yeo networks
ind_Yeo = V_Yeo > 0 & V_Yeo <= N_Yeo;
Yeo_vox = double(V_Yeo(ind_Yeo));
Parcel_vox = double(V_Parcels(ind_Yeo));

% Number of voxels of each parcel in each network, and of each network
n_vox = zeros(n_areas,N_Yeo);
for n = 1:n_areas
    n_vox(n,:) = histcounts(Yeo_vox(Parcel_vox == n),0.5:1:N_Yeo+0.5);
end
n_vox_net = histcounts(Yeo_vox,0.5:1:N_Yeo+0.5);

dice_V_yeo7 = zeros(length(rangeK),max(rangeK),N_Yeo);
overlap_yeo7 = false(length(rangeK),max(rangeK),N_Yeo);

for k = 1:length(rangeK)
    for c = 1:rangeK(k)
        % Parcels with positive values in the centroid vector
        V = centroids{k}.C(c,1:n_areas);
        if sum(V) > 0 % for eigenvectors with sum of entries > 0
            V = -V;
        end
        active = V > 0;

        % Dice coefficient with each network
        n_overlap = sum(n_vox(active,:),1);
        n_state = sum(n_vox(active,:),'all');
        dice_V_yeo7(k,c,:) = 2*n_overlap./(n_state + n_vox_net);

        % Network with the largest overlap
        [max_dice,net] = max(dice_V_yeo7(k,c,:));
        overlap_yeo7(k,c,net) = max_dice >= dice_threshold;
    end
end
