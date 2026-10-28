function [dice_V_yeo7,overlap_yeo7,dice_max] = Overlap_LEiDA_Yeo(parcellation,n_areas,centroids,rangeK,relative_threshold,min_dice)
%
% Overlap between the LEiDA centroids (FC states) and the 7 resting-state
% networks of Yeo et al. (2011), measured with a weighted Dice coefficient
% in voxel space:
%   - each centroid is mapped to MNI 2mm space: each voxel gets the value
%     of its parcel in the centroid vector, keeping only positive values
%     (the parcels rendered in the brain plots) and dividing by the largest
%     one, so that the most phase-shifted parcels weigh 1 and weakly
%     shifted parcels weigh less;
%   - within the cortical voxels of the Yeo parcellation (networks 1 to 7),
%     weighted Dice = 2 sum(w in network) / (sum(w) + voxels of network),
%     which equals the usual Dice coefficient if all weights are 0 or 1.
% Subcortical and cerebellar voxels (labels 8 and 9 of V_Yeo7) are not
% considered, as they are not part of the 7 networks.
%
% Coarse parcels limit how well a FC state can overlap some networks (e.g.
% with the first 90 areas of AAL116 the Dorsal Attention network cannot
% reach a Dice of 0.3). Therefore, each FC state is assigned to the network
% with the largest weighted Dice only if this Dice is at least
% relative_threshold times the largest Dice that any set of parcels of the
% parcellation can reach for that network (dice_max), and at least
% min_dice.
%
% INPUT:
% parcellation        parcellation used to segment the brain (atlas in
%                     ParcelsMNI2mm.mat)
% n_areas             number of areas of the parcellation considered
% centroids           centroids obtained from K-means (Kmeans_results)
% rangeK              range of clustering solutions
% relative_threshold  minimum Dice relative to dice_max (default 0.4)
% min_dice            minimum Dice (default 0.15)
%
% OUTPUT:
% dice_V_yeo7     (n_K x max(K) x 7) weighted Dice coefficient between each
%                 FC state and each Yeo network
% overlap_yeo7    (n_K x max(K) x 7) true for the network assigned to each
%                 FC state (FC states colored with that network)
% dice_max        (1 x 7) largest Dice that any set of parcels can reach
%                 for each network with this parcellation
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

if nargin < 5 || isempty(relative_threshold)
    relative_threshold = 0.4;
end
if nargin < 6 || isempty(min_dice)
    min_dice = 0.15;
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

% Number of voxels of each parcel in each network, of each parcel, and of
% each network
n_vox = zeros(n_areas,N_Yeo);
for n = 1:n_areas
    n_vox(n,:) = histcounts(Yeo_vox(Parcel_vox == n),0.5:1:N_Yeo+0.5);
end
n_vox_parcel = sum(n_vox,2);
n_vox_net = histcounts(Yeo_vox,0.5:1:N_Yeo+0.5);

% Largest Dice that any set of parcels can reach for each network: the
% best set contains the parcels with the largest fraction of their voxels
% in the network, so it is found among the sets of the first parcels
% sorted by that fraction
dice_max = zeros(1,N_Yeo);
for net = 1:N_Yeo
    [~,order] = sort(n_vox(:,net)./max(n_vox_parcel,1),'descend');
    dice_sets = 2*cumsum(n_vox(order,net))./(cumsum(n_vox_parcel(order)) + n_vox_net(net));
    dice_max(net) = max(dice_sets);
end

dice_V_yeo7 = zeros(length(rangeK),max(rangeK),N_Yeo);
overlap_yeo7 = false(length(rangeK),max(rangeK),N_Yeo);

for k = 1:length(rangeK)
    for c = 1:rangeK(k)
        % Strength of the phase shift of each parcel (positive values of
        % the centroid vector, divided by the largest)
        V = centroids{k}.C(c,1:n_areas);
        if sum(V) > 0 % for eigenvectors with sum of entries > 0
            V = -V;
        end
        w = max(V,0);
        if max(w) > 0
            w = w/max(w);
        end

        % Weighted Dice coefficient with each network
        dice_V_yeo7(k,c,:) = 2*(w*n_vox)./(w*n_vox_parcel + n_vox_net);

        % Network with the largest overlap, if large enough relative to
        % what the parcellation allows
        [max_dice,net] = max(dice_V_yeo7(k,c,:));
        overlap_yeo7(k,c,net) = max_dice >= relative_threshold*dice_max(net) && max_dice >= min_dice;
    end
end
