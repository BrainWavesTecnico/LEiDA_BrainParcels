function [dice_V_yeo7,overlap_yeo7,dice_max,frac_V_yeo7] = Overlap_LEiDA_Yeo(parcellation,n_areas,centroids,rangeK,relative_threshold,min_dice,frac_threshold)
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
% Each FC state is assigned to a network (and colored with it) if:
%   1) its weighted Dice with the network of largest Dice is at least
%      relative_threshold times the largest Dice that any set of parcels of
%      the parcellation can reach for that network (dice_max; coarse
%      parcels limit the overlap with some networks), and at least
%      min_dice; or otherwise
%   2) at least frac_threshold of the FC state (weighted, over all its
%      voxels in the brain) lies within one network. The Dice coefficient
%      is low when a FC state is much smaller than the network (e.g. a few
%      occipital parcels inside the Visual network), although the state is
%      clearly part of it; this criterion assigns such focal states.
% FC states mostly outside the 7 cortical networks (e.g. subcortical)
% are not assigned (shown in black).
%
% INPUT:
% parcellation        parcellation used to segment the brain (atlas in
%                     ParcelsMNI2mm.mat)
% n_areas             number of areas of the parcellation considered
% centroids           centroids obtained from K-means (Kmeans_results)
% rangeK              range of clustering solutions
% relative_threshold  minimum Dice relative to dice_max (default 0.4)
% min_dice            minimum Dice (default 0.15)
% frac_threshold      minimum fraction of the FC state within a network
%                     (default 0.5)
%
% OUTPUT:
% dice_V_yeo7     (n_K x max(K) x 7) weighted Dice coefficient between each
%                 FC state and each Yeo network
% overlap_yeo7    (n_K x max(K) x 7) true for the network assigned to each
%                 FC state (FC states colored with that network)
% dice_max        (1 x 7) largest Dice that any set of parcels can reach
%                 for each network with this parcellation
% frac_V_yeo7     (n_K x max(K) x 7) fraction of each FC state (weighted)
%                 within each network
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, University of Minho, miguel.farinha@ccabraga.org

if nargin < 5 || isempty(relative_threshold)
    relative_threshold = 0.4;
end
if nargin < 6 || isempty(min_dice)
    min_dice = 0.15;
end
if nargin < 7 || isempty(frac_threshold)
    frac_threshold = 0.5;
end
N_Yeo = 7;

% Parcels of the parcellation and Yeo networks in MNI 2mm space
Atlas = Load_atlas(parcellation);
V_Parcels = Atlas.V;
S = load('ParcelsMNI2mm.mat','V_Yeo7');
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
% Number of voxels of each parcel in the whole brain (including voxels
% outside the 7 cortical networks)
n_vox_brain = histcounts(double(V_Parcels(V_Parcels > 0 & V_Parcels <= n_areas)),0.5:1:n_areas+0.5)';

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
frac_V_yeo7 = zeros(length(rangeK),max(rangeK),N_Yeo);

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

        % Fraction of the FC state within each network
        frac_V_yeo7(k,c,:) = (w*n_vox)/max(w*n_vox_brain,eps);

        % 1) Network with the largest Dice, if large enough relative to
        % what the parcellation allows
        [max_dice,net] = max(dice_V_yeo7(k,c,:));
        if max_dice >= relative_threshold*dice_max(net) && max_dice >= min_dice
            overlap_yeo7(k,c,net) = true;
        else
            % 2) Network containing most of the FC state
            [max_frac,net] = max(frac_V_yeo7(k,c,:));
            overlap_yeo7(k,c,net) = max_frac >= frac_threshold;
        end
    end
end
