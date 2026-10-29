function Export_centroids_nifti(P)
%
% Save the centroids (FC states) of each K as NIfTI files in MNI space, to
% view them in other software (e.g. FSLeyes, MRIcroGL) or use them in
% other analyses: one 4D file per K, K<K>_centroids.nii.gz, with one volume
% per FC state, in folder Centroids_NIfTI of the results of the run.
%
% Each voxel takes the value of its area in the centroid vector (sign as in
% the figures: positive values are the areas of the FC state, which are
% out of phase with the rest of the brain). The files are in the grid of
% the atlas when the atlas is a NIfTI file, or in the FSL MNI152 2mm grid
% for the atlases of ParcelsMNI2mm.mat.
%
% INPUT:
% P         parameters set in LEiDA_Main (after LEiDA_Start)
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

Atlas = Load_atlas(P.Parcellation);
if isempty(Atlas)
    disp('Centroids not saved as NIfTI: the parcellation is not a NIfTI file nor an atlas of ParcelsMNI2mm.mat')
    return
end
load([P.leida_results 'LEiDA_Clusters.mat'],'Kmeans_results','rangeK');
n_areas = size(Kmeans_results{1}.C,2);

% Header of the output files: grid of the atlas
if endsWith(lower(P.Parcellation),{'.nii','.nii.gz'})
    info = niftiinfo(P.Parcellation);
else
    info = niftiinfo('MNI152_T1_2mm_brain_mask.nii');
end
labels = Atlas.V_native;
labels(labels > n_areas) = 0;
in_atlas = labels > 0;

save_dir = [P.leida_results 'Centroids_NIfTI' filesep];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

disp(' ')
disp(['Saving the centroids as NIfTI files in ' save_dir])
for k = 1:length(rangeK)
    C = Kmeans_results{k}.C;
    Vol = zeros([size(labels) rangeK(k)],'single');
    for c = 1:rangeK(k)
        Vc = C(c,:);
        if sum(Vc) > 0 % same sign as in the figures (sum of entries < 0)
            Vc = -Vc;
        end
        vol_c = zeros(size(labels),'single');
        vol_c(in_atlas) = Vc(labels(in_atlas));
        Vol(:,:,:,c) = vol_c;
    end

    info_K = info;
    info_K.ImageSize = size(Vol,1:4);
    info_K.PixelDimensions = [info.PixelDimensions(1:3) 1];
    info_K.Datatype = 'single';
    info_K.BitsPerPixel = 32;
    info_K.MultiplicativeScaling = 1;
    info_K.AdditiveOffset = 0;
    info_K.Description = sprintf('LEiDA centroids K=%d (%s)',rangeK(k),Atlas.name);
    niftiwrite(Vol,[save_dir 'K' num2str(rangeK(k)) '_centroids'],info_K,'Compressed',true);
end
disp(['- Centroids for K = ' num2str(rangeK(1)) ' to ' num2str(rangeK(end)) ' saved as K<K>_centroids.nii.gz'])
