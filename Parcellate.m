function Parcellate(data_dir,file_pattern,save_dir,parcellation,n_areas)
%
% Parcellate fMRI data in MNI space: the time series of each area is the
% mean over the voxels of the area that have signal. One .mat file
% (variable data, n_areas x time) is saved per scan in save_dir, keeping
% the subfolder structure of data_dir (e.g. sub-XXX/ses-YYY/). Scans
% already parcellated are not computed again.
%
% The atlas is resampled to the voxel grid of each scan (nearest
% neighbour), using the affine transformations in the NIfTI headers of the
% data and of the atlas. The data can therefore have any resolution, field
% of view or orientation, as long as it is in the same space as the atlas
% (MNI). The data itself is not interpolated.
%
% INPUT:
% data_dir      directory with the 4D NIfTI files (.nii or .nii.gz); files
%               are searched in this directory and all its subfolders
% file_pattern  pattern of the file names to use, e.g. '*bold.nii.gz'
% save_dir      directory to save the parcellated data
% parcellation  atlas: NIfTI file with integer labels in MNI space (at its
%               own resolution, as delivered), or name of an atlas in
%               ParcelsMNI2mm.mat (FSL MNI152 2mm grid); see Load_atlas
% n_areas       number of areas of the atlas to consider ([] for all)
%
% Voxels of an area with constant signal (e.g. 0 outside the brain mask of
% the scan) or NaN are excluded from the mean. Areas without any voxel with
% signal are NaN (and replaced by low-amplitude noise in LEiDA_data).
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARCELLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ')

Atlas = Load_atlas(parcellation);
if isempty(Atlas)
    error('LEiDA:atlas','To parcellate data in MNI space, Parcellation must be a NIfTI atlas or an atlas of ParcelsMNI2mm.mat');
end
if isempty(n_areas)
    n_areas = Atlas.N_parcels;
end
[~,atlas_name] = fileparts(regexprep(parcellation,'\.nii(\.gz)?$',''));

% Files to parcellate
files = dir(fullfile(data_dir,'**',file_pattern));
files = files(~[files.isdir]);
if isempty(files)
    error('LEiDA:nodata','No files matching ''%s'' were found in %s', file_pattern, data_dir);
end
data_root = [regexprep(data_dir,'[\\/]+$','') filesep];

% Atlas in the grid of the data (recomputed only if the grid changes)
grid = [];
V_data = [];

for s = 1:numel(files)

    in_file = fullfile(files(s).folder,files(s).name);
    % Same subfolders as in data_dir
    out_dir = fullfile(save_dir,erase([files(s).folder filesep],data_root));
    base_name = regexprep(files(s).name,'\.nii(\.gz)?$','');
    out_file = fullfile(out_dir,[base_name '_' atlas_name '.mat']);
    if isfile(out_file)
        disp(['Already parcellated: ' out_file])
        continue
    end
    if ~exist(out_dir,'dir')
        mkdir(out_dir);
    end

    disp(['Parcellating ' in_file ' using parcellation ' atlas_name]);
    info = niftiinfo(in_file);
    A_data = Nifti_affine(info);
    size_data = info.ImageSize(1:3);

    % Atlas in the voxel grid of this scan
    if ~isequal(grid,{size_data,A_data})
        grid = {size_data,A_data};
        V_data = Resample_labels(Atlas.V_native,Atlas.A_native,size_data,A_data);
        V_data(V_data > n_areas) = 0;
        missing_areas = setdiff(1:n_areas,unique(V_data(V_data > 0)));
        disp(['      - Voxel size ' mat2str(info.PixelDimensions(1:3),3) ' mm, grid ' mat2str(size_data)])
        if ~isempty(missing_areas)
            disp(['      - Attention: areas without voxels in the field of view of the data: ' mat2str(missing_areas)])
        end
        ind_atlas = find(V_data > 0);
        label_atlas = V_data(ind_atlas);
    end

    % Time series of the voxels in the atlas (single precision to save memory)
    fMRI = niftiread(info);
    T = size(fMRI,4);
    fMRI = single(reshape(fMRI,[],T));
    fMRI = fMRI(ind_atlas,:);
    if isfield(info,'MultiplicativeScaling') && info.MultiplicativeScaling ~= 0 && info.MultiplicativeScaling ~= 1
        fMRI = fMRI*info.MultiplicativeScaling + info.AdditiveOffset;
    end

    % Voxels with signal: not constant (e.g. outside the brain mask) and
    % not NaN
    valid = std(fMRI,0,2) > 0 & all(isfinite(fMRI),2);

    % Mean time series of the voxels with signal in each area
    data = NaN(n_areas,T);
    for n = 1:n_areas
        vox = label_atlas == n & valid;
        if any(vox)
            data(n,:) = mean(fMRI(vox,:),1);
        end
    end
    data = double(data);
    if any(isnan(data(:,1)))
        disp(['      - Areas without voxels with signal: ' mat2str(find(isnan(data(:,1)))')])
    end

    save(out_file,'data')
end
disp(' ')
