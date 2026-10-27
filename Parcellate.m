function Parcellate(data_dir,file_pattern,save_dir,parcellation,n_areas)
%
% Parcellate fMRI data in MNI 2mm space according to a given atlas from
% ParcelsMNI2mm.mat. The time series of each area is the mean over its
% voxels. One .mat file (variable data, n_areas x time) is saved per scan
% in save_dir, keeping the subfolder structure of data_dir (e.g.
% sub-XXX/ses-YYY/). Scans already parcellated are not computed again.
%
% INPUT:
% data_dir      directory with the 4D NIfTI files (.nii or .nii.gz); files
%               are searched in this directory and all its subfolders
% file_pattern  pattern of the file names to use, e.g. '*bold.nii.gz'
% save_dir      directory to save the parcellated data
% parcellation  atlas in ParcelsMNI2mm.mat (e.g. 'AAL116', 'glasser360')
% n_areas       number of areas of the atlas to consider ([] for all)
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARCELLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ')

N_parcels = Atlas_Nparcels(parcellation);
if isempty(N_parcels)
    error('LEiDA:atlas','To parcellate data in MNI space, Parcellation must be an atlas available in ParcelsMNI2mm.mat');
end
if isempty(n_areas)
    n_areas = N_parcels;
end

% Atlas in MNI 2mm space
S = load('ParcelsMNI2mm.mat',['V_' parcellation]);
V_Parcel = S.(['V_' parcellation]);
sz = size(V_Parcel);

% Files to parcellate
files = dir(fullfile(data_dir,'**',file_pattern));
files = files(~[files.isdir]);
if isempty(files)
    error('LEiDA:nodata','No files matching ''%s'' were found in %s', file_pattern, data_dir);
end
data_root = [regexprep(data_dir,'[\\/]+$','') filesep];

for s = 1:numel(files)

    in_file = fullfile(files(s).folder,files(s).name);
    % Same subfolders as in data_dir
    out_dir = fullfile(save_dir,erase([files(s).folder filesep],data_root));
    base_name = regexprep(files(s).name,'\.nii(\.gz)?$','');
    out_file = fullfile(out_dir,[base_name '_' parcellation '.mat']);
    if isfile(out_file)
        disp(['Already parcellated: ' out_file])
        continue
    end
    if ~exist(out_dir,'dir')
        mkdir(out_dir);
    end

    disp(['Parcellating ' in_file ' using parcellation ' parcellation]);
    fMRI = double(niftiread(in_file));
    T = size(fMRI,4);

    % Check if the data is in the MNI 2mm grid of the atlas
    if ~isequal([size(fMRI,1) size(fMRI,2) size(fMRI,3)],sz)
        disp(['      - Data size ' mat2str(size(fMRI,1:3)) ' differs from the atlas ' mat2str(sz) ...
              ' -> data resized to the atlas grid (check that the data is in MNI space)'])
        fMRI_resized = zeros([sz T]);
        for t = 1:T
            fMRI_resized(:,:,:,t) = imresize3(fMRI(:,:,:,t),sz);
        end
        fMRI = fMRI_resized;
        clear fMRI_resized
    end

    % Mean time series of the voxels in each area
    fMRI = reshape(fMRI,[],T);
    data = zeros(n_areas,T);
    for n = 1:n_areas
        data(n,:) = mean(fMRI(V_Parcel(:) == n,:),1,'omitnan');
    end

    save(out_file,'data')
end
disp(' ')
