function Atlas = Load_atlas(parcellation)
%
% Load a parcellation (atlas), either one of the atlases included in
% ParcelsMNI2mm.mat (e.g. 'AAL116', 'Desikan', 'Schaefer100', 'glasser360')
% or a NIfTI file with integer labels in MNI space, at any resolution (e.g.
% the atlas used to parcellate the data, as it is delivered).
%
% For a NIfTI atlas, the names of the areas are read, if present, from a
% text file with the same name and extension .txt, .tsv or .csv (e.g.
% atlas.nii.gz -> atlas.txt): one area per line, in the order of the labels,
% either 'name' or 'label name' (the first column is ignored if it is the
% label number).
%
% INPUT:
% parcellation  name of an atlas in ParcelsMNI2mm.mat, or full path of a
%               NIfTI file (.nii or .nii.gz)
%
% OUTPUT:
% Atlas         struct with:
%   V           labels in the MNI 2mm grid of the toolbox (91 x 109 x 91,
%               as MNI152_T1_2mm_brain_mask.nii), used for the brain figures
%   V_native    labels in the grid of the atlas
%   A_native    affine of the atlas grid (Matlab indices to mm)
%   N_parcels   largest label
%   labels      names of the areas ({} if not available)
%   name        short name of the atlas
% Empty ([]) if the atlas is not found.
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

persistent cache
if isempty(cache)
    cache = containers.Map();
end
Atlas = [];
if isempty(parcellation)
    return
end
if isKey(cache,parcellation)
    Atlas = cache(parcellation);
    return
end

% MNI 2mm grid of the toolbox
info_MNI = niftiinfo('MNI152_T1_2mm_brain_mask.nii');
A_MNI = Nifti_affine(info_MNI);
size_MNI = info_MNI.ImageSize(1:3);

if endsWith(lower(parcellation),{'.nii','.nii.gz'})
    % Atlas as a NIfTI file
    if ~isfile(parcellation)
        warning('LEiDA:atlas','Atlas file %s not found',parcellation);
        return
    end
    info = niftiinfo(parcellation);
    Atlas.V_native = round(double(niftiread(info)));
    Atlas.V_native = Atlas.V_native(:,:,:,1);
    Atlas.A_native = Nifti_affine(info);
    Atlas.V = Resample_labels(Atlas.V_native,Atlas.A_native,size_MNI,A_MNI);
    [~,Atlas.name] = fileparts(regexprep(parcellation,'\.gz$',''));
    Atlas.labels = Read_labels(regexprep(parcellation,'\.nii(\.gz)?$',''));
else
    % Atlas included in ParcelsMNI2mm.mat (MNI 2mm grid)
    if ~ismember(['V_' parcellation],who('-file','ParcelsMNI2mm.mat'))
        return
    end
    S = load('ParcelsMNI2mm.mat',['V_' parcellation]);
    Atlas.V = double(S.(['V_' parcellation]));
    Atlas.V_native = Atlas.V;
    Atlas.A_native = A_MNI;
    Atlas.name = parcellation;
    Atlas.labels = {};
end
Atlas.N_parcels = max(Atlas.V_native(:));
cache(parcellation) = Atlas;
end


function labels = Read_labels(base)
% Names of the areas from base.txt, base.tsv or base.csv, if present
labels = {};
for ext = {'.txt','.tsv','.csv'}
    if isfile([base ext{1}])
        lines = cellstr(strtrim(readlines([base ext{1}])));
        lines = lines(~cellfun(@isempty,lines));
        % Remove the label number in the first column, if present
        parts = regexp(lines,'^\s*(\d+)[\s,;]+(.*)$','tokens','once');
        has_number = ~cellfun(@isempty,parts);
        lines(has_number) = cellfun(@(t) t{2},parts(has_number),'UniformOutput',false);
        % Remove a header line (without label number), if present
        if numel(lines) > 1 && ~has_number(1) && all(has_number(2:end))
            lines = lines(2:end);
        end
        labels = strtrim(lines(:));
        return
    end
end
end
