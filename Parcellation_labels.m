function labels = Parcellation_labels(parcellation,n_areas)
%
% Names of the areas of a parcellation (from ParcelsMNI2mm.mat, or from the
% text file given with a NIfTI atlas, see Load_atlas). For
% parcellations without names, the areas are named 'Area 1', 'Area 2', ...
%
% INPUT:
% parcellation  name of the parcellation (e.g. 'AAL116', 'Desikan')
% n_areas       number of areas considered (the first n_areas)
%
% OUTPUT:
% labels        (n_areas x 1) cell array with the name of each area
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

labels = arrayfun(@(n) sprintf('Area %d',n),(1:n_areas)','UniformOutput',false);

% Names given with a NIfTI atlas (see Load_atlas)
Atlas = Load_atlas(parcellation);
if ~isempty(Atlas) && ~isempty(Atlas.labels)
    n = min(n_areas,numel(Atlas.labels));
    labels(1:n) = Atlas.labels(1:n);
    return
end

switch lower(parcellation)
    case 'aal116'
        name = 'label116';
    case 'aal120'
        name = 'label120';
    case 'schaefer100'
        name = 'labelSchaefer100P7Networks';
    case 'desikan'
        name = 'labelDesikan70';
    otherwise
        return
end
S = load('ParcelsMNI2mm.mat',name);
names = strtrim(cellstr(string(S.(name))));
n = min(n_areas,numel(names));
labels(1:n) = names(1:n);
