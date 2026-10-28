function N_parcels = Atlas_Nparcels(parcellation)
%
% Number of parcels of an atlas (included in ParcelsMNI2mm.mat, or a NIfTI
% file; see Load_atlas)
%
% INPUT:
% parcellation  name of the atlas (e.g. 'AAL116', 'Desikan', 'glasser360'),
%               or full path of a NIfTI atlas
%
% OUTPUT:
% N_parcels     number of parcels in the atlas; empty if parcellation is
%               empty or not available
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

N_parcels = [];
if isempty(parcellation)
    return
end

Atlas = Load_atlas(parcellation);
if isempty(Atlas)
    atlases = erase(who('-file','ParcelsMNI2mm.mat'),'V_');
    atlases = atlases(~startsWith(atlases,'label') & ~ismember(atlases,{'Yeo7','Subcortical'}));
    warning('LEiDA:atlas', ['Parcellation ''%s'' is not a NIfTI file nor an atlas of ParcelsMNI2mm.mat ' ...
            '(%s). The FC states will not be rendered in the brain.'], ...
            parcellation, strjoin(atlases',', '));
    return
end
N_parcels = Atlas.N_parcels;
