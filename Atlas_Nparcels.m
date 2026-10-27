function N_parcels = Atlas_Nparcels(parcellation)
%
% Number of parcels of an atlas available in ParcelsMNI2mm.mat
%
% INPUT:
% parcellation  name of the atlas (e.g. 'AAL116', 'Desikan', 'glasser360')
%
% OUTPUT:
% N_parcels     number of parcels in the atlas; empty if parcellation is
%               empty or not available in ParcelsMNI2mm.mat
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

N_parcels = [];
if isempty(parcellation)
    return
end

available = who('-file','ParcelsMNI2mm.mat');
if ~ismember(['V_' parcellation], available)
    atlases = erase(available(startsWith(available,'V_')),'V_');
    warning('LEiDA:atlas', ['Parcellation ''%s'' is not available in ParcelsMNI2mm.mat ' ...
            '(available: %s). The FC states will not be rendered in the brain.'], ...
            parcellation, strjoin(atlases',', '));
    return
end

S = load('ParcelsMNI2mm.mat',['V_' parcellation]);
V = S.(['V_' parcellation]);
N_parcels = double(max(V(:)));
