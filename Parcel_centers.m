function MNI_coord = Parcel_centers(parcellation,n_areas)
%
% Center of gravity of each area of a parcellation in MNI 2mm voxel space,
% as [y x z] voxel coordinates (the axes used in the brain renderings).
%
% INPUT:
% parcellation  name of the parcellation (atlas in ParcelsMNI2mm.mat)
% n_areas       number of areas considered (the first n_areas)
%
% OUTPUT:
% MNI_coord     (n_areas x 3) center of each area
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

S = load('ParcelsMNI2mm.mat',['V_' parcellation]);
Volume = S.(['V_' parcellation]);
MNI_coord = zeros(n_areas,3);
for n = 1:n_areas
    [X,Y,Z] = ind2sub(size(Volume),find(Volume == n));
    MNI_coord(n,:) = mean([Y,X,Z],1);
end
