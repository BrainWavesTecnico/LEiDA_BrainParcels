function V_dst = Resample_labels(V_src,A_src,size_dst,A_dst)
%
% Resample a volume of labels (e.g. a parcellation) to another voxel grid
% with nearest-neighbour interpolation, using the affine transformations
% of both grids (see Nifti_affine), so that different resolutions, fields
% of view and orientations of the same space (e.g. MNI) are handled.
%
% INPUT:
% V_src     (3D) labels in the source grid
% A_src     (4 x 4) affine of the source grid (Matlab indices to mm)
% size_dst  size of the destination grid [nx ny nz]
% A_dst     (4 x 4) affine of the destination grid
%
% OUTPUT:
% V_dst     labels in the destination grid (0 outside the source volume)
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

[I,J,K] = ndgrid(1:size_dst(1),1:size_dst(2),1:size_dst(3));
% Destination voxels -> world (mm) -> source voxels
src = A_src \ (A_dst*[I(:) J(:) K(:) ones(numel(I),1)]');
src = round(src(1:3,:));
inside = all(src >= 1,1) & src(1,:) <= size(V_src,1) & src(2,:) <= size(V_src,2) & src(3,:) <= size(V_src,3);

V_dst = zeros(size_dst,'like',V_src);
V_dst(inside) = V_src(sub2ind(size(V_src),src(1,inside),src(2,inside),src(3,inside)));
