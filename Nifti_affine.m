function A = Nifti_affine(info)
%
% Affine transformation of a NIfTI image from Matlab voxel indices (i,j,k,
% starting at 1, as returned by niftiread) to world (mm) coordinates:
% [x y z 1]' = A * [i j k 1]'.
%
% INPUT:
% info      header of the image (niftiinfo)
%
% OUTPUT:
% A         (4 x 4) affine transformation
%
% The sform is used if present (sform_code > 0), otherwise the qform.
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

if isfield(info,'raw') && isfield(info.raw,'sform_code') && info.raw.sform_code > 0
    % NIfTI sform, which applies to voxel indices starting at 0
    M = double([info.raw.srow_x(:)'; info.raw.srow_y(:)'; info.raw.srow_z(:)'; 0 0 0 1]);
else
    % qform (or scaling only), as computed by niftiinfo (row-vector form)
    M = double(info.Transform.T');
end
% From indices starting at 1 to indices starting at 0
A = M*[eye(3) -ones(3,1); 0 0 0 1];
