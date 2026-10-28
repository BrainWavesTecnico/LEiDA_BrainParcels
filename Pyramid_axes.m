function ax = Pyramid_axes(L,row,col,n_span,gap)
%
% Create the axes in the cell (row, col) of a pyramid layout.
%
% INPUT:
% L         layout from Pyramid_layout
% row, col  row and column of the cell (row 1 at the top)
% n_span    number of columns spanned by the axes (default 1)
% gap       [horizontal vertical] empty space around the axes, as a
%           fraction of the cell width and height (default [0.1 0.1])
%
% OUTPUT:
% ax        axes handle
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

if nargin < 4 || isempty(n_span)
    n_span = 1;
end
if nargin < 5 || isempty(gap)
    gap = [0.1 0.1];
end
gap = gap.*[1 1];   % same gap in both directions if a scalar is given

x = L.left + (col-1)*L.w + gap(1)*L.w/2;
y = L.bottom + (L.n_rows-row)*L.h + gap(2)*L.h/2;
ax = axes('Position',[x y n_span*L.w-gap(1)*L.w L.h*(1-gap(2))]);
