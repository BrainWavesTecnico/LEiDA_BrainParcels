function L = Pyramid_layout(n_rows,n_cols)
%
% Layout of the pyramid figures (one row per K, one or more columns per FC
% state), with a margin on the left for the row labels (K = ...) and an
% extra space at the bottom for the x tick labels and the column labels
% (c = ...). Use with Pyramid_axes and Pyramid_labels.
%
% INPUT:
% n_rows    number of rows (number of values of K)
% n_cols    number of columns of the grid
%
% OUTPUT:
% L         struct with the layout (normalized figure units)
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

L.n_rows = n_rows;
L.n_cols = n_cols;
L.left = 0.035;                       % space for the K labels
L.right = 0.005;
L.top = 0.005;
L.h = (1 - L.top)/(n_rows + 1.5);     % height of each row
L.bottom = 1.5*L.h;                   % space for x tick labels and c labels
L.w = (1 - L.left - L.right)/n_cols;  % width of each column
