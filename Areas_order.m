function Order = Areas_order(parcellation,n_areas)
%
% Order of the areas in the vector and matrix plots of the centroids: in
% the AAL parcellations, where left and right areas alternate, the left
% areas are shown first and the right areas after them in reverse order,
% so that homologous areas are symmetric. Other parcellations are shown in
% their own order.
%
% INPUT:
% parcellation  name of the parcellation
% n_areas       number of areas
%
% OUTPUT:
% Order         (1 x n_areas) order of the areas
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

if startsWith(lower(parcellation),'aal')
    Order = [1:2:n_areas n_areas-mod(n_areas,2):-2:2];
else
    Order = 1:n_areas;
end
