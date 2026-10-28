function Draw_links_in_cortex(Vc, MNI_coord, scortex, n_links)
%
% Plot the glass cortex, the n_links strongest links between regions with
% positive value in the centroid vector and the nodes colored by their
% contribution, in the current axes. The view and camlight are set by the
% caller.
%
% INPUT:
% Vc         centroid vector (scaled by its largest absolute value)
% MNI_coord  (n_areas x 3) center of each area in MNI 2mm voxel space
% scortex    smoothed brain mask
% n_links    number of links to plot
%
% Authors: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%          Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

hold on

% Cortex patch
cortexpatch = patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.1);
reducepatch(cortexpatch,0.5);
isonormals(scortex,cortexpatch);

% Plot the strongest links between regions with positive value in centroid
% vector, where the strength of the link between areas n and p is Vc(n)*Vc(p)
n_strong = find(Vc > 0);

if numel(n_strong) > 1
    pairs = nchoosek(n_strong,2);
    link_strength = Vc(pairs(:,1)).*Vc(pairs(:,2));
    [~, order] = sort(link_strength,'descend');
    pairs = pairs(order(1:min(n_links,numel(order))),:);
    for l = 1:size(pairs,1)
        c1 = MNI_coord(pairs(l,1),:);
        c2 = MNI_coord(pairs(l,2),:);
        plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)],'Color',[0.8 0 0],'LineWidth',1);
    end
end

a = 2;
[x,y,z] = sphere;
x = a*x;
y = a*y;
z = a*z;

% Scale positive and negative values separately so that the most positive
% node is pure red and the most negative node is pure blue
max_pos = max([Vc(Vc > 0) eps]);
max_neg = max([-Vc(Vc < 0) eps]);

for n = 1:length(Vc)
    if Vc(n) > 0 % color nodes with positive contribution in a gradient of red colors
        surf(x+MNI_coord(n,1), y+MNI_coord(n,2),z+MNI_coord(n,3),'FaceColor',[Vc(n)/max_pos 0 0],'EdgeColor','none','FaceAlpha',0.7);
    elseif Vc(n) < 0 % color nodes with negative contribution in a gradient of blue colors
        surf(x+MNI_coord(n,1), y+MNI_coord(n,2),z+MNI_coord(n,3),'FaceColor',[0 0 -Vc(n)/max_neg],'EdgeColor','none','FaceAlpha',0.7);
    end
end

axis off;
axis equal;

material dull;
lighting gouraud;
daspect([1 1 1])
xlim([5 105]); ylim([5 85]); zlim([0 80])

hold off

end
