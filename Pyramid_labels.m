function Pyramid_labels(L,rangeK,cols_per_state)
%
% Add the labels of the rows (K = ...) on the left and of the FC states
% (c = ...) at the bottom of a pyramid figure.
%
% INPUT:
% L               layout from Pyramid_layout
% rangeK          values of K (one per row)
% cols_per_state  number of grid columns used by each FC state
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

for k = 1:length(rangeK)
    annotation('textbox',[0 L.bottom+(L.n_rows-k)*L.h L.left L.h],'String',['K=' num2str(rangeK(k))],...
               'EdgeColor','none','HorizontalAlignment','right','VerticalAlignment','middle',...
               'FontSize',9,'FontWeight','bold','FitBoxToText','off');
end
for c = 1:rangeK(end)
    annotation('textbox',[L.left+(c-1)*cols_per_state*L.w 0.22*L.bottom cols_per_state*L.w 0.2*L.bottom],...
               'String',['c=' num2str(c)],'EdgeColor','none','HorizontalAlignment','center',...
               'VerticalAlignment','middle','FontSize',9,'FontWeight','bold','FitBoxToText','off');
end
