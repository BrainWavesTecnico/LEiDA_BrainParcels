function Pyramid_legend(L,lines,rangeK)
%
% Add the legend of a pyramid figure in its top right corner (the empty
% part of the pyramid), followed by the colors of the significance levels.
%
% INPUT:
% L         layout from Pyramid_layout
% lines     cell array with the lines of text of the legend
% rangeK    range of K (for the significance thresholds)
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

% Escape the TeX special characters of the text (e.g. _ in condition names)
lines = regexprep(cellstr(lines),'([_^{}\\])','\\$1');

% Colors of the significance levels (see Significance_color.m)
sig = @(p,K) sprintf('\\color[rgb]{%g,%g,%g}',Significance_color(p,K,rangeK));
lines{end+1} = ['Significance: ' sig(0.05,Inf) 'p<0.05\color{black}, ' sig(0.05/rangeK(end),rangeK(end)) ...
                'p<0.05/K\color{black}, ' sig(0,1) 'p<0.05/' num2str(sum(rangeK)) '\color{black} (all K and FC states)'];

% Top right corner, over the first rows (empty in the pyramid)
n_rows = min(6,L.n_rows-1);
annotation('textbox',[0.5 1-L.top-n_rows*L.h 0.495 n_rows*L.h],'String',lines,'FitBoxToText','off',...
           'HorizontalAlignment','right','VerticalAlignment','top','EdgeColor','none',...
           'Interpreter','tex','FontSize',11);
