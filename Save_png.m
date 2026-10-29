function Save_png(Fig,file)
%
% Save a figure as a .png image at 600 dpi.
%
% INPUT:
% Fig       figure handle
% file      full path of the .png file
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

if exist('exportgraphics','file')
    exportgraphics(Fig,file,'Resolution',600);
else
    print(Fig,file,'-dpng','-r600');
end
