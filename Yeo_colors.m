function [YeoColor,Yeo_names] = Yeo_colors
%
% Colors and names of the 7 resting-state networks defined by Yeo et al.
% (2011), used in all the LEiDA plots. Change the colors here to change
% them in all plots.
%
% OUTPUT:
% YeoColor      (7 x 3) RGB color of each network (values in [0,1])
% Yeo_names     (1 x 7) name of each network
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

YeoColor = [180 115 208; 76 112 210; 100 174 117; 230 150 200; 207 225 185; 250 221 107; 227 120 91]./256;
Yeo_names = {'Visual','Somatomotor','Dorsal Att.','Ventral Att.','Limbic','Frontoparietal','Default Mode'};
