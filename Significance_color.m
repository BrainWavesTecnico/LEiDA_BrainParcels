function color = Significance_color(pval,K,rangeK)
%
% Color of a p-value according to the significance levels used in all the
% LEiDA figures.
%
% INPUT:
% pval      p-value
% K         number of FC states of the clustering solution
% rangeK    range of K tested (for the correction over all FC states)
%
% OUTPUT:
% color     dark blue  p <= 0.05/sum(rangeK)
%           dark green p <= 0.05/K
%           dark red   p <= 0.05
%           black      otherwise
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

if pval <= 0.05/sum(rangeK)
    color = [0.07, 0.62, 1.00];   % dark blue
elseif pval <= 0.05/K
    color = [0.3, 0.5, 0];        % dark green
elseif pval <= 0.05
    color = [0.6, 0, 0];          % dark red
else
    color = [0, 0, 0];            % black
end
