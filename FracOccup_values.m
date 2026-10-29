function [values,subjects] = FracOccup_values(P,Index_Conditions,Subjects,stats_level,n_Cond,k,c)
%
% Fractional occupancy of one FC state in each condition, for the same
% observations used in the statistical tests: one value per subject
% (average of the scans of the subject in the condition) if stats_level is
% 'subject', or one value per scan if it is 'scan'.
%
% INPUT:
% P                 (n_scans x n_K x max(K)) fractional occupancy
% Index_Conditions  (1 x n_scans) condition of each scan (0: none)
% Subjects          (n_scans x 1) subject ID of each scan
% stats_level       'subject' or 'scan'
% n_Cond            number of conditions
% k, c              index of K in rangeK and FC state
%
% OUTPUT:
% values            (1 x n_Cond) cell with the values of each condition
% subjects          (1 x n_Cond) cell with the subject of each value
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

values = cell(1,n_Cond);
subjects = cell(1,n_Cond);
for j = 1:n_Cond
    in_cond = Index_Conditions == j;
    v = P(in_cond,k,c);
    s = string(Subjects(in_cond));
    if strcmpi(stats_level,'subject')
        [s,~,g] = unique(s,'stable');
        v = accumarray(g(:),v(:),[],@mean);
    end
    values{j} = v(:);
    subjects{j} = s(:);
end
