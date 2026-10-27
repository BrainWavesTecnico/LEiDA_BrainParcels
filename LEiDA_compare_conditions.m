function results = LEiDA_compare_conditions(X,Index_Conditions,subjects,cond,pair,stats_level,n_permutations,n_bootstraps)
%
% Compare a set of measures (e.g. fractional occupancy of each FC state)
% between each pair of conditions using permutation tests
% (see permutation_test.m).
%
% INPUT:
% X                 (n_scans x M) value of M measures in each scan
% Index_Conditions  (1 x n_scans) condition of each scan (0: none)
% subjects          (n_scans x 1) subject ID of each scan
% cond              (1 x n_Cond) name of each condition
% pair              0, unpaired samples; 1, paired samples (the scans of
%                   the same subject in both conditions are paired)
% stats_level       'subject': the values of all the scans of a subject in
%                   a condition are averaged, so that each subject counts
%                   once (recommended when subjects have several scans);
%                   'scan': each scan is an independent observation
% n_permutations    number of permutations
% n_bootstraps      number of bootstrap samples for the confidence interval
%                   of the effect size (0 to skip)
%
% OUTPUT:
% results.pval          (n_pairs x M) smallest one-sided p-value
% results.pval2sided    (n_pairs x M) two-sided p-value
% results.tstat         (n_pairs x M) t statistic (> 0: 1st condition higher)
% results.effectsize    (n_pairs x M) absolute effect size (Hedges' g)
% results.effectsize_signed  (n_pairs x M) signed effect size
% results.effectsize_ci (n_pairs x M x 2) 95% bootstrap CI of the signed
%                       effect size
% results.n             (n_pairs x 2) number of observations in each
%                       condition (subjects/scans or pairs)
% results.cond_pairs    (n_pairs x 2) conditions compared in each pair
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

n_Cond = numel(cond);
n_pairs = n_Cond*(n_Cond-1)/2;
M = size(X,2);
subjects = string(subjects(:));

results.pval = ones(n_pairs,M);
results.pval2sided = ones(n_pairs,M);
results.tstat = zeros(n_pairs,M);
results.effectsize = zeros(n_pairs,M);
results.effectsize_signed = zeros(n_pairs,M);
results.effectsize_ci = NaN(n_pairs,M,2);
results.n = zeros(n_pairs,2);
results.cond_pairs = zeros(n_pairs,2);

cond_pair = 1;
for cond1 = 1:n_Cond-1
    for cond2 = cond1+1:n_Cond

        [A,subj1] = Select_condition(X,Index_Conditions == cond1,subjects,stats_level);
        [B,subj2] = Select_condition(X,Index_Conditions == cond2,subjects,stats_level);

        if pair
            % Pair the observations of the same subject in both conditions
            if numel(unique(subj1)) < numel(subj1) || numel(unique(subj2)) < numel(subj2)
                error('LEiDA:stats',['Some subjects have more than one scan in ''%s'' or ''%s''; ' ...
                      'use Stats_level = ''subject'' for paired tests'], cond{cond1}, cond{cond2});
            end
            [~,i1,i2] = intersect(subj1,subj2,'stable');
            n_out = numel(subj1) + numel(subj2) - 2*numel(i1);
            A = A(i1,:);
            B = B(i2,:);
            msg = [num2str(numel(i1)) ' paired subjects'];
            if n_out
                msg = [msg ' (' num2str(n_out) ' subjects without both conditions excluded)']; %#ok<AGROW>
            end
        else
            msg = [num2str(size(A,1)) ' vs ' num2str(size(B,1)) ' ' stats_level 's'];
            if any(ismember(subj1,subj2))
                msg = [msg ' (Attention: some subjects are in both conditions; ' ...
                       'observations are not independent, consider paired tests)']; %#ok<AGROW>
            end
        end
        disp(['- ' cond{cond1} ' vs ' cond{cond2} ': ' msg])

        stats = permutation_test(A,B,pair,n_permutations,n_bootstraps);

        results.pval(cond_pair,:) = min(stats.pvals,[],1);
        results.pval2sided(cond_pair,:) = stats.pvals_2sided;
        results.tstat(cond_pair,:) = stats.tvals;
        results.effectsize(cond_pair,:) = abs(stats.eff);
        results.effectsize_signed(cond_pair,:) = stats.eff;
        results.effectsize_ci(cond_pair,:,:) = permute(stats.eff_ci,[3 2 1]);
        results.n(cond_pair,:) = stats.n;
        results.cond_pairs(cond_pair,:) = [cond1 cond2];
        cond_pair = cond_pair + 1;
    end
end
end


function [Xc,subj] = Select_condition(X,in_cond,subjects,stats_level)
% Values of the scans in one condition, averaged per subject if
% stats_level is 'subject'
switch lower(stats_level)
    case 'subject'
        [subj,~,g] = unique(subjects(in_cond),'stable');
        Xc = splitapply(@(x) mean(x,1,'omitnan'),X(in_cond,:),g);
    case 'scan'
        subj = subjects(in_cond);
        Xc = X(in_cond,:);
    otherwise
        error('LEiDA:stats','Stats_level must be ''subject'' or ''scan''');
end
end
