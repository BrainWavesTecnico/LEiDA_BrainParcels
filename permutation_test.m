function stats = permutation_test(A,B,paired,n_permutations,n_bootstraps)
%
% Permutation test comparing the means of two conditions, for paired or
% unpaired samples, with bootstrap confidence intervals of the effect size.
% Several measures (columns) are tested at once; the same permutations are
% used for all measures. NaN values are ignored.
%
% Unpaired samples (paired = 0):
%   Statistic: Welch's t (does not assume equal variances, so no test of
%   equality of variances is needed). Null distribution: the condition
%   labels are randomly permuted across all observations.
%   Effect size: Hedges' g (pooled standard deviation, small-sample
%   corrected).
%
% Paired samples (paired = 1):
%   A(i,:) and B(i,:) must belong to the same subject.
%   Statistic: paired t (t of the differences A-B). Null distribution: the
%   condition labels are randomly swapped within each subject, i.e. the
%   sign of each difference is randomly flipped.
%   Effect size: Hedges' g_z (mean difference / std of the differences,
%   small-sample corrected).
%
% The p-values are computed directly from the permutation distribution
% (no fitting or smoothing of the null distribution), as
% p = (1 + number of permutations at least as extreme) / (1 + n_permutations)
% so they are never 0 and their precision is limited by n_permutations
% (smallest possible p-value = 1/(n_permutations+1)).
%
% INPUT:
% A               (n1 x M) values of M measures in condition 1
% B               (n2 x M) values of M measures in condition 2
% paired          0, unpaired samples; 1, paired samples (n1 == n2)
% n_permutations  number of permutations (e.g. 10000)
% n_bootstraps    number of bootstrap samples for the 95% confidence
%                 interval of the effect size (0 to skip)
%
% OUTPUT:
% stats.tvals         (1 x M) t statistic; positive if condition 1 > 2
% stats.pvals         (2 x M) one-sided p-values: 1st row, condition 1 > 2;
%                     2nd row, condition 1 < 2
% stats.pvals_2sided  (1 x M) two-sided p-values
% stats.eff           (1 x M) effect size (Hedges' g or g_z); positive if
%                     condition 1 > 2
% stats.eff_ci        (2 x M) 95% bootstrap confidence interval of eff
%                     (percentile method); NaN if n_bootstraps = 0
% stats.diffs         (1 x M) difference between the means
% stats.n             number of observations [n1 n2] (pairs if paired)
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

if paired && size(A,1) ~= size(B,1)
    error('LEiDA:stats','Paired samples must have the same number of observations');
end

% Tolerance for comparing permuted and observed statistics (ties)
tol = 1e-12;

if paired
    D = A - B;
    stats.tvals = paired_t(D);
    stats.diffs = mean(D,1,'omitnan');
    stats.n = [size(D,1) size(D,1)];
    n = size(D,1);

    count = zeros(3,size(D,2));
    for iter = 1:n_permutations
        signs = 2*(rand(n,1) > 0.5) - 1;
        count = count + exceed(paired_t(D.*signs), stats.tvals, tol);
    end
else
    X = [A; B];
    n1 = size(A,1);
    N = size(X,1);
    stats.tvals = welch_t(A,B);
    stats.diffs = mean(A,1,'omitnan') - mean(B,1,'omitnan');
    stats.n = [n1 size(B,1)];

    count = zeros(3,size(X,2));
    for iter = 1:n_permutations
        Xp = X(randperm(N),:);
        count = count + exceed(welch_t(Xp(1:n1,:),Xp(n1+1:end,:)), stats.tvals, tol);
    end
end

pvals = (count + 1)./(n_permutations + 1);
stats.pvals = pvals(1:2,:);
stats.pvals_2sided = pvals(3,:);
% Measures with no variability at all (e.g. a state never visited) are not
% different between conditions
no_test = isnan(stats.tvals);
stats.pvals(:,no_test) = 1;
stats.pvals_2sided(no_test) = 1;

% Effect size and bootstrap confidence interval
if paired
    stats.eff = hedges_gz(D);
else
    stats.eff = hedges_g(A,B);
end
stats.eff_ci = NaN(2,numel(stats.eff));
if n_bootstraps > 1
    eff_boot = zeros(n_bootstraps,numel(stats.eff));
    for b = 1:n_bootstraps
        if paired
            eff_boot(b,:) = hedges_gz(D(randi(n,n,1),:));
        else
            eff_boot(b,:) = hedges_g(A(randi(size(A,1),size(A,1),1),:),B(randi(size(B,1),size(B,1),1),:));
        end
    end
    stats.eff_ci = prctile(eff_boot,[2.5 97.5],1);
end
end


function c = exceed(t_perm,t_obs,tol)
% Count permuted statistics at least as extreme as the observed ones:
% rows: right tail, left tail, two-sided
c = [t_perm >= t_obs - tol; t_perm <= t_obs + tol; abs(t_perm) >= abs(t_obs) - tol];
end


function t = welch_t(A,B)
% Welch's t statistic for each column, ignoring NaN
[m1,v1,n1] = mean_var(A);
[m2,v2,n2] = mean_var(B);
t = safe_ratio(m1 - m2, sqrt(v1./n1 + v2./n2));
end


function t = paired_t(D)
% Paired t statistic (one-sample t of the differences), ignoring NaN
[m,v,n] = mean_var(D);
t = safe_ratio(m, sqrt(v./n));
end


function g = hedges_g(A,B)
% Hedges' g for two independent samples
[m1,v1,n1] = mean_var(A);
[m2,v2,n2] = mean_var(B);
s_pooled = sqrt(((n1-1).*v1 + (n2-1).*v2)./(n1+n2-2));
J = 1 - 3./(4*(n1+n2) - 9);
g = J.*safe_ratio(m1 - m2, s_pooled);
g(isnan(g)) = 0;
end


function g = hedges_gz(D)
% Hedges' g_z for paired samples (differences D)
[m,v,n] = mean_var(D);
J = 1 - 3./(4*(n-1) - 1);
g = J.*safe_ratio(m, sqrt(v));
g(isnan(g)) = 0;
end


function [m,v,n] = mean_var(X)
% Mean, variance and number of non-NaN values of each column
n = sum(~isnan(X),1);
m = mean(X,1,'omitnan');
v = var(X,0,1,'omitnan');
end


function r = safe_ratio(num,den)
% num./den, with 0/0 = NaN (no variability) and x/0 = +-Inf
r = num./den;
r(den == 0 & num == 0) = NaN;
end
