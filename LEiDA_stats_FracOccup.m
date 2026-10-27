function LEiDA_stats_FracOccup(data_dir,Index_Conditions,cond,pair,n_permutations,n_bootstraps,stats_level)
%
% For each clustering solution compute the fractional occupancy of each FC
% state for all scans. Perform permutation tests to check for differences
% in the probability of occurrence of FC states between conditions
% (see LEiDA_compare_conditions.m and permutation_test.m).
%
% INPUT:
% data_dir          directory where the LEiDA results are saved
% Index_Conditions  (1 x n_scans) condition of each scan (0: none)
% cond              tags of each condition considered in the experiment
% pair              0, unpaired samples (Welch's t); 1, paired samples
%                   (same subjects across conditions)
% n_permutations    number of permutations
% n_bootstraps      number of bootstrap samples for the 95% confidence
%                   interval of the effect size (0 to skip)
% stats_level       'subject' (default): average the scans of each subject
%                   in each condition; 'scan': each scan is independent
%
% OUTPUT (saved in LEiDA_Stats_FracOccup.mat):
% P             (n_scans x n_K x max(K)) fractional occupancy of each FC
%               state in each clustering solution for all scans
% P_pval        smallest one-sided p-value
% P_pval2sided  two-sided p-value from testing whether the mean fractional
%               occupancy of a given FC state differs between conditions
% P_tstat       t statistic (> 0: first condition of the pair higher)
% effectsize    absolute Hedges' effect size
% effectsize_signed, effectsize_ci  signed effect size and its 95%
%               bootstrap confidence interval
%
% Author: Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

if nargin < 7 || isempty(stats_level)
    stats_level = 'subject';
end

% File with leading eigenvectors (output from LEiDA_data.m)
file_V1 = 'LEiDA_EigenVectors.mat';
% File with K-means results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';

% Load required data:
load([data_dir file_V1], 'Time_sessions', 'idx_data', 'Scan_info');
load([data_dir file_cluster], 'Kmeans_results', 'rangeK');

% Condition and subject of each scan used to compute V1
Index_Conditions = Index_Conditions(idx_data);
Subjects = Scan_info.subject(idx_data);

% Number of scans considered to compute V1
N_scans = length(idx_data);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FRACTIONAL OCCUPANCY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ')
%% CALCULATE THE FRACTIONAL OCCUPANCY OF EACH STATE FOR K-MEANS CLUSTERING SOLUTIONS

% For every fMRI scan calculate probability of each state c.
P = zeros(N_scans,length(rangeK),rangeK(end));

for k = 1:length(rangeK)
    for s = 1:N_scans

        T = Time_sessions == idx_data(s);
        % Vector of the cluster to which each observation belongs to for scan s
        Ctime = Kmeans_results{k}.IDX(T);

        for c = 1:rangeK(k)
            % Probability of occurrence
            P(s,k,c) = mean(Ctime == c);
        end
    end
end

clear Kmeans_results

%% PERMUTATION TESTS

if pair
    test_name = 'paired t';
else
    test_name = 'Welch''s t';
end
disp(['Testing differences in fractional occupancy between conditions (' test_name ', ' ...
      num2str(n_permutations) ' permutations, observations: ' stats_level 's):'])

results = LEiDA_compare_conditions(reshape(P,N_scans,[]),Index_Conditions,Subjects,cond,pair,...
                                   stats_level,n_permutations,n_bootstraps);

% Back to the format (condition pair x K x FC state)
dims = [size(results.pval,1) length(rangeK) rangeK(end)];
P_pval = reshape(results.pval,dims);
P_pval2sided = reshape(results.pval2sided,dims);
P_tstat = reshape(results.tstat,dims);
effectsize = reshape(results.effectsize,dims);
effectsize_signed = reshape(results.effectsize_signed,dims);
effectsize_ci = reshape(results.effectsize_ci,[dims 2]);
n_obs = results.n;

% Name of the file to save output
save_file = 'LEiDA_Stats_FracOccup.mat';

% Save K-means clustering solutions results:
save(fullfile(data_dir,save_file),'P','P_pval','P_pval2sided','P_tstat','effectsize',...
     'effectsize_signed','effectsize_ci','n_obs','cond','rangeK','file_cluster','file_V1',...
     'Index_Conditions','Subjects','pair','stats_level','n_permutations','n_bootstraps')

disp(' ')
disp(['Fractional occupancy values and results from permutation tests saved successfully as ' save_file])
disp(' ')
