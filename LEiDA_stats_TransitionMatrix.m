function LEiDA_stats_TransitionMatrix(data_dir,save_dir,selectedK)
%
% For the selected K compute the transition probability matrix for each
% participant. Perform hypothesis tests to check for differences in the
% state-to-state transition probabilities between conditions.
%
% INPUT:
% data_dir       directory where the LEiDA results are saved
% save_dir       directory to save the results from the hypothesis tests
%                of the comparison of the mean transition probabilities for
%                the selected K
% selectedK      K defined by the user
%
% OUTPUT:
% TM             transition matrix for optimal K for all participants
% TMnorm         transition probability matrix (TPM) for optimal K for all
%                participants
% TM_pval2sided  two-sided p-value from testing whether the mean transition
%                probability of a given FC state differs between conditions      
% effectsize    absolute Hedges' effect size
%
% The statistical options (paired/unpaired, number of permutations and
% bootstraps, subject or scan level) are the ones used in
% LEiDA_stats_FracOccup (see LEiDA_compare_conditions.m)
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

% File with leading eigenvectors (output from LEiDA_data.m)
file_V1 = 'LEiDA_EigenVectors.mat';
% File with the Kmeans results (output from LEiDA_cluster.m)
file_cluster = 'LEiDA_Clusters.mat';
% File with results for the fractional occupancy (output from LEiDA_stats_FracOccup.m)
file_P = 'LEiDA_Stats_FracOccup.mat';

% Load required data:
if isfile([data_dir file_V1])
    load([data_dir file_V1], 'Time_sessions', 'idx_data', 'Scan_info');
end
if isfile([data_dir file_cluster])
    load([data_dir file_cluster], 'Kmeans_results', 'rangeK');
end
if isfile([data_dir file_P])
    load([data_dir file_P], 'cond', 'P', 'Index_Conditions', 'pair', ...
         'stats_level', 'n_permutations', 'n_bootstraps');
end
Subjects = Scan_info.subject(idx_data);

% Number of scans considered to compute V1
N_scans = length(idx_data);

% Number of conditions of the experiment
n_Cond = size(cond,2);

% Code below computes the fractional occupancy and runs hypothesis tests
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSITION PROBABILITIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ')

%% CALCULATE THE TRANSITION MATRIX FOR THE OPTIMAL K

% Matrix with probability of occurrence of each FC state for chosen K
P_K = squeeze(P(:,rangeK == selectedK,1:selectedK));

% Transition probability matrix
TM = zeros(N_scans, selectedK, selectedK);
TMnorm = zeros(N_scans, selectedK, selectedK);

for s = 1:N_scans

    T = Time_sessions == idx_data(s);
    Ctime = Kmeans_results{rangeK == selectedK}.IDX(T);
    % alpha corresponds to departure state;
    for alpha = 1:selectedK
        % beta corresponds to arrival state;
        for beta = 1:selectedK
            % count number transition from alpha to beta
            alpha2beta = 0;
            for i = 1:length(Ctime)-1
                % if alpha = beta then we have 1 transition
                if isequal(Ctime(i),alpha) && isequal(Ctime(i+1),beta)
                    alpha2beta = alpha2beta + 1;
                end
            end
            TM(s,alpha,beta) = alpha2beta;
        end
    end
    % normalise by T-1 as in Vohryzek et al., (2020)
    TM(s,:,:) = TM(s,:,:)./(size(Ctime,2)-1);
    TMnorm(s,:,:) = squeeze(TM(s,:,:))./P_K(s,:)';
end

% Set all Nan values to zero
% If a FC state has probability of occurence 0 then there should be a
% probability of 0 to transition to that state
% numNaN = sum(isnan(TMnorm),'all');
% TMnorm(isnan(TMnorm)) = 0;

clear Kmeans_results

%% PERMUTATION TESTS

disp('Testing differences in transition probabilities between conditions:')
results = LEiDA_compare_conditions(reshape(TMnorm,N_scans,[]),Index_Conditions,Subjects,cond,pair,...
                                   stats_level,n_permutations,n_bootstraps);

% Back to the format (condition pair x departure state x arrival state)
dims = [size(results.pval,1) selectedK selectedK];
TM_pval = reshape(results.pval,dims);
TM_pval2sided = reshape(results.pval2sided,dims);
effectsize = reshape(results.effectsize,dims);
effectsize_signed = reshape(results.effectsize_signed,dims);
effectsize_ci = reshape(results.effectsize_ci,[dims 2]);

% Name of the file to save output
save_file = 'LEiDA_Stats_TransitionMatrix.mat';

save(fullfile(save_dir,save_file),'TM','TMnorm','TM_pval','TM_pval2sided','effectsize',...
     'effectsize_signed','effectsize_ci','cond','rangeK','file_V1','file_cluster','file_P','Index_Conditions')

disp(' ')
disp(['Transition probability values and results from permutation tests saved successfully as ' save_file])
disp(' ')