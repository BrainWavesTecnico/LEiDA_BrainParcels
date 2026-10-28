function [mean_P,sem_P] = FracOccup_mean_sem(P,Index_Conditions,Subjects,stats_level,n_Cond)
%
% Mean and standard error of the mean of the fractional occupancy of each
% FC state in each condition, computed over the same observations used in
% the statistical tests: subjects (averaging the scans of each subject in
% a condition) if stats_level is 'subject', or scans if it is 'scan'.
%
% INPUT:
% P                 (n_scans x n_K x max(K)) fractional occupancy
% Index_Conditions  (1 x n_scans) condition of each scan (0: none)
% Subjects          (n_scans x 1) subject ID of each scan
% stats_level       'subject' or 'scan'
% n_Cond            number of conditions
%
% OUTPUT:
% mean_P, sem_P     (n_Cond x n_K x max(K)) mean and SEM in each condition
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

dims = [size(P,2) size(P,3)];
mean_P = zeros([n_Cond dims]);
sem_P = zeros([n_Cond dims]);
for j = 1:n_Cond
    P_j = reshape(P(Index_Conditions == j,:,:),nnz(Index_Conditions == j),[]);
    if strcmpi(stats_level,'subject')
        [~,~,g] = unique(Subjects(Index_Conditions == j));
        P_j = splitapply(@(x) mean(x,1),P_j,g);
    end
    mean_P(j,:,:) = reshape(mean(P_j,1),dims);
    sem_P(j,:,:) = reshape(std(P_j,0,1)/sqrt(size(P_j,1)),dims);
end
