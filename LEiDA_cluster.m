function [Kmeans_results,rangeK] = LEiDA_cluster(data_dir,rangeK)
%
% Cluster all leading eigenvectors into a set of K clusters
% using the K-means algorithm. 
% This function returns an optimal solution for each value of K.
%
% INPUT:
% data_dir         directory where the LEiDA_EigenVectors file is saved
% rangeK           range of values of K to compute (default 2:20)
%
% OUTPUT:
% Kmeans_results   clustering solutions for the range of values K
% rangeK           range of values of K considered for clustering
%
% Author: Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLUSTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

% Default range of the number of clusters
if nargin < 2 || isempty(rangeK)
    rangeK = 2:20;
end

% File with leading eigenvectors (output from LEiDA_data.m)
file_V1 = 'LEiDA_EigenVectors.mat';

% Load the leading eigenvectors
load([data_dir file_V1], 'V1_all');

SignV1=sum(V1_all,2);

V1_all(SignV1>0,:)=-V1_all(SignV1>0,:);

% Cell to store results for each K using the K-means algorithm
Kmeans_results = cell(size(rangeK));

% Run the replicates in parallel if the Parallel Computing Toolbox is available
use_parallel = ~isempty(ver('parallel')) && license('test','Distrib_Computing_Toolbox');

% Number of times to run the clustering
% (it is convenient to run more replicates for larger samples)
replicates = 25; %ceil(20 + size(V1_all,1)/400);

disp(' ');
disp('Clustering eigenvectors into:')

% For each K do:
for K = 1:length(rangeK)
    disp(['- ' num2str(rangeK(K)) ' FC states'])
    [IDX, C, SUMD, D] = kmeans(V1_all,rangeK(K),'Distance','Cosine','Replicates',replicates,'MaxIter',500,'OnlinePhase','off',...,
                              'Display','off','Options',statset('UseParallel',use_parallel));
    
                         
    % ind_sort sorts the clusters in descending order of occupancy
    [~, ind_sort] = sort(hist(IDX,1:rangeK(K)),'descend');
    [~, idx_sort] = sort(ind_sort,'ascend');
    Kmeans_results{K}.IDX = idx_sort(IDX);     % Cluster time course - numeric collumn vector
    Kmeans_results{K}.C = C(ind_sort,:);       % Cluster centroids (FC patterns)
    Kmeans_results{K}.SUMD = SUMD(ind_sort);   % Within-cluster sums of point-to-centroid distances
    Kmeans_results{K}.D = D(:,ind_sort);       % Distance from each point to every centroid
end

% Name of the file to save output
save_file = 'LEiDA_Clusters.mat';

save(fullfile(data_dir,save_file), 'Kmeans_results', 'rangeK')
disp(['K-means clustering completed and results saved as ' save_file])
disp(' ');