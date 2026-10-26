function LEiDA_Start

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function serves to run LEiDA on any given dataset. 
%
% Since we do not know a priori the optimal number of FC states that 
% differentiante between conditions, this function analyzes the solutions 
% obtained for different numbers of FC states (K).
%
% After analyszing the output figures saved in folder LEiDA_Results, 
% the user can then choose the optimal number of FC states for subsequent 
% analysis using the functions LEiDA_AnalysysK.m and LEiDA_AnalysisCentroid.m
% 
% This function contains three sections: 
%      (A) User defines the parameters and properties of the study.
%      (B) Run LEiDA and statistics.
%      (C) Generate and save figures.
%
% Start by reading the README.md file.
%
% A: User input parameters
% B: Run Leading Eigenvector Dynamics Analysis:
%    - Compute the leading eigenvectors for all participants
%    - Cluster the leading eigenvectors of all participants
%    - Compute statistics to compare across conditions.
% C: Figures
%    - Analysis of Fractional Occupancy values
%    - Analysis of Dwell Time values
%    - Pyramid of FC states
%
% Tutorial: README.md
% Version:  V1.2, June 2023
% Authors:  Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

%% A: STUDY PARAMETERS

% Directory of the LEiDA toolbox folder:
LEiDA_directory = '/Users/user/Documents/Research/CognitiveDecline/DONAU/LEiDA_Toolbox_May2025/';
% Directory of the folder with the parcellated neuroimaging data:
Data_directory = '/Users/user/Documents/Research/CognitiveDecline/DONAU/parcellated_data/';
% Give a name of the run to be used to create the folder to save the results:



run_names = {'DONAU_1_Schaefer_100_cortical';
             'DONAU_2_Schaefer_119_cortical_subcortical';
             'DONAU_3_Desikan_70_cortical';
             'DONAU_4__Desikan_89_cortical_subcortical';
             'DONAU_5_Glasser_360_cortical';
             'DONAU_6_Glasser_379_cortical_subcortical'};

for run=6

run_name = run_names{run};

if run == 1 || run == 2
    Parcellation = 'Schaefer100';
    N_areas=100;
elseif run == 3 || run == 4
    Parcellation = 'Desikan';
    N_areas=70;
elseif run == 5 || run == 6
    Parcellation = 'Glasser';
    N_areas=360;
end

if run == 2 || run == 4 || run == 6
    Subcortical = 1;
else
    Subcortical = 0;
end

% Parcellation applied to the imaging data (see tutorial):
% Parcellation = 'Schaefer100'; 
% Parcellation = 'glasser360'; 


% Repetition time (TR) of the fMRI data (if unknown set to 1):
TR = 3;
% Maximum number of TRs for all fMRI sessions:
Tmax = 240;
% Apply temporal filtering to data (0: no; 1: yes)
apply_filter = 0; 
% Lowpass frequency of filter (default 0.1):
flp = 0.1;
% Highpass frequency of filter (default 0.01):
fhi = 0.01;
% For the statistics: 
% Choose 0 (unpaired) if subjects in different conditions are not the
% same; or 1 (paired) if subjects are the same across conditions.
Paired_tests = 0;
% Number of permutations. For the first analysis to be relatively quick, 
% choose like 500, but then increase to 10000 to increase the precision of
% the final statistical results for publication. 
n_permutations=10000;
% Number of bootstraps. For the first analysis to be relatively quick, 
% choose 1, but then increase for final results to be more precise. 
n_bootstraps=1;
% Direction to plot the FC states/brain ('SideView' or 'TopView'):
CortexDirection = 'TopView';

% AFTER FILLING IN THE INPUT PARAMETERS:
% ||||||||||||||||||||||||||||||| CLICK RUN |||||||||||||||||||||||||||||||

% Go to the directory containing the LEiDA functions
cd(LEiDA_directory)
% Create a directory to store the results from the current LEiDA run
leida_results = [LEiDA_directory 'LEiDA_Results_' run_name '/'];
if ~exist(leida_results, 'dir')
    mkdir(leida_results);
end
addpath(genpath(LEiDA_directory))

%% B: RUN LEADING EIGENVECTOR DYNAMICS ANALYSIS

% Compute the leading eigenvectors of the data
LEiDA_data(Data_directory,Parcellation,Subcortical,leida_results,N_areas,Tmax,apply_filter,flp,fhi,TR);

SelectSample345(leida_results)

% Cluster the leading eigenvectors of all subjects
LEiDA_cluster(leida_results);

[Index_Conditions, Condition_tags] = Get_Diagnosis_from_CSV(leida_results);

% Compute the fractional occupancy and perform hypothesis tests
LEiDA_stats_FracOccup(leida_results,Index_Conditions,Condition_tags,Paired_tests,n_permutations,n_bootstraps);

%% C: MAKE FIGURES
% Plot the centroids obtained using LEiDA and their overlap with Yeo nets
Plot_Centroid_Pyramid(leida_results,Condition_tags,Parcellation,Subcortical,N_areas,CortexDirection)

%CortexDirection = 'SideView';
%Plot_Centroid_Pyramid(leida_results,Condition_tags,Parcellation,Subcortical,N_areas,CortexDirection)

% Generate and save the p-value and barplot plots for fractional occupancy
Plot_FracOccup(leida_results)

end

% % 
% %% 
% selectedK=8;
% Plot_K_overlap_yeo_nets(LEiDA_directory,leida_results,selectedK,Parcellation)
% %