function LEiDA_AnalysisK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to analyse the set of K FC states (PL states) chosen according
% to the analysis from LEiDA_Start. The user should select a value for K
% and this function will return summary plots of the set of FC states.
% This function contains two sections: (A) user input parameters; and (B)
% code to plot the figures for the selected K. The user should only need
% to adapt section A.
%
% The parcellation and the statistics are read from the results of
% LEiDA_Start (LEiDA_Parameters.mat, LEiDA_Clusters.mat and the
% LEiDA_Stats_FracOccup.mat of the study design).
%
% B: Analysis plots for selected K:
%    - FC states as links in cortex
%    - FC states as nodes in cortex
%    - FC states in 3D glass brain, colored by their Yeo network
%    - FC states in matrix format
%    - FC states in vector format with number of each parcel
%    - Fractional occupancy of each FC state in each condition
%    - Summary of the repertoire of K FC states
%    - Overlap of the FC states with the Yeo networks (Yeo et al. 2011)
%
% Tutorial: README.md
% Version:  V2.0, October 2026
% Authors:  Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

%% A: USER INPUT PARAMETERS

% Folder with the results and name of the run (as in LEiDA_Start):
Results_directory = '/Users/user/Documents/Research/MIGRAINE/LEiDA_Results/';
run_name = 'MIGRAINE_AAL90';
% Study design whose statistics are shown (Design_name in LEiDA_Start).
% Leave '' to use the statistics of the last design run (in the results
% folder of the run):
Design_name = 'Migraine_cycle';
% Number of FC states K to analyse:
SelectK = 6;

% AFTER FILLING IN THE INPUT PARAMETERS:
% ||||||||||||||||||||||||||||||| CLICK RUN |||||||||||||||||||||||||||||||

% Add the LEiDA functions (the folder containing this file) to the path
LEiDA_directory = fileparts(mfilename('fullpath'));
addpath(genpath(LEiDA_directory))

% Results of the run and parameters used in LEiDA_Start
leida_results = [fullfile(Results_directory,['LEiDA_Results_' run_name]) filesep];
if ~isfile([leida_results 'LEiDA_Parameters.mat'])
    error('LEiDA:analysisK','LEiDA_Parameters.mat not found in %s; run LEiDA_Start first',leida_results);
end
load([leida_results 'LEiDA_Parameters.mat'],'Parcellation');
if isempty(Atlas_Nparcels(Parcellation))
    error('LEiDA:analysisK','The figures need a parcellation available in ParcelsMNI2mm.mat');
end

% Statistics of the study design, and folder to save the figures
if isempty(Design_name)
    design_dir = leida_results;
else
    design_dir = [leida_results 'Design_' Design_name filesep];
end
stats_file = [design_dir 'LEiDA_Stats_FracOccup.mat'];
if ~isfile(stats_file)
    warning('LEiDA:analysisK','%s not found: the fractional occupancy is not shown',stats_file);
    stats_file = '';
end
K_dir = [design_dir 'K' num2str(SelectK) filesep];
if ~exist(K_dir,'dir')
    mkdir(K_dir);
end

%% B: FIGURES FOR THE K FC STATES

disp(' ')
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES FOR K = ' num2str(SelectK) ' FC STATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

% FC states as links in cortex
Plot_K_links_in_cortex(leida_results,K_dir,SelectK,Parcellation);

% FC states as nodes in cortex
Plot_K_nodes_in_cortex(leida_results,K_dir,SelectK,Parcellation);

% FC states in 3D glass brain
Plot_K_3Dbrain(leida_results,K_dir,SelectK,Parcellation);

% FC states in matrix format
Plot_K_matrix(leida_results,K_dir,SelectK,Parcellation);

% FC states in vector format with number of each parcel
Plot_K_vector_numbered(leida_results,K_dir,SelectK,Parcellation);

% Fractional occupancy of each FC state in each condition
if ~isempty(stats_file)
    Plot_K_boxplot_FO(leida_results,stats_file,K_dir,SelectK,Parcellation);
end

% Summary of the repertoire of K FC states
Plot_K_repertoire(leida_results,stats_file,K_dir,SelectK,Parcellation);

% Overlap between FC states and Yeo networks
Plot_K_overlap_yeo_nets(leida_results,K_dir,SelectK,Parcellation)
