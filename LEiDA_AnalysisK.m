function LEiDA_AnalysisK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to analyse the PL states chosen according to the analysis from
% LEiDA_Start. The user should select a value for K and this function will
% return summary plots of the set of PL states. This function contains two
% sections: (A) user input parameters; and (B) code to plot the figures for
% the select K. The user should only need to adapt section A.
%
% Start by reading the README.md file.
%
% A: User input parameters
% B: Analysis plots for selected K:
%    - Plot PL states as links in cortex
%    - Plot PL states as nodes in cortex
%    - Plot PL states in 3D glass brain
%    - Plot PL states in matrix format
%    - Plot PL states in vector format with number of each parcel
%    - Plot boxplot of fractional occupancy values
%    - Plot of summary information for the repertoire of K PL states
%    - Plot overlap of PL states with RSNs defined by Yeo et al., (2011)
%
% Tutorial: README.md
% Version:  V1.0, May 2022
% Authors:  Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

%% A: USER INPUT PARAMETERS

% Define K value, i.e., K returning the most significant differences between conditions:
SelectK = 9;

% Name of the run to be used to create the folder to save the data:
run = 3;

% Directory of the LEiDA toolbox folder:
LEiDA_directory = '/Users/user/Documents/Research/CognitiveDecline/DONAU/LEiDA_Toolbox_May2025/';

run_names = {'DONAU_1_Schaefer_100_cortical';
             'DONAU_2_Schaefer_119_cortical_subcortical';
             'DONAU_3_Desikan_70_cortical';
             'DONAU_4__Desikan_89_cortical_subcortical';
             'DONAU_5_Glasser_360_cortical';
             'DONAU_6_Glasser_379_cortical_subcortical'};

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

addpath(genpath(LEiDA_directory))

% AFTER FILLING IN THE INPUT PARAMETERS:
% ||||||||||||||||||||||||||||||| CLICK RUN |||||||||||||||||||||||||||||||

% Go to the directory containing the LEiDA functions
cd(LEiDA_directory)

% Directory with the results from LEiDA
leida_res = [LEiDA_directory 'LEiDA_Results_' run_name '/'];
% Create a directory to store results for defined value of K
if ~exist([leida_res 'K' num2str(SelectK) '/'], 'dir')
    mkdir([leida_res 'K' num2str(SelectK) '/']);
end
K_dir = [leida_res 'K' num2str(SelectK) '/'];


%% ANALYSE K CENTROIDS SELECTED ACCORDING TO OUTPUT FROM LEiDA_Start


disp(' ')
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES FOR K = ' num2str(SelectK) ' PL STATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

% Plot PL states as links in cortex
Plot_K_links_in_cortex(leida_res,K_dir,SelectK,Parcellation,Subcortical);

% Plot PL states as nodes in cortex
Plot_K_nodes_in_cortex(leida_res,K_dir,SelectK,Parcellation);

% Plot PL states in 3D glass brain
Plot_K_3Dbrain(leida_res,K_dir,SelectK,Parcellation);

% Plot PL states in matrix format
Plot_K_matrix(leida_res,K_dir,SelectK);

% Plot PL states in vector format with number of each parcel
Plot_K_vector_numbered(leida_res,K_dir,SelectK);

% Plot boxplot of fractional occupancy values
Plot_K_boxplot_FO(leida_res,K_dir,SelectK,Parcellation);

% Plot of summary information for the repertoire of K PL states
Plot_K_repertoire(leida_res,K_dir,SelectK,Parcellation);

% Plot overlap between PL states and Yeo RSNs
Plot_K_overlap_yeo_nets(leida_res,K_dir,SelectK,Parcellation)
