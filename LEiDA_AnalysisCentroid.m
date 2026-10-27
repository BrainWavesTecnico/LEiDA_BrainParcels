function LEiDA_AnalysisCentroid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to analyse one PL state from the set of K PL states chosen
% according to the analysis from LEiDA_Start. This script will provide
% detailed figures with relevant information regarding the selected
% PL state. This function contains two sections: (A) user input parameters;
% and (B) code to plot the figures for the selected PL state. The user
% should only need to adapt section A.
%
% Start by reading the README.md file.
%
% A: User input parameters
% B: Analysis plots for selected centroid:
%    - Plot PL state as vector with labels of parcellation
%    - Plot PL state as vector with areas ordered by contribution
%    - Plot boxplot of fractional occupancy values by condition
%    - Plot of summary information for the selected PL state
%
% Tutorial: README.md
% Version:  V1.2, May 2025
% Authors:  Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

%% A: USER INPUT PARAMETERS

% %Define K value, i.e., K returning the most significant differences between conditions:
% SelectK = 16;
% % Define the PL state to be studied (1 <= PL state <= SelectK):
% Centroid = 16;
% % Name of the run to be used to create the folder to save the data:
% run = 4;

% % Name of the run to be used to create the folder to save the data:
% run = 6;
% % Define K value, i.e., K returning the most significant differences between conditions:
% SelectK = 18;
% % Define the PL state to be studied (1 <= PL state <= SelectK):
% Centroid = 10;


% % Define K value, i.e., K returning the most significant differences between conditions:
SelectK = 4;
% Define the PL state to be studied (1 <= PL state <= SelectK):
Centroid = 2;
% Name of the run to be used to create the folder to save the data:
run = 2;

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



% AFTER FILLING IN THE INPUT PARAMETERS:
% ||||||||||||||||||||||||||||||| CLICK RUN |||||||||||||||||||||||||||||||


% Go to the directory containing the LEiDA functions
cd(LEiDA_directory)
addpath(genpath(LEiDA_directory))

% Directory with the results from LEiDA
leida_res = [LEiDA_directory 'LEiDA_Results_' run_name '/'];

% Check if the value for centroid is within the allowed range
if Centroid < 1 || Centroid > SelectK
    error(['Select a cluster c between 1 and ' num2str(SelectK)])
end

% Create a directory to store results for defined value of K
if ~exist([leida_res 'K' num2str(SelectK) 'C' num2str(Centroid) '/'], 'dir')
    mkdir([leida_res 'K' num2str(SelectK) 'C' num2str(Centroid) '/']);
end
Vc_dir = [leida_res 'K' num2str(SelectK) 'C' num2str(Centroid) '/'];


%% ANALYSE PL STATE OF INTEREST

% Close all open figures
close all;


disp(' ')
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES FOR PL STATE ' num2str(Centroid) ' (K = ' num2str(SelectK) ') %%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

% % Plot PL state in vector format with label of each parcel
% Plot_C_vector_labelled(leida_res,Vc_dir,SelectK,Centroid,Parcellation);
% 
% % Plot PL state in vector format with areas organised by contribution
% Plot_C_vector_ordered(leida_res,Vc_dir,SelectK,Centroid,Parcellation);
% 
% % Plot boxplot of the fractional occupancy values for the selected PL state


Plot_centroid_links_FO(leida_res,SelectK,Centroid, Parcellation,Subcortical)



% Plot of summary information for the selected PL state
%Plot_C_summary(leida_res,Vc_dir,SelectK,Centroid,Parcellation,N_areas,Subcortical);

%Plot_K_links_in_cortex(leida_res,Vc_dir,SelectK,Parcellation,Subcortical)