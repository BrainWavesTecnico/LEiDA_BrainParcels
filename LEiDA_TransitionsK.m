function LEiDA_TransitionsK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to compute the transition probability matrix for each subject
% and compare the transition probabilities between conditions for a value
% of K chosen according to the analysis from LEiDA_Start. This function
% contains two sections: (A) user input parameters; and (B) code to compute
% the transition probability matrix and intergroup differences in
% transition probabilities. The user should only need to adapt section A.
%
% Start by reading the README.md file.
%
% A: User input parameters
% B: Analysis for selected K:
%    - Compute Transition Probability Matrix (TPM) for each participant
%    - Compare state-to-state transition probabilities between conditions
%    - Plot the mean TPM for each condition
%    - Plot the intergroup differences between transition probabilties
%
% Tutorial: README.md
% Version:  V1.0, May 2022
% Authors:  Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

%% A: USER INPUT PARAMETERS

% Define K value, i.e., K returning the most significant differences between conditions:
SelectK = 15;

% Directory of the LEiDA toolbox folder:
LEiDA_directory = 'D:/LEiDA_Toolbox/';
% Name of the run to be used to create the folder to save the data:
run_name = 'ABIDE_dparsf_AAL120';

% AFTER FILLING IN THE INPUT PARAMETERS:
% ||||||||||||||||||||||||||||||| CLICK RUN |||||||||||||||||||||||||||||||

%% COMPUTE TPM AND ANALYSE INTERGROUP DIFFERENCES IN TRANSITION PROBABILITIES FOR SELECTED K


% Go to the directory containing the LEiDA functions
cd(LEiDA_directory)

% Directory with the results from LEiDA
leida_res = [LEiDA_directory 'LEiDA_Results_' run_name '/'];

% Create a directory to store results for defined value of K
if ~exist([leida_res 'K' num2str(SelectK) '/'], 'dir')
    mkdir([leida_res 'K' num2str(SelectK) '/']);
end
K_dir = [leida_res 'K' num2str(SelectK) '/'];

disp(' ')
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSITIONS FOR K = ' num2str(SelectK) ' CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

% Compute the transition probability matrix and perform hypothesis tests
% LEiDA_stats_TransitionMatrix(leida_res,K_dir,SelectK);

% Plot the mean transition probability matrix
% Plot_K_tpm(K_dir,SelectK);

% Plot summary of differences in transition probabilities between conditions
Plot_K_diffs_transitions(K_dir,SelectK);