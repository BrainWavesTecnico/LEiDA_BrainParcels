function LEiDA_AnalysisCentroid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to analyse one FC state (PL state) from the set of K FC states
% obtained with LEiDA_Start, e.g. a FC state whose fractional occupancy
% differs between conditions in the centroid pyramid. This function
% contains two sections: (A) user input parameters; and (B) the figure for
% the selected FC state. The user should only need to adapt section A.
%
% The parcellation, number of areas and statistics are read from the
% results of LEiDA_Start (LEiDA_Parameters.mat, LEiDA_Clusters.mat and
% LEiDA_Stats_FracOccup.mat).
%
% The figure (Plot_Centroid_Summary) shows:
%    - the FC state rendered in the brain (side and top views), colored
%      with its Yeo network
%    - the strongest links between the areas of the FC state
%    - its overlap with the 7 Yeo networks (weighted Dice and fraction)
%    - the areas of the FC state, ordered by their value in the centroid
%    - its fractional occupancy in each condition, with the p-value, effect
%      size and confidence interval of each comparison
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
% Number of FC states K and FC state to analyse (1 <= Centroid <= SelectK):
SelectK = 6;
Centroid = 4;

% AFTER FILLING IN THE INPUT PARAMETERS:
% ||||||||||||||||||||||||||||||| CLICK RUN |||||||||||||||||||||||||||||||

% Add the LEiDA functions (the folder containing this file) to the path
LEiDA_directory = fileparts(mfilename('fullpath'));
addpath(genpath(LEiDA_directory))

% Results of the run and parameters used in LEiDA_Start
leida_results = [fullfile(Results_directory,['LEiDA_Results_' run_name]) filesep];
if ~isfile([leida_results 'LEiDA_Parameters.mat'])
    error('LEiDA:centroid','LEiDA_Parameters.mat not found in %s; run LEiDA_Start first',leida_results);
end
load([leida_results 'LEiDA_Parameters.mat'],'Parcellation','N_areas');

% Statistics of the study design, and folder to save the figure
if isempty(Design_name)
    design_dir = leida_results;
else
    design_dir = [leida_results 'Design_' Design_name filesep];
end
stats_file = [design_dir 'LEiDA_Stats_FracOccup.mat'];
if ~isfile(stats_file)
    warning('LEiDA:centroid','%s not found: the fractional occupancy is not shown',stats_file);
    stats_file = '';
end
save_dir = [design_dir 'K' num2str(SelectK) 'C' num2str(Centroid) filesep];

%% B: FIGURE FOR THE SELECTED FC STATE

disp(' ')
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%% FC STATE ' num2str(Centroid) ' (K = ' num2str(SelectK) ') %%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

Plot_Centroid_Summary(leida_results,stats_file,save_dir,SelectK,Centroid,Parcellation,N_areas);
