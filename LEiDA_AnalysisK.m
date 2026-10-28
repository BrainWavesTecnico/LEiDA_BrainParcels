function LEiDA_AnalysisK(P,SelectK)
%
% Third step of LEiDA (called from LEiDA_Main): figures of the set of K FC
% states for a value of K chosen from the results of the study design
% (e.g. the centroid pyramid), saved in the folder K<SelectK> of the
% design (P.design_dir).
%
% Figures:
%    - FC states as links in cortex
%    - FC states as nodes in cortex
%    - FC states in 3D glass brain, colored by their Yeo network
%    - FC states in matrix format
%    - FC states in vector format with number of each parcel
%    - Fractional occupancy of each FC state in each condition
%    - Summary of the repertoire of K FC states
%    - Overlap of the FC states with the Yeo networks (Yeo et al. 2011)
%
% INPUT:
% P         parameters set in LEiDA_Main (after LEiDA_Start)
% SelectK   number of FC states K
%
% Tutorial: README.md
% Version:  V2.0, October 2026
% Authors:  Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

if nargin < 2
    error('LEiDA:analysisK','Set the parameters and run the analysis from LEiDA_Main (Analyse_K)');
end
if ~any(P.rangeK == SelectK)
    error('LEiDA:analysisK','K = %d is not in the range of K computed (%d to %d)',SelectK,P.rangeK(1),P.rangeK(end));
end
if isempty(Atlas_Nparcels(P.Parcellation))
    error('LEiDA:analysisK','The figures need a parcellation available in ParcelsMNI2mm.mat');
end

% Statistics of the study design, and folder to save the figures
stats_file = [P.design_dir 'LEiDA_Stats_FracOccup.mat'];
if ~isfile(stats_file)
    warning('LEiDA:analysisK','%s not found (run the design with Run_design = 1): the fractional occupancy is not shown',stats_file);
    stats_file = '';
end
K_dir = [P.design_dir 'K' num2str(SelectK) filesep];
if ~exist(K_dir,'dir')
    mkdir(K_dir);
end

disp(' ')
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES FOR K = ' num2str(SelectK) ' FC STATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

% FC states as links in cortex
Plot_K_links_in_cortex(P.leida_results,K_dir,SelectK,P.Parcellation);

% FC states as nodes in cortex
Plot_K_nodes_in_cortex(P.leida_results,K_dir,SelectK,P.Parcellation);

% FC states in 3D glass brain
Plot_K_3Dbrain(P.leida_results,K_dir,SelectK,P.Parcellation);

% FC states in matrix format
Plot_K_matrix(P.leida_results,K_dir,SelectK,P.Parcellation);

% FC states in vector format with number of each parcel
Plot_K_vector_numbered(P.leida_results,K_dir,SelectK,P.Parcellation);

% Fractional occupancy of each FC state in each condition
if ~isempty(stats_file)
    Plot_K_boxplot_FO(P.leida_results,stats_file,K_dir,SelectK,P.Parcellation);
end

% Summary of the repertoire of K FC states
Plot_K_repertoire(P.leida_results,stats_file,K_dir,SelectK,P.Parcellation);

% Overlap between FC states and Yeo networks
Plot_K_overlap_yeo_nets(P.leida_results,K_dir,SelectK,P.Parcellation)
