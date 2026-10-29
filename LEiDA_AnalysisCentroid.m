function LEiDA_AnalysisCentroid(P,SelectK,Centroid)
%
% Fourth step of LEiDA (called from LEiDA_Main): summary figure of one FC
% state (centroid), e.g. a FC state whose fractional occupancy differs
% between conditions, saved in the folder K<SelectK>C<Centroid> of the
% design (P.design_dir). See Plot_Centroid_Summary:
%    - the FC state rendered in the brain (side and top views), colored
%      with its Yeo network
%    - the strongest links between the areas of the FC state
%    - its overlap with the 7 Yeo networks (weighted Dice and fraction)
%    - the areas of the FC state, ordered by their value in the centroid
%    - its fractional occupancy in each condition, with the p-value, effect
%      size and confidence interval of each comparison
%
% INPUT:
% P         parameters set in LEiDA_Main (after LEiDA_Start)
% SelectK   number of FC states K
% Centroid  FC state to analyse (1 <= Centroid <= SelectK)
%
% Tutorial: README.md
% Version:  V2.0, October 2026
% Authors:  Joana Cabral, Universidade do Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

if nargin < 3
    error('LEiDA:centroid','Set the parameters and run the analysis from LEiDA_Main (Analyse_states)');
end

% Statistics of the study design, and folder to save the figure
stats_file = [P.design_dir 'LEiDA_Stats_FracOccup.mat'];
if ~isfile(stats_file)
    warning('LEiDA:centroid','%s not found (run the design with Run_design = 1): the fractional occupancy is not shown',stats_file);
    stats_file = '';
end
save_dir = [P.design_dir 'K' num2str(SelectK) 'C' num2str(Centroid) filesep];

disp(' ')
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%% FC STATE ' num2str(Centroid) ' (K = ' num2str(SelectK) ') %%%%%%%%%%%%%%%%%%%%%%%%%%%%'])

Plot_Centroid_Summary(P.leida_results,stats_file,save_dir,SelectK,Centroid,P.Parcellation,P.N_areas);
