function LEiDA_Start

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function serves to run LEiDA on any given dataset.
%
% Since we do not know a priori the optimal number of FC states that
% differentiate between conditions, this function analyzes the solutions
% obtained for different numbers of FC states (K).
%
% After analysing the output figures saved in the results folder, the user
% can then choose the optimal number of FC states for subsequent analysis
% using the functions LEiDA_AnalysisK.m and LEiDA_AnalysisCentroid.m
%
% This function contains three sections:
%      (A) User defines the parameters and properties of the study.
%      (B) Compute the leading eigenvectors of all scans and cluster them.
%      (C) Assign the scans to conditions, compute statistics and figures.
%
% Start by reading the README.md file.
%
% A: User input parameters
% B: Run Leading Eigenvector Dynamics Analysis:
%    - (Optional) Parcellate data in MNI space
%    - Compute the leading eigenvectors for all scans
%    - Cluster the leading eigenvectors of all scans
% C: Study design, statistics and figures
%    - Assign each scan to a condition (from file names or metadata table)
%    - Fractional occupancy and hypothesis tests between conditions
%    - Pyramid of FC states
%
% Tutorial: README.md
% Version:  V2.0, October 2026
% Authors:  Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

%% A: STUDY PARAMETERS

% ------------------------------ DIRECTORIES ------------------------------
% Folder with the input data. Files are searched in this folder and in all
% its subfolders (e.g. sub-XXX/ses-YYY/):
Data_directory = '/Users/user/Documents/Research/MIGRAINE/func_task-rest/data/parcellated_data/';
% Folder where the results will be saved (a subfolder named after run_name
% is created inside it):
Results_directory = '/Users/user/Documents/Research/MIGRAINE/LEiDA_Results/';
% Name of the run, used to name the folder with the results:
run_name = 'MIGRAINE_AAL90';

% ------------------------------ INPUT DATA -------------------------------
% Space of the input data:
%   'parcellated' : one file per scan with the time series of each area
%   'MNI'         : one 4D NIfTI file per scan in MNI 2mm space, which is
%                   first parcellated using the atlas in Parcellation
Data_space = 'parcellated';
% Files to use (* is a wildcard). Only files matching this pattern are
% considered, e.g. '*BOLDpreprocessed_rp_mo_csf_wm_nui-parcellatedAAL90.mat'
% for parcellated data or '*.nii.gz' for data in MNI space:
File_pattern = '*BOLDpreprocessed_rp_mo_csf_wm_nui-parcellatedAAL90.mat';
% Orientation of the parcellated data ('AreasxTime', 'TimexAreas' or
% 'auto' to detect it from the number of areas in the atlas):
Data_orientation = 'AreasxTime';
% Name of the variable with the time series inside the .mat files
% ('' to detect it automatically):
Data_variable = '';
% Optional text file with the IDs of the subjects to include (one per
% line, e.g. sub-control019). Leave '' to use all files matching
% File_pattern:
Subject_list = '';

% ----------------------------- PARCELLATION ------------------------------
% Atlas used to parcellate the data. Must be one of the atlases available
% in ParcelsMNI2mm.mat to parcellate data in MNI space and to render the FC
% states in the brain: 'AAL116', 'AAL120', 'dbs80', 'Desikan',
% 'Schaefer100', 'glasser360'. Set to '' if the data were parcellated with
% another atlas (the brain renderings are then skipped).
Parcellation = 'AAL116';
% Number of areas to consider (the first N_areas areas of the data). Leave
% empty [] to use all areas in the data (or all parcels of the atlas for
% data in MNI space). E.g. 90 to use only the cerebrum areas of AAL116.
N_areas = 90;

% --------------------------- SIGNAL PROCESSING ---------------------------
% Repetition time (TR) of the fMRI data in seconds (if unknown set to 1):
TR = 2;
% Apply temporal filtering to data (0: no; 1: yes)
apply_filter = 0;
% Lowpass frequency of filter (default 0.1):
flp = 0.1;
% Highpass frequency of filter (default 0.01):
fhi = 0.01;

% ------------------------------ CLUSTERING -------------------------------
% Range of numbers of FC states (K) to compute:
rangeK = 2:20;

% ----------------------------- STUDY DESIGN ------------------------------
% How to assign each scan to a condition:
%   'filename' : a scan belongs to the condition whose tag (in
%                Condition_tags) is found in its file path (relative to
%                Data_directory, e.g. sub-patient001/ses-ictal/...), e.g.
%                Condition_tags = {'sub-control','sub-patient'} or
%                Condition_tags = {'sub-control','ses-interictal','ses-ictal'}
%                A condition can group several tags, given as a cell, e.g.
%                {{'ses-interictal','ses-preictal','ses-postictal'},'ses-ictal'}
%   'metadata' : the condition is read from a table (.csv or .xlsx) with
%                one row per subject (or per subject and session)
%   'none'     : no conditions (only eigenvectors, clustering and centroids)
Condition_source = 'filename';
% Metadata table (only used if Condition_source = 'metadata'). In
% spreadsheets with several sheets, all sheets are read and the name of
% the sheet is available as a column named 'sheet':
Metadata_file = '/Users/user/Documents/Research/MIGRAINE/Demographic_Clinical_Data.xlsx';
% Column with the subject IDs (as found in the folder/file names, e.g.
% sub-control019):
Metadata_subject_column = 'subject';
% Column with the session IDs (e.g. ses-midcycle). Leave '' if there is
% only one scan per subject:
Metadata_session_column = 'session';
% Column with the condition of each scan:
Condition_column = 'sheet';
% Conditions to compare. With 'filename', the tags found in the file
% paths; with 'metadata', the values of Condition_column. Scans not
% belonging to any of these conditions are left out of the statistics.
% Leave {} with 'metadata' to use all values found in Condition_column.
% Examples:
%   controls vs patients ('metadata' with Condition_column = 'sheet'):
%       Condition_tags = {'controls','patients'};
%   controls vs patients interictal vs patients ictal ('filename'):
%       Condition_tags = {'sub-control','ses-interictal','ses-ictal'};
%       Condition_names = {'Controls','Interictal','Ictal'};
%   patients without vs with migraine attack ('filename'):
%       Condition_tags = {{'ses-interictal','ses-preictal','ses-postictal'},'ses-ictal'};
%       Condition_names = {'Non-ictal','Ictal'};
%   phases of the migraine cycle in patients ('filename'):
%       Condition_tags = {'ses-interictal','ses-preictal','ses-ictal','ses-postictal'};
%       Condition_names = {'Interictal','Pre-ictal','Ictal','Post-ictal'};
Condition_tags = {'ses-interictal','ses-preictal','ses-ictal','ses-postictal'};
% Names of the conditions in the figures (same order as Condition_tags).
% Leave {} to use Condition_tags.
Condition_names = {'Interictal','Pre-ictal','Ictal','Post-ictal'};
% Name of this study design. The statistics, conditions and figures of the
% design are also copied to a subfolder Design_<Design_name> of the results
% folder, so that they are kept when another design is run ('' to skip).
Design_name = 'Migraine_cycle';

% ------------------------------ STATISTICS -------------------------------
% Differences between conditions are tested with permutation tests:
%   0 (unpaired): subjects in different conditions are not the same.
%      Statistic: Welch's t (does not assume equal variances).
%   1 (paired): the same subjects are scanned in both conditions; each
%      subject's scans in the two conditions are paired by subject ID.
%      Statistic: paired t.
%   'auto': for each pair of conditions, paired if they share subjects
%      (e.g. interictal vs ictal sessions of the same patients) and
%      unpaired otherwise (e.g. controls vs patients).
Paired_tests = 'auto';
% Observations used in the tests:
%   'subject': the scans of each subject in a condition are averaged, so
%              each subject counts once (recommended when subjects have
%              more than one scan per condition, and required for paired
%              tests if they do)
%   'scan'   : each scan is treated as an independent observation
Stats_level = 'subject';
% Number of permutations. The smallest p-value that can be obtained is
% 1/(n_permutations+1). To reach the strictest threshold, 0.05/sum(rangeK)
% (= 2.4e-4 for K = 2:20), more than 4200 permutations are
% needed. For a first quick analysis choose like 1000, and 10000 or more
% for the final results.
n_permutations = 20000;
% Number of bootstrap samples for the 95% confidence interval of the effect
% sizes (0 to skip; e.g. 1000 for the final results).
n_bootstraps = 1000;

% -------------------------------- FIGURES --------------------------------
% Direction to plot the FC states/brain ('SideView' or 'TopView'):
CortexDirection = 'TopView';

% ------------------------------ STEPS TO RUN -----------------------------
% Set to 0 to skip steps already computed in a previous run (e.g. to test
% a new study design without recomputing the eigenvectors and clusters):
Run_eigenvectors = 0;
Run_clustering = 0;

% AFTER FILLING IN THE INPUT PARAMETERS:
% ||||||||||||||||||||||||||||||| CLICK RUN |||||||||||||||||||||||||||||||

% Add the LEiDA functions (the folder containing this file) to the path
LEiDA_directory = fileparts(mfilename('fullpath'));
addpath(genpath(LEiDA_directory))

% Create a directory to store the results from the current LEiDA run
leida_results = [fullfile(Results_directory, ['LEiDA_Results_' run_name]) filesep];
if ~exist(leida_results, 'dir')
    mkdir(leida_results);
end

% Number of parcels in the atlas ([] if not available in ParcelsMNI2mm)
N_parcels = Atlas_Nparcels(Parcellation);
if isempty(N_areas) && strcmpi(Data_space,'MNI')
    N_areas = N_parcels;
end

% Save the parameters of this run for reference and for later steps
save([leida_results 'LEiDA_Parameters.mat'], 'Data_directory', 'Data_space', ...
     'File_pattern', 'Data_orientation', 'Data_variable', 'Subject_list', ...
     'Parcellation', 'N_areas', 'TR', 'apply_filter', 'flp', 'fhi', 'rangeK', ...
     'Condition_source', 'Metadata_file', 'Metadata_subject_column', ...
     'Metadata_session_column', 'Condition_column', 'Condition_tags', 'Condition_names', 'Design_name', ...
     'Paired_tests', 'Stats_level', 'n_permutations', 'n_bootstraps', 'CortexDirection');

%% B: RUN LEADING EIGENVECTOR DYNAMICS ANALYSIS

if Run_eigenvectors
    % Parcellate the data in MNI space with the selected atlas
    if strcmpi(Data_space,'MNI')
        parcellated_dir = [leida_results 'Parcellated_' Parcellation filesep];
        Parcellate(Data_directory,File_pattern,parcellated_dir,Parcellation,N_areas);
        Data_directory = parcellated_dir;
        File_pattern = ['*_' Parcellation '.mat'];
        Data_orientation = 'AreasxTime';
        Data_variable = 'data';
    end

    % Compute the leading eigenvectors of the data
    LEiDA_data(Data_directory,File_pattern,leida_results,N_areas,N_parcels,...
               Data_orientation,Data_variable,Subject_list,apply_filter,flp,fhi,TR);
end

% Cluster the leading eigenvectors of all scans
if Run_clustering
    LEiDA_cluster(leida_results,rangeK);
end

%% C: STUDY DESIGN, STATISTICS AND FIGURES

% Start time of this step (to keep a copy of the files of this design)
t_design = now;

% Assign each scan to a condition
if strcmpi(Condition_source,'none')
    Condition_tags = {};
else
    [Index_Conditions,Condition_tags] = LEiDA_conditions(leida_results,Condition_source,...
        Condition_tags,Metadata_file,Metadata_subject_column,Metadata_session_column,Condition_column);
    if ~isempty(Condition_names)
        if numel(Condition_names) ~= numel(Condition_tags)
            error('LEiDA:conditions','Condition_names must have one name per condition in Condition_tags');
        end
        Condition_tags = Condition_names(:)';
    end

    % Compute the fractional occupancy and perform hypothesis tests
    LEiDA_stats_FracOccup(leida_results,Index_Conditions,Condition_tags,Paired_tests,n_permutations,n_bootstraps,Stats_level);

    % Print the 10 largest effect sizes of each pair of conditions, with
    % their 95% bootstrap confidence intervals and p-values
    S = load([leida_results 'LEiDA_Stats_FracOccup.mat'],'effectsize_signed','effectsize_ci','P_pval2sided','rangeK');
    [cond2,cond1] = find(triu(ones(numel(Condition_tags)),1)');
    for pair = 1:numel(cond1)
        g_pair = reshape(S.effectsize_signed(pair,:,:),[],1);
        [~,idx] = sort(abs(g_pair),'descend');
        idx = idx(1:min(10,numel(idx)));
        [k,c] = ind2sub([numel(S.rangeK) S.rangeK(end)],idx);
        ci = reshape(S.effectsize_ci(pair,:,:,:),[],2);
        p = reshape(S.P_pval2sided(pair,:,:),[],1);
        disp(['Largest effect sizes, ' Condition_tags{cond1(pair)} ' vs ' Condition_tags{cond2(pair)} ...
              ' (g > 0: ' Condition_tags{cond1(pair)} ' higher):'])
        disp(table(S.rangeK(k)',c,g_pair(idx),ci(idx,1),ci(idx,2),p(idx),...
                   'VariableNames',{'K','c','g','CI_low','CI_high','p'}))
    end

    % Generate and save the p-value and barplot plots for fractional occupancy
    Plot_FracOccup(leida_results)
end

% Plot the centroids obtained using LEiDA and their overlap with Yeo nets
Plot_Centroid_Pyramid(leida_results,Condition_tags,Parcellation,N_areas,CortexDirection)

% Keep a copy of the results of this design in its own subfolder
if ~strcmpi(Condition_source,'none') && ~isempty(Design_name)
    design_dir = [leida_results 'Design_' Design_name filesep];
    if ~exist(design_dir,'dir')
        mkdir(design_dir);
    end
    design_files = [dir([leida_results 'LEiDA_Stats_FracOccup.mat']); dir([leida_results 'LEiDA_Conditions.mat']);
                    dir([leida_results 'FracOccup_*']); dir([leida_results 'Centroid_Pyramid_*'])];
    % Only the files written by this design (not those of previous designs)
    design_files = design_files([design_files.datenum] >= t_design - 1/86400);
    for f = 1:numel(design_files)
        copyfile(fullfile(design_files(f).folder,design_files(f).name),design_dir);
    end
    copyfile([leida_results 'LEiDA_Parameters.mat'],design_dir);
    disp(['Results of this design copied to ' design_dir])
end
