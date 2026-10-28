function LEiDA_Main

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  LEADING EIGENVECTOR DYNAMICS ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Main script of the LEiDA toolbox: all the parameters of the analysis are
% set here, and the analysis is run in steps:
%
%   1) LEiDA_Start: leading eigenvectors of all scans and their clustering
%      into K = rangeK FC states (computed once per run; set
%      Run_eigenvectors / Run_clustering to 1 to compute them again)
%   2) LEiDA_Design: for the selected study design, assign the scans to
%      conditions, compare the fractional occupancy of the FC states
%      between conditions and plot the results for all K (p-values, effect
%      sizes, barplots and the centroid pyramid)
%   3) LEiDA_AnalysisK: figures of all the FC states of selected values of K
%   4) LEiDA_AnalysisCentroid: summary figure of selected FC states
%
% Results are saved in Results_directory/LEiDA_Results_<run_name>/, and the
% results of each study design in its subfolder Design_<design name>/
% (figures of K and of single FC states in K<K>/ and K<K>C<c>/ inside it).
%
% Several study designs can be defined below; Design_name selects the one
% analysed. To analyse another design, change Design_name and run again
% (the eigenvectors and clusters are reused).
%
% Start by reading the README.md file.
%
% Tutorial: README.md
% Version:  V2.0, October 2026
% Authors:  Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

%% 1: DATA AND RUN

% ------------------------------ DIRECTORIES ------------------------------
% Folder with the input data. Files are searched in this folder and in all
% its subfolders (e.g. sub-XXX/ses-YYY/):
P.Data_directory = '/Users/user/Documents/Research/MIGRAINE/func_task-rest/data/parcellated_data/';
% Folder where the results will be saved (a subfolder LEiDA_Results_<run_name>
% is created inside it):
P.Results_directory = '/Users/user/Documents/Research/MIGRAINE/LEiDA_Results/';
% Name of the run, used to name the folder with the results:
P.run_name = 'MIGRAINE_AAL90';

% ------------------------------ INPUT DATA -------------------------------
% Space of the input data:
%   'parcellated' : one file per scan with the time series of each area
%   'MNI'         : one 4D NIfTI file per scan in MNI 2mm space, which is
%                   first parcellated using the atlas in Parcellation
P.Data_space = 'parcellated';
% Files to use (* is a wildcard). Only files matching this pattern are
% considered, e.g. '*BOLDpreprocessed_rp_mo_csf_wm_nui-parcellatedAAL90.mat'
% for parcellated data or '*.nii.gz' for data in MNI space:
P.File_pattern = '*BOLDpreprocessed_rp_mo_csf_wm_nui-parcellatedAAL90.mat';
% Orientation of the parcellated data ('AreasxTime', 'TimexAreas' or
% 'auto' to detect it from the number of areas in the atlas):
P.Data_orientation = 'AreasxTime';
% Name of the variable with the time series inside the .mat files
% ('' to detect it automatically):
P.Data_variable = '';
% Optional text file with the IDs of the subjects to include (one per
% line, e.g. sub-control019). Leave '' to use all files matching
% File_pattern:
P.Subject_list = '';

% ----------------------------- PARCELLATION ------------------------------
% Atlas used to parcellate the data, needed to parcellate data in MNI
% space and to render the FC states in the brain. Either:
%   - the full path of the atlas as a NIfTI file with integer labels in MNI
%     space, at its own resolution (e.g. the same file used to parcellate
%     the data). Names of the areas can be given in a text file with the
%     same name and extension .txt/.tsv/.csv (one area per line, see
%     Load_atlas). E.g. '/Users/user/Atlases/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz'
%   - one of the atlases included in ParcelsMNI2mm.mat (FSL MNI152 2mm
%     grid): 'AAL116', 'AAL120', 'dbs80', 'Desikan', 'Schaefer100',
%     'glasser360'
% Set to '' if no atlas is available (the brain figures are then skipped).
P.Parcellation = 'AAL116';
% Number of areas to consider (the first N_areas areas of the data). Leave
% empty [] to use all areas in the data (or all parcels of the atlas for
% data in MNI space). E.g. 90 to use only the cerebrum areas of AAL116.
P.N_areas = 90;

% --------------------------- SIGNAL PROCESSING ---------------------------
% Repetition time (TR) of the fMRI data in seconds (if unknown set to 1):
P.TR = 2;
% Apply temporal filtering to data (0: no; 1: yes)
P.apply_filter = 0;
% Lowpass frequency of filter (default 0.1):
P.flp = 0.1;
% Highpass frequency of filter (default 0.01):
P.fhi = 0.01;

% ------------------------------ CLUSTERING -------------------------------
% Range of numbers of FC states (K) to compute:
P.rangeK = 2:20;

% -------------------------------- FIGURES --------------------------------
% Direction to plot the FC states/brain in the pyramid ('SideView' or
% 'TopView'):
P.CortexDirection = 'TopView';

%% 2: STUDY DESIGNS

% ---------------------- SETTINGS COMMON TO ALL DESIGNS -------------------
% (each design below can override any of them)
%
% Differences between conditions are tested with permutation tests:
%   0 (unpaired): subjects in different conditions are not the same.
%      Statistic: Welch's t (does not assume equal variances).
%   1 (paired): the same subjects are scanned in both conditions; each
%      subject's scans in the two conditions are paired by subject ID.
%      Statistic: paired t.
%   'auto': for each pair of conditions, paired if they share subjects
%      (e.g. two sessions of the same patients) and unpaired otherwise
%      (e.g. controls vs patients).
Common.Paired_tests = 'auto';
% Observations used in the tests:
%   'subject': the scans of each subject in a condition are averaged, so
%              each subject counts once (recommended when subjects have
%              more than one scan per condition)
%   'scan'   : each scan is treated as an independent observation
Common.Stats_level = 'subject';
% Number of permutations. The smallest p-value that can be obtained is
% 1/(n_permutations+1). To reach the strictest threshold, 0.05/sum(rangeK)
% (= 2.4e-4 for K = 2:20), more than 4200 permutations are needed. For a
% first quick analysis choose like 1000, and 10000 or more for the final
% results.
Common.n_permutations = 20000;
% Number of bootstrap samples for the 95% confidence interval of the effect
% sizes (0 to skip; e.g. 1000 for the final results):
Common.n_bootstraps = 1000;
% How to assign each scan to a condition:
%   'filename' : a scan belongs to the condition whose tag (in
%                Condition_tags) is found in its file path (relative to
%                Data_directory, e.g. sub-patient001/ses-ictal/...). A
%                condition can group several tags, given as a cell, e.g.
%                {{'ses-interictal','ses-preictal','ses-postictal'},'ses-ictal'}
%   'metadata' : the condition is read from a table (.csv or .xlsx) with
%                one row per subject (or per subject and session), in the
%                column Condition_column; Condition_tags are its values
%   'none'     : no conditions (only the centroid pyramid)
Common.Condition_source = 'filename';
% Metadata table (only used if Condition_source = 'metadata'). In
% spreadsheets with several sheets, all sheets are read and the name of
% the sheet is available as a column named 'sheet':
Common.Metadata_file = '/Users/user/Documents/Research/MIGRAINE/Demographic_Clinical_Data.xlsx';
% Columns with the subject IDs and session IDs (as found in the folder/file
% names, e.g. sub-control019 and ses-midcycle; session '' if one scan per
% subject):
Common.Metadata_subject_column = 'subject';
Common.Metadata_session_column = 'session';

% -------------------------------- DESIGNS --------------------------------
% Each design is Designs.<name>, with:
%   Condition_tags   conditions to compare (tags in the file paths, or
%                    values of Condition_column); scans not in any of them
%                    are left out
%   Condition_names  names of the conditions in the figures (same order)
%   Condition_column (only for 'metadata') column with the condition

d = 'Controls_vs_Patients';
Designs.(d).Condition_source = 'metadata';
Designs.(d).Condition_column = 'sheet';
Designs.(d).Condition_tags = {'controls','patients'};
Designs.(d).Condition_names = {'Controls','Patients'};

d = 'Controls_Interictal_Ictal';
Designs.(d).Condition_tags = {'sub-control','ses-interictal','ses-ictal'};
Designs.(d).Condition_names = {'Controls','Interictal','Ictal'};

d = 'NonIctal_vs_Ictal';
Designs.(d).Condition_tags = {{'ses-interictal','ses-preictal','ses-postictal'},'ses-ictal'};
Designs.(d).Condition_names = {'Non-ictal','Ictal'};

d = 'Migraine_cycle';
Designs.(d).Condition_tags = {'ses-interictal','ses-preictal','ses-ictal','ses-postictal'};
Designs.(d).Condition_names = {'Interictal','Pre-ictal','Ictal','Post-ictal'};

% Design analysed in this run:
Design_name = 'Migraine_cycle';

%% 3: STEPS TO RUN

% Compute the eigenvectors and the clusters again (1), or reuse those
% computed in a previous run (0; they are computed if not found):
Run_eigenvectors = 0;
Run_clustering = 0;
% Save the centroids of all K as NIfTI files (1) in Centroids_NIfTI/ (one
% 4D file per K, one volume per FC state), e.g. to view them in FSLeyes:
Save_centroids_nifti = 1;
% Statistics and figures of the study design for all K (1), or reuse
% those computed in a previous run (0):
Run_design = 1;
% Values of K whose FC states are analysed (LEiDA_AnalysisK), e.g. [6 9];
% [] for none:
Analyse_K = 6;
% FC states analysed individually (LEiDA_AnalysisCentroid), one row
% [K c] per FC state, e.g. [6 4; 9 2]; [] for none:
Analyse_states = [6 4];

% AFTER FILLING IN THE INPUT PARAMETERS:
% ||||||||||||||||||||||||||||||| CLICK RUN |||||||||||||||||||||||||||||||

% Add the LEiDA functions (the folder containing this file) to the path
addpath(genpath(fileparts(mfilename('fullpath'))))

% Study design: common settings and those of the selected design
if ~isfield(Designs,Design_name)
    error('LEiDA:main','Design ''%s'' is not defined. Designs: %s',Design_name,strjoin(fieldnames(Designs),', '));
end
P.Design = Common;
for f = fieldnames(Designs.(Design_name))'
    P.Design.(f{1}) = Designs.(Design_name).(f{1});
end
P.Design_name = Design_name;
P.Run_eigenvectors = Run_eigenvectors;
P.Run_clustering = Run_clustering;

% Folders of the results of the run and of the design
P.leida_results = [fullfile(P.Results_directory,['LEiDA_Results_' P.run_name]) filesep];
P.design_dir = [P.leida_results 'Design_' Design_name filesep];

% 1) Eigenvectors and clustering
P = LEiDA_Start(P);
if Save_centroids_nifti
    Export_centroids_nifti(P);
end

% 2) Statistics and figures of the study design for all K
if Run_design
    LEiDA_Design(P);
end

% 3) FC states of selected values of K
for K = Analyse_K(:)'
    LEiDA_AnalysisK(P,K);
end

% 4) Selected FC states
for s = 1:size(Analyse_states,1)
    LEiDA_AnalysisCentroid(P,Analyse_states(s,1),Analyse_states(s,2));
end
