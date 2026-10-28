# Leading Eigenvector Dynamics Analysis Toolbox - Matlab

1. [Installation](#installation)
2. [Workflow: LEiDA_Main](#workflow-leida_main)
3. [Parameters of LEiDA_Main](#parameters-of-leida_main)
4. [Outputs](#outputs)
5. [Statistics](#statistics)
6. [Format of the inputs](#format-of-the-inputs)
7. [Functions](#functions)
8. [Not yet integrated in LEiDA_Main](#not-yet-integrated-in-leida_main)

---

## Installation

Clone or download this repository into a folder. *LEiDA_Main* adds the folder with the toolbox functions to the Matlab path automatically.

> Requires MATLAB R2020b or later, with the Statistics and Machine Learning Toolbox. The Image Processing Toolbox is needed only to parcellate data in MNI space that is not in the 2mm grid of the atlases, and the Parallel Computing Toolbox (optional) speeds up the clustering.

---

## Workflow: LEiDA_Main

All the parameters of the analysis are set in ***LEiDA_Main.m***, which runs the analysis in four steps. Fill in its three parameter sections and click *Run*.

| Step | Function | What it does | Runs when |
|---|---|---|---|
| 1 | *LEiDA_Start* | Leading eigenvectors of the phase coherence of every scan (*LEiDA_data*), and their clustering into K FC states for each K in *rangeK* (*LEiDA_cluster*) | Once per run: reused afterwards unless *Run_eigenvectors* / *Run_clustering* = 1 |
| 2 | *LEiDA_Design* | For the selected study design: conditions of each scan, fractional occupancy of each FC state, permutation tests between conditions, and figures for all K (p-values, effect sizes, barplots, centroid pyramid) | *Run_design* = 1 |
| 3 | *LEiDA_AnalysisK* | Figures of all the FC states of each K in *Analyse_K* | *Analyse_K* not empty |
| 4 | *LEiDA_AnalysisCentroid* | Summary figure of each FC state in *Analyse_states* | *Analyse_states* not empty |

A typical sequence:

1. Run once with the data parameters set to compute the eigenvectors and clusters (they are computed automatically the first time).
2. Look at the centroid pyramid and the p-value plots of the design to choose values of K and FC states of interest.
3. Set *Analyse_K* and *Analyse_states* (and *Run_design* = 0 to reuse the statistics) and run again.
4. To test another study design, change *Design_name* and run again: the eigenvectors and clusters are reused, and each design keeps its own results folder.

---

## Parameters of LEiDA_Main

### 1. Data and run

- ***Data_directory***: folder with the input data. Files are searched in this folder and in all its subfolders (e.g. *sub-XXX/ses-YYY/*);
- ***Results_directory***: folder where the results are saved, in a subfolder *LEiDA_Results_<run_name>*;
- ***run_name***: name of the run;
- ***Data_space***: `'parcellated'` if each file contains the time series of each brain area, or `'MNI'` if the files are 4D NIfTI images in MNI 2mm space, which are first parcellated with the atlas in *Parcellation* (*Parcellate.m*);
- ***File_pattern***: pattern of the file names to use (`*` is a wildcard), e.g. `'*BOLDpreprocessed_rp_mo_csf_wm_nui-parcellatedAAL90.mat'` or `'*.nii.gz'`;
- ***Data_orientation***: `'AreasxTime'`, `'TimexAreas'` or `'auto'`;
- ***Data_variable***: name of the variable (or struct field) with the time series in .mat files; `''` to detect it automatically;
- ***Subject_list***: optional text file with the IDs of the subjects to include (one per line, e.g. `sub-control019`); `''` to use all files;
- ***Parcellation***: atlas of the data. To parcellate data in MNI space and to render the FC states in the brain it must be one of the atlases in *ParcelsMNI2mm.mat*: `'AAL116'`, `'AAL120'`, `'dbs80'`, `'Desikan'`, `'Schaefer100'`, `'glasser360'`; `''` for other atlases (the brain figures are then skipped);
- ***N_areas***: number of areas to consider (the first *N_areas* of the data, e.g. 90 for the cerebrum areas of AAL116); `[]` for all;
- ***TR***, ***apply_filter***, ***flp***, ***fhi***: repetition time (s) and optional band-pass filter (0/1, lowpass and highpass frequencies);
- ***rangeK***: numbers of FC states K to compute (default `2:20`);
- ***CortexDirection***: view of the brains in the pyramid (`'TopView'` or `'SideView'`).

### 2. Study designs

Settings common to all designs (*Common.<field>*), which each design can override:

- ***Paired_tests***: 0 (unpaired, Welch's t), 1 (paired t, subjects paired by ID) or `'auto'` (for each pair of conditions, paired if they share subjects, e.g. two sessions of the same patients, and unpaired otherwise, e.g. controls vs patients);
- ***Stats_level***: `'subject'` to average the scans of each subject within a condition, so that each subject counts once (recommended when subjects have several scans per condition), or `'scan'`;
- ***n_permutations***: number of permutations. The smallest possible p-value is 1/(n_permutations+1), so more than ~4200 permutations are needed to reach the strictest threshold 0.05/sum(rangeK) for K = 2:20 (e.g. 1000 for a first quick analysis, 10000 or more for final results);
- ***n_bootstraps***: number of bootstrap samples for the 95% confidence interval of the effect sizes (0 to skip, e.g. 1000 for final results);
- ***Condition_source***: how to assign each scan to a condition: `'filename'` (a scan belongs to the condition whose tag is found in its file path), `'metadata'` (from a .csv/.xlsx table) or `'none'` (only the centroid pyramid);
- ***Metadata_file***, ***Metadata_subject_column***, ***Metadata_session_column***: metadata table and its columns with the subject and session IDs as in the folder/file names. The header is the first row containing the subject column (rows above it, e.g. units, are ignored); all sheets of a spreadsheet are read and the sheet name is available as a column `'sheet'`.

Each design is defined as *Designs.<name>* with:

- ***Condition_tags***: conditions to compare: tags in the file paths (`'filename'`) or values of *Condition_column* (`'metadata'`). A condition can group several tags, e.g. `{{'ses-interictal','ses-preictal','ses-postictal'},'ses-ictal'}`. Scans not in any condition are left out of the statistics;
- ***Condition_names***: names of the conditions in the figures, in the same order;
- ***Condition_column*** (`'metadata'` only): column with the condition.

***Design_name*** selects the design analysed in the run.

### 3. Steps to run

- ***Run_eigenvectors***, ***Run_clustering***: 1 to compute the eigenvectors / clusters again; 0 to reuse those of a previous run (computed anyway if not found);
- ***Run_design***: 1 to compute the statistics and figures of the design for all K; 0 to reuse them;
- ***Analyse_K***: values of K whose FC states are analysed (e.g. `[6 9]`; `[]` for none);
- ***Analyse_states***: FC states analysed individually, one row `[K c]` per state (e.g. `[6 4; 9 2]`; `[]` for none).

---

## Outputs

```
Results_directory/LEiDA_Results_<run_name>/
    LEiDA_Parameters.mat        parameters of the run
    LEiDA_EigenVectors.mat      leading eigenvectors (V1_all), scan of each (Time_sessions), Scan_info (file, subject, session of each scan)
    LEiDA_Clusters.mat          k-means results for each K (Kmeans_results, rangeK)
    Parcellated_<atlas>/        parcellated data (only if Data_space = 'MNI')
    Design_<design name>/
        LEiDA_Parameters.mat        parameters of the run and the design
        LEiDA_Conditions.mat        condition of each scan, Scan_info with the metadata
        LEiDA_Stats_FracOccup.mat   fractional occupancy, p-values, effect sizes and CIs
        FracOccup_pvalues           p-values across K for each pair of conditions
        FracOccup_effetcsize        effect sizes across K for each pair of conditions
        FracOccup_Barplot_*         pyramids of barplots (all conditions, and each pair)
        Centroid_Pyramid_*          FC states for all K, colored by Yeo network, with the fractional occupancy of each condition
        K<K>/                       figures of LEiDA_AnalysisK
        K<K>C<c>/                   summary figure of LEiDA_AnalysisCentroid
```

All figures are saved as .png (600 dpi) and .fig, and stay open after saving. The console also prints the number of scans in each condition, the test used for each pair of conditions and the 10 largest effect sizes of each pair, with their confidence intervals and p-values.

---

## Statistics

Differences between each pair of conditions are tested for each FC state and each K with permutation tests (*permutation_test.m*, called through *LEiDA_compare_conditions.m*):

- **Unpaired** samples: Welch's t statistic; the null distribution is obtained by randomly permuting the condition labels;
- **Paired** samples: paired t statistic; the null distribution is obtained by randomly swapping the condition labels within each subject (sign flips of the differences);
- The p-values are counted directly from the permutation distribution, p = (1 + number of permuted statistics at least as extreme as the observed one) / (1 + n_permutations), without fitting the null distribution;
- Effect size: Hedges' g (unpaired) or g_z (paired), with a 95% bootstrap (percentile) confidence interval if *n_bootstraps* > 0;
- Significance levels shown in the figures (colors): p < 0.05 (dark red), p < 0.05/K (dark green), p < 0.05/sum(rangeK) (dark blue).
- Overlap with the Yeo et al. (2011) networks (*Overlap_LEiDA_Yeo.m*): each FC state is mapped to MNI space, each voxel weighted by the (positive) centroid value of its parcel divided by the largest, and compared with each of the 7 cortical networks with a weighted Dice coefficient, so that strongly phase-shifted parcels weigh more. A state is colored with the network of largest Dice if it reaches at least 40% of the largest Dice that any set of parcels of the parcellation can reach for that network (coarse parcels limit the overlap with some networks) and at least 0.15, or if at least half of the state (weighted) lies within one network (focal states have a low Dice with large networks); otherwise it is shown in black (e.g. subcortical states).

---

## Format of the inputs

The input data should be contained in the folder specified in *Data_directory* (or its subfolders), with one file per scan. Subject and session IDs are read from BIDS-like paths (`sub-XXX`, `ses-YYY`); if the path has no `sub-XXX`, the file name is used as the subject ID. The time series data can be supplied in the following formats:

- A (no. of brain areas $\times$ no. of time points) or (no. of time points $\times$ no. of brain areas) matrix (see *Data_orientation*) for each scan. The file can be either a Matlab file (.mat), a text file (.txt) or a .1D file. If the .mat file contains several variables or a struct, set *Data_variable* to the name of the variable/field with the time series (by default, a variable or field named `data` is used);
- A 4D NIfTI file (.nii or .nii.gz) in MNI 2mm space for each scan, with *Data_space* set to `'MNI'`.

---

## Functions

Main steps (called from *LEiDA_Main*):

| Function | Description |
|---|---|
| *LEiDA_Main* | Parameters of the analysis and calls to the steps below |
| *LEiDA_Start* | Step 1: eigenvectors and clustering (computed once, reused afterwards) |
| *LEiDA_Design* | Step 2: conditions, statistics and figures of a study design for all K |
| *LEiDA_AnalysisK* | Step 3: figures of the FC states of one K |
| *LEiDA_AnalysisCentroid* | Step 4: summary figure of one FC state |

Analysis:

| Function | Description |
|---|---|
| *Parcellate* | Mean time series of each parcel of an atlas from 4D NIfTI data in MNI 2mm space |
| *LEiDA_data* | Reads the scans (any orientation and file format), repairs missing or flat areas, and computes the leading eigenvector of the phase coherence at each TR |
| *TemporalFiltering* | Band-pass Butterworth filter of the time series |
| *LEiDA_cluster* | k-means clustering (cosine distance) of all eigenvectors for each K; states ordered by occupancy |
| *LEiDA_conditions* | Condition of each scan from tags in the file paths or a metadata table |
| *LEiDA_stats_FracOccup* | Fractional occupancy of each FC state in each scan, and tests between conditions |
| *LEiDA_compare_conditions* | Tests between each pair of conditions (paired/unpaired, per subject or scan) |
| *permutation_test* | Permutation test (Welch's t or paired t), exact permutation p-values, Hedges' g and bootstrap CI |
| *Overlap_LEiDA_Yeo* | Weighted Dice and fraction of each FC state within the 7 Yeo networks, and network assigned |

Figures:

| Function | Description |
|---|---|
| *Plot_FracOccup* | p-values and effect sizes across K, pyramids of barplots of the fractional occupancy |
| *Plot_Centroid_Pyramid* | FC states for all K colored by Yeo network, with the fractional occupancy of each condition, the 7 networks for reference and a legend |
| *Plot_Centroid_Summary* | Summary of one FC state: brain renderings, links, Yeo overlap, areas, fractional occupancy and statistics |
| *Plot_FO_boxplot* | Fractional occupancy of one FC state per subject in each condition, with the significant pairs |
| *Plot_K_links_in_cortex*, *Plot_K_nodes_in_cortex*, *Plot_K_3Dbrain*, *Plot_K_matrix*, *Plot_K_vector_numbered*, *Plot_K_boxplot_FO*, *Plot_K_repertoire*, *Plot_K_overlap_yeo_nets* | Figures of the FC states of one K (see LEiDA_AnalysisK) |

Helpers: *Atlas_Nparcels* (parcels of an atlas), *Parcellation_labels* (area names), *Parcel_centers* (area centers), *Areas_order* (left/right order of AAL areas), *Yeo_colors* (colors and names of the Yeo networks), *Significance_color* (colors of the significance levels), *FracOccup_mean_sem* and *FracOccup_values* (fractional occupancy per subject or scan), *Draw_links_in_cortex*, *Plot_Yeo_reference*, *Pyramid_layout*, *Pyramid_axes*, *Pyramid_labels*, *Pyramid_legend* (layout of the pyramid figures), *Save_png* (600 dpi), *subplot_tight*.

---

## Not yet integrated in LEiDA_Main

These functions have not been revised yet and still have their own parameters (and paths of a previous study) in their first section. They will be integrated as additional steps of *LEiDA_Main*.

### LEiDA Transitions K

The function ***LEiDA_TransitionsK.m*** should be used to analyse the state-to-state transitions considering a state space of dimension to be defined by the user, i.e., by the value selected for K. This function is divided into two sections. The user is expected to only need to change the inputs in *Section A*.

Input parameters:

- ***SelectK***: define the number of PL states (K) to consider for further analyses;
- ***LEiDA_directory***: main directory of the LEiDA Toolbox; this should be a folder containing a folder (for example, designated LEiDA_Functions) which contains all the functions available in the toolbox;
- ***run_name***: name of the current run to be used to create a folder to save the results from the LEiDA analyses;

After filling the input parameters click *Run* on your Matlab dashboard.

Output:

- Computes the Transition Probability Matrix (TPM) for each participant and performs permutation tests to compre the state-to-state transition probabilities between conditions. The results are saved in the file *LEiDA_Stats_TransitionMatrix.mat* to the folder corresponding to selected K;
- Plots the mean TPM for each condition;
- Plots the intergroup differences between state-to-state transition probabilties;

---

### LEiDA State Time Courses

The function ***LEiDA_StateTime.m*** should be used to analyse the state time courses defined for a specified value of K across participants and for a given subject. This function is divided into two sections. The user is expected to only need to change the inputs in *Section A*.

Input parameters:

- ***SelectK***: define the number of PL states (K) to consider for further analyses;
- ***Subject***: define the subject whose state time course should be plotted. The user should provide the filename or an ID/number that uniquely identified the subject;
- ***LEiDA_directory***: main directory of the LEiDA Toolbox; this should be a folder containing a folder (for example, designated LEiDA_Functions) which contains all the functions available in the toolbox;
- ***run_name***: name of the current run to be used to create a folder to save the results from the LEiDA analyses;

After filling the input parameters click *Run* on your Matlab dashboard.

Output:

- Plots the state time courses considering K PL states for all participants and separates them by condition;
- Plots the state time course for a given subject as cluster blocks;
- Plots the state time course for a given subject as stairs plot;
### Other functions

- *LEiDA_diam_cluster*, *pdfc_diametrical_clustering*, *pdfc_diametrical_clustering_plusplus*: diametrical clustering, an alternative to k-means;
- *cluster_performance*, *cluster_stability*, *dunns*, *rand_index*: quality and stability of the clustering solutions;
- *EigenVectors_VoxelSpace*, *Plot_K_V1_VoxelSpace*, *Plot_K_V1_VoxelSpace_Slice*, *VIDEO_fMRI_Slices*: voxel-space visualisations (see also the LEiDA_Voxel repository);
- *Parcellate_ABIDE_func_preproc*: parcellation of the ABIDE dataset;
- *Plot_C_summary*, *Plot_C_vector_labelled*, *Plot_C_vector_ordered*, *Plot_C_boxplot_FO*, *Plot_K_vector_labelled*: previous single-state figures, replaced by *Plot_Centroid_Summary*.
