# Leading Eigenvector Dynamics Analysis Toolbox - Matlab

1. [Installation](#installation-of-leida-toolbox)
2. [LEiDA Start](#leida-start)
3. [LEiDA Analysis K](#leida-analysis-k)
4. [LEiDA Analysis Centroid](#leida-analysis-centroid)
5. [LEiDA Transitions K](#leida-transitions-k)
6. [LEiDA State Time Courses](#leida-state-time-courses)

---

## Installation of LEiDA Toolbox

Start by creating a new folder where all the LEiDA code and analyses will be saved. 
Then clone or download this repository into that folder. LEiDA_Start adds the folder with the toolbox functions to the Matlab path automatically.


> Note that the LEiDA Toolbox requires the Statistics and Machine Learning Matlab Toolbox.

---

## LEiDA Start

The function ***LEiDA_Start.m*** runs the first stage of the LEiDA analyses. The user is expected to only need to change the inputs in *Section A*. All parameters of the run are saved in *LEiDA_Parameters.mat* in the results folder.

> Requires MATLAB R2020b or later.

Input parameters:

*Directories*

- ***Data_directory***: folder with the input data. Files are searched in this folder and in all its subfolders (e.g. *sub-XXX/ses-YYY/*);
- ***Results_directory***: folder where the results are saved. A subfolder *LEiDA_Results_<run_name>* is created inside it. It does not need to be inside the toolbox folder;
- ***run_name***: name of the current run, used to name the results folder;

*Input data*

- ***Data_space***: `'parcellated'` if each file contains the time series of each brain area, or `'MNI'` if the files are 4D NIfTI images in MNI 2mm space, which are first parcellated with the atlas defined in *Parcellation* (see *Parcellate.m*);
- ***File_pattern***: pattern of the file names to use (`*` is a wildcard), e.g. `'*BOLDpreprocessed_rp_mo_csf_wm_nui-parcellatedAAL90.mat'` or `'*.nii.gz'`. Only files matching this pattern are used;
- ***Data_orientation***: `'AreasxTime'`, `'TimexAreas'` or `'auto'` (detect it from the number of areas);
- ***Data_variable***: name of the variable (or struct field) with the time series in .mat files; `''` to detect it automatically;
- ***Subject_list***: optional text file with the IDs of the subjects to include (one per line, e.g. `sub-control019`); `''` to use all files matching *File_pattern*;

*Parcellation*

- ***Parcellation***: atlas used to parcellate the data. To parcellate data in MNI space and to render the FC states in the brain it must be one of the atlases in *ParcelsMNI2mm.mat*: `'AAL116'`, `'AAL120'`, `'dbs80'`, `'Desikan'`, `'Schaefer100'`, `'glasser360'`. Set to `''` for other atlases (the brain renderings are skipped);
- ***N_areas***: number of areas to consider (the first *N_areas* areas of the data, e.g. 90 to use only the cerebrum areas of AAL116). Leave `[]` to use all the areas in the data (or all the parcels of the atlas for data in MNI space);

*Signal processing*

- ***TR***: TR of the fMRI data in seconds;
- ***apply_filter***: apply temporal filtering to the parcellated data (0: no; 1: yes);
- ***flp***: lowpass frequency of filter (default 0.1);
- ***fhi***: highpass frequency of filter (default 0.01);

*Clustering*

- ***rangeK***: range of numbers of FC states (K) to compute (default `2:20`);

*Study design*

- ***Condition_source***: how to assign each scan to a condition:
    - `'filename'`: a scan belongs to the condition whose tag in *Condition_tags* is found in its file path, e.g. `{'sub-control','sub-patient'}`;
    - `'metadata'`: the condition is read from a table (.csv or .xlsx), see below;
    - `'none'`: no conditions (only eigenvectors, clustering and centroids);
- ***Metadata_file***: table (.csv or .xlsx) with one row per subject, or per subject and session. The header is the first row containing *Metadata_subject_column*, so rows above it (e.g. units) are ignored. All sheets of a spreadsheet are read, and the name of each sheet is available as a column named `sheet`;
- ***Metadata_subject_column***: column with the subject IDs, as they appear in the folder/file names (e.g. `sub-control019`);
- ***Metadata_session_column***: column with the session IDs (e.g. `ses-midcycle`); `''` if there is only one scan per subject;
- ***Condition_column***: column with the condition of each scan (or `'sheet'`);
- ***Condition_tags***: conditions to compare (tags in the file names, or values of *Condition_column*). Scans not belonging to any of these conditions are left out of the statistics. With `'metadata'`, leave `{}` to use all values found;
- ***Condition_names***: names of the conditions shown in the figures, in the same order as *Condition_tags* (`{}` to use *Condition_tags*), e.g. `{'Controls','Interictal','Ictal'}` for `Condition_tags = {'sub-control','ses-interictal','ses-ictal'}`;

*Statistics and figures*

- ***Paired_tests***: experimental paradigm. 0 (unpaired): subjects in different conditions are not the same; the test statistic is Welch's t, which does not assume equal variances. 1 (paired): the same subjects are scanned in all conditions; the scans of each subject are paired by subject ID and the test statistic is the paired t. `'auto'`: for each pair of conditions, paired if the two conditions share subjects (e.g. interictal vs ictal sessions of the same patients) and unpaired otherwise (e.g. controls vs patients);
- ***Stats_level***: `'subject'` to average the scans of each subject within a condition, so that each subject counts once (recommended when subjects have several scans per condition), or `'scan'` to treat each scan as an independent observation;
- ***n_permutations***: number of permutations. The smallest possible p-value is 1/(n_permutations+1), so at least ~4200 permutations are needed to reach the strictest threshold 0.05/sum(rangeK) for K = 2:20 (e.g. 1000 for a first quick analysis, 10000 for final results);
- ***n_bootstraps***: number of bootstrap samples for the 95% confidence interval of the effect sizes (0 to skip, e.g. 1000 for final results);
- ***CortexDirection***: direction to plot the FC states/brain (`'SideView'` or `'TopView'`);
- ***Run_eigenvectors***, ***Run_clustering***: set to 0 to reuse the eigenvectors and clusters computed in a previous run (e.g. to test a new study design).

After filling the input parameters click *Run* on your Matlab dashboard.

Output:

- Creates a directory to store the results from the LEiDA analyses;
- (If *Data_space* is `'MNI'`) Parcellates the data and saves it in the folder *Parcellated_<Parcellation>* inside the results folder;
- Computes the leading eigenvectors of the data and saves them in the file *LEiDA_EigenVectors.mat*, together with *Scan_info*, a table with the file, subject and session of each scan;
- Clusters the leading eigenvectors of all scans and saves the clusters in the file *LEiDA_Clusters.mat*;
- Assigns each scan to a condition and saves it in *LEiDA_Conditions.mat* (with *Scan_info* including all the columns of the metadata table);
- Computes the fractional occupancy, performs hypothesis (permutation) tests to compare the mean fractional occupancy across conditions and saves the results in the file *LEiDA_Stats_FracOccup.mat*;
- Plots figure with two-sided p-values obtained from the comparison of the mean fractional occupancy, plots figures with barplots of the fractional occupancy values and saves the figures in png and fig format;
- Plots the PL states (centroids) obtained by clustering the leading eigenvectors and saves the figure in png and fig format;

---

## Statistics

Differences between each pair of conditions are tested for each FC state and each K with permutation tests (*permutation_test.m*, called through *LEiDA_compare_conditions.m*):

- **Unpaired** samples: Welch's t statistic; the null distribution is obtained by randomly permuting the condition labels;
- **Paired** samples: paired t statistic; the null distribution is obtained by randomly swapping the condition labels within each subject (sign flips of the differences);
- The p-values are counted directly from the permutation distribution, p = (1 + number of permuted statistics at least as extreme as the observed one) / (1 + n_permutations), without fitting the null distribution;
- Effect size: Hedges' g (unpaired) or g_z (paired), with a 95% bootstrap (percentile) confidence interval if *n_bootstraps* > 0;
- Significance levels shown in the figures (colors): p < 0.05 (dark red), p < 0.05/K (dark green), p < 0.05/sum(rangeK) (dark blue).
- Overlap with the Yeo et al. (2011) networks (*Overlap_LEiDA_Yeo.m*): each FC state is mapped to MNI space (voxels of the parcels with positive values in the centroid) and compared with each of the 7 cortical networks with the Dice coefficient, as in LEiDA_Voxel. A state is colored with the network of largest Dice if it is at least 0.3, otherwise in black.

---

## Format of the inputs

The input data should be contained in the folder specified in *Data_directory* (or its subfolders), with one file per scan. Subject and session IDs are read from BIDS-like paths (`sub-XXX`, `ses-YYY`); if the path has no `sub-XXX`, the file name is used as the subject ID. The time series data can be supplied in the following formats:

- A (no. of brain areas $\times$ no. of time points) or (no. of time points $\times$ no. of brain areas) matrix (see *Data_orientation*) for each scan. The file can be either a Matlab file (.mat), a text file (.txt) or a .1D file. If the .mat file contains several variables or a struct, set *Data_variable* to the name of the variable/field with the time series (by default, a variable or field named `data` is used);
- A 4D NIfTI file (.nii or .nii.gz) in MNI 2mm space for each scan, with *Data_space* set to `'MNI'`.

---

## LEiDA Analysis K

The function ***LEiDA_AnalysisK.m*** plots a number of relevant figures to analyse the LEiDA results for a specific number of clusters. Based on the based on the outputs from the previous analyses (*LEiDA_Start.m*) the user should select one or more values for K to explore with greater detail. This function is divided into two sections. The user is expected to only need to change the inputs in *Section A*.

Input parameters:

- ***SelectK***: define the number of PL states to consider for further analyses;
- ***LEiDA_directory***: main directory of the LEiDA Toolbox; this should be a folder containing a folder (for example, designated LEiDA_Functions) which contains all the functions available in the toolbox;
- ***run_name***: name of the current run to be used to create a folder to save the results from the LEiDA analyses;
- ***Parcellation***: parcellation template applied to the segment the fMRI data;

After filling the input parameters click *Run* on your Matlab dashboard.

Output:

- Plots several plots of the selected K:
    - PL states as links in cortex between the brain areas with positive values in the centroid vector;
    - PL states as nodes in cortex, where each brain area is coloured and scaled according to its value in the centroid vector;
    - PL states rendered in a 3D glass brain indicating the overlap with the functional brain networks defined by Yeo et al., (2011);
    - PL states in matrix format;
    - PL states in vector format with numbered areas;  
    - Boxplot of the fractional occupancy values of all K PL states;
    - PL states rendered in 3D glass brain, in vector format and boxplots of the fractional occupancy values;
    - Overlap between the K PL states and the 7 RSNs defined by Yeo et al., (2011);

The plots are saved in new directory created according to the value defined for K.

---

## LEiDA Analysis Centroid

The function ***LEiDA_AnalysisCentroid.m*** should be used to analyse a PL state for a given value of K selected by the user. This function plots a number of relevant figures to analyse this PL state. This function is divided into two sections. The user is expected to only need to change the inputs in *Section A*.

Input parameters:

- ***SelectK***: define the number of PL states (K) to consider for further analyses;
- ***Centroid***: define the PL state to consider for further analyses (1 <= Centroid <= SelectK);
- ***LEiDA_directory***: main directory of the LEiDA Toolbox; this should be a folder containing a folder (for example, designated LEiDA_Functions) which contains all the functions available in the toolbox;
- ***run_name***: name of the current run to be used to create a folder to save the results from the LEiDA analyses;
- ***Parcellation***: parcellation template applied to the segment the fMRI data;

After filling the input parameters click *Run* on your Matlab dashboard.

Output:

- Plots several plots of the selected PL state:
    - PL state in vector format with label of each brain area according to the specified parcellation;
    - PL state in vector format with brain areas organised by contribution;
    - Boxplot of the fractional occupancy values for the selected PL state;
    - PL state in vector format, rendered in a 3D glass brain, in matrix format and boxplots of the fractional occupancy values;

The plots are saved in new directory created according to the value of K and the specified PL state.

---

## LEiDA Transitions K

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

## LEiDA State Time Courses

The function ***LEiDA_TransitionsK.m*** should be used to analyse the state time courses defined for a specified value of K across participants and for a given subject. This function is divided into two sections. The user is expected to only need to change the inputs in *Section A*.

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