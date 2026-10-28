function P = LEiDA_Start(P)
%
% First step of LEiDA (called from LEiDA_Main, where all the parameters
% are set): compute the leading eigenvectors of the phase coherence of all
% scans and cluster them into K FC states for K = rangeK.
%
% The eigenvectors (LEiDA_EigenVectors.mat) and clusters
% (LEiDA_Clusters.mat) are computed only if P.Run_eigenvectors or
% P.Run_clustering are 1, or if they are not found in the results folder;
% otherwise those of a previous run are reused.
%
% INPUT:
% P         parameters set in LEiDA_Main
%
% OUTPUT:
% P         parameters, with N_areas set to the number of areas used
%
% Saves the parameters of the run in LEiDA_Parameters.mat.
%
% Tutorial: README.md
% Version:  V2.0, October 2026
% Authors:  Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%           Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

if nargin == 0
    error('LEiDA:start','Set the parameters and run the analysis from LEiDA_Main');
end

% Create a directory to store the results from the current LEiDA run
if ~exist(P.leida_results, 'dir')
    mkdir(P.leida_results);
end

% Number of parcels in the atlas ([] if not available in ParcelsMNI2mm)
N_parcels = Atlas_Nparcels(P.Parcellation);
if isempty(P.N_areas) && strcmpi(P.Data_space,'MNI')
    P.N_areas = N_parcels;
end

%% LEADING EIGENVECTORS

file_V1 = [P.leida_results 'LEiDA_EigenVectors.mat'];
if P.Run_eigenvectors || ~isfile(file_V1)
    Data_directory = P.Data_directory;
    File_pattern = P.File_pattern;
    Data_orientation = P.Data_orientation;
    Data_variable = P.Data_variable;

    % Parcellate the data in MNI space with the selected atlas
    if strcmpi(P.Data_space,'MNI')
        parcellated_dir = [P.leida_results 'Parcellated_' P.Parcellation filesep];
        Parcellate(Data_directory,File_pattern,parcellated_dir,P.Parcellation,P.N_areas);
        Data_directory = parcellated_dir;
        File_pattern = ['*_' P.Parcellation '.mat'];
        Data_orientation = 'AreasxTime';
        Data_variable = 'data';
    end

    % Compute the leading eigenvectors of the data
    LEiDA_data(Data_directory,File_pattern,P.leida_results,P.N_areas,N_parcels,...
               Data_orientation,Data_variable,P.Subject_list,P.apply_filter,P.flp,P.fhi,P.TR);
    % The clusters must be computed again for the new eigenvectors
    P.Run_clustering = 1;
else
    disp(['Using the leading eigenvectors computed before (' file_V1 ')'])
end
% Number of areas used
S = load(file_V1,'N_areas');
P.N_areas = S.N_areas;

%% CLUSTERING

file_cluster = [P.leida_results 'LEiDA_Clusters.mat'];
if P.Run_clustering || ~isfile(file_cluster)
    LEiDA_cluster(P.leida_results,P.rangeK);
else
    S = load(file_cluster,'rangeK');
    if ~isequal(S.rangeK,P.rangeK)
        warning('LEiDA:start',['The clusters found were computed for K = %d to %d; set ' ...
                'Run_clustering = 1 to compute them for the rangeK set'],S.rangeK(1),S.rangeK(end));
        P.rangeK = S.rangeK;
    end
    disp(['Using the clusters computed before (' file_cluster ')'])
end

% Save the parameters of this run for reference
save([P.leida_results 'LEiDA_Parameters.mat'],'-struct','P');
