function [V1_all,Time_sessions,Data_info] = LEiDA_data(data_dir,file_pattern,save_dir,n_areas,n_parcels,orientation,data_variable,subject_list,filter,flp,fhi,tr)
%
% For each scan compute the leading eigenvector of the phase coherence
% matrix calculated at each recording frame.
%
% INPUT:
% data_dir      directory where the parcellated fMRI data are saved; files
%               are searched in this directory and all its subfolders
% file_pattern  pattern of the file names to use (* is a wildcard), e.g.
%               '*parcellatedAAL90.mat'; one file per scan; files can be
%               .mat (matrix, or struct with the time series in a field),
%               .txt or .1D
% save_dir      directory where the leading eigenvectors will be saved
% n_areas       number of brain areas to consider for analysis (the first
%               n_areas of the data); [] to use all areas
% n_parcels     number of parcels of the atlas ([] if unknown); only used
%               to detect the orientation of the data when 'auto'
% orientation   'AreasxTime', 'TimexAreas' or 'auto'
% data_variable name of the variable/field with the time series in .mat
%               files ('' to detect it automatically)
% subject_list  text file with the IDs of the subjects to include ('' to
%               include all scans)
% filter        0, no temporal filtering; 1, temporal filtering
% flp           lowpass frequency of filter
% fhi           highpass frequency of filter
% tr            TR of fMRI data
%
% OUTPUT:
% V1_all        (n_timepoints x n_areas) leading eigenvectors of all scans
%                at each time point
% Time_sessions (1 x n_timepoints) scan number of each leading eigenvector
% Data_info      file information of each scan (output of dir)
%
% Also saves Scan_info, a table with the file, subject and session of each
% scan (subject and session are read from BIDS-like paths, sub-XXX/ses-YYY)
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LEADING EIGENVECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ')

% Get the files containing the data in data_dir and all its subfolders
Data_info = dir(fullfile(data_dir, '**', file_pattern));
Data_info = Data_info(~[Data_info.isdir]);
if isempty(Data_info)
    error('LEiDA:nodata','No files matching ''%s'' were found in %s', file_pattern, data_dir);
end
% Sort by full path so that the order of the scans is reproducible
[~,ind_sort] = sort(fullfile({Data_info.folder},{Data_info.name}));
Data_info = Data_info(ind_sort);

% Subject and session of each scan
Scan_info = Get_scan_info(Data_info,data_dir);

% Keep only the subjects in subject_list
if ~isempty(subject_list)
    keep_subjs = strtrim(readlines(subject_list));
    keep_subjs = keep_subjs(keep_subjs ~= "");
    keep_scan = ismember(Scan_info.subject, keep_subjs);
    disp(['Subjects not in ' subject_list ' are excluded: ' num2str(nnz(~keep_scan)) ' scans removed'])
    missing_subjs = setdiff(keep_subjs, Scan_info.subject);
    if ~isempty(missing_subjs)
        disp(['Attention: no scans were found for ' strjoin(missing_subjs,', ')])
    end
    Data_info = Data_info(keep_scan);
    Scan_info = Scan_info(keep_scan,:);
end

% Total number of scans that will be read
n_scans = size(Data_info,1);
disp(['Total number of scans: ' num2str(n_scans)])
disp(' ')

% Leading eigenvectors of each scan
V1_scans = cell(n_scans,1);

discarded = [];
idx_data = [];
% Scans in which some areas had NaN/Inf or 0s and were repaired
% (columns: scan index, list of repaired areas)
repaired = {};
% Fixed seed so that the noise added to empty areas is reproducible
rng(0,'twister');
for s = 1:n_scans

    disp(['Computing the leading eigenvectors for scan ' num2str(s) ' ' Scan_info.file{s}]);

    signal = Load_signal(fullfile(Data_info(s).folder,Data_info(s).name),data_variable);

    % Put the signal in the format areas x time
    signal = Orient_signal(signal,orientation,n_areas,n_parcels,s == 1);

    % Selecting only the areas specified by the user
    if isempty(n_areas)
        n_areas = size(signal,1);
        disp(['      - Using all ' num2str(n_areas) ' areas in the data'])
    end
    if size(signal,1) < n_areas
        disp(['      - Only ' num2str(size(signal,1)) ' areas (expected ' num2str(n_areas) ') -> Discarded']);
        discarded = cat(2,discarded,s);
        continue
    end
    signal = signal(1:n_areas,:);

    % Instead of discarding scans with NaN, Inf or rows of 0s, repair them:
    % - areas with only some NaN/Inf time points are linearly interpolated
    % - areas that are entirely NaN/Inf, 0 or constant are replaced by
    %   low-amplitude Gaussian noise
    [signal,bad_areas] = RepairSignal(signal);
    if ~isempty(bad_areas)
        disp(['      - NaN, Inf or rows of 0s were found in areas ' mat2str(bad_areas) ' -> Repaired']);
        repaired = cat(1,repaired,{s,bad_areas});
    end

    if any(~isfinite(signal(:)))
        disp('      - Signal could not be repaired -> Discarded');
        discarded = cat(2,discarded,s);
        continue
    end

    % Add index of scan to idx_data
    idx_data = cat(2,idx_data,s);

    % De-meaning the fMRI signal
    for n = 1:size(signal,1)
        signal(n,:) = detrend(signal(n,:) - mean(signal(n,:)));
    end

    % Apply temporal filtering to the parcellated data
    if filter
        signal = TemporalFiltering(signal,flp,fhi,tr);
    end

    % Get the fMRI signal phase using the Hilbert transform
    for seed = 1:size(signal,1)
        signal(seed,:) = angle(hilbert(signal(seed,:)));
    end

    % Compute leading eigenvector of each phase coherence matrix
    % (exclude 1st and last TR of each fMRI signal)
    T = size(signal,2);
    V1_scans{s} = zeros(T-2,n_areas);
    for t = 2:T-1
        [v1,~] = eigs(cos(signal(:,t)-signal(:,t)'),1);
        V1_scans{s}(t-1,:) = v1;
    end
end

% Leading eigenvectors of all scans and the scan to which each belongs
V1_all = cat(1,V1_scans{:});
Time_sessions = repelem(1:n_scans, cellfun(@(v) size(v,1), V1_scans)');
clear V1_scans
N_areas = n_areas;

disp(' ')
disp(['Total number of scans used to compute the leading eigenvectors: ' num2str(length(idx_data))]);

for d = 1:length(discarded)
    disp(['Attention: ' Scan_info.file{discarded(d)} ' was discarded'])
end
for d = 1:size(repaired,1)
    disp(['Attention: ' Scan_info.file{repaired{d,1}} ' was repaired in areas ' mat2str(repaired{d,2})])
end
disp(' ')

% Name of the file to save output
save_file = 'LEiDA_EigenVectors.mat';

% Use the v7.3 format if the eigenvectors are too large for the v7 format
if numel(V1_all)*8 > 1.9e9
    save_format = '-v7.3';
else
    save_format = '-v7';
end
save(fullfile(save_dir,save_file), 'V1_all','Time_sessions','Data_info','Scan_info',...
     'idx_data','repaired','discarded','N_areas',save_format)

disp(['fMRI Phase Leading Eigenvectors saved successfully as ' save_file])
disp(' ')
end


function Scan_info = Get_scan_info(Data_info,data_dir)
%
% Table with the file (relative to data_dir), subject and session of each
% scan. Subject and session are read from BIDS-like paths (sub-XXX and
% ses-YYY). If the path has no sub-XXX, the file name is used as subject.
%

n_scans = numel(Data_info);
data_root = [regexprep(data_dir,'[\\/]+$','') filesep];
file = cell(n_scans,1);
subject = strings(n_scans,1);
session = strings(n_scans,1);
for s = 1:n_scans
    full_file = fullfile(Data_info(s).folder,Data_info(s).name);
    file{s} = erase(full_file, data_root);
    sub = regexp(file{s}, 'sub-[A-Za-z0-9]+', 'match', 'once');
    if isempty(sub)
        [~,sub] = fileparts(Data_info(s).name);
    end
    subject(s) = sub;
    session(s) = regexp(file{s}, 'ses-[A-Za-z0-9]+', 'match', 'once');
end
Scan_info = table(file,subject,session);
end


function signal = Load_signal(file,data_variable)
%
% Load the time series of one scan from a .mat, .txt or .1D file
%

[~,~,ext] = fileparts(file);
if strcmpi(ext,'.mat')
    S = load(file);
    names = fieldnames(S);
    if ~isempty(data_variable)
        if ~isfield(S,data_variable)
            error('LEiDA:variable','Variable ''%s'' not found in %s (variables: %s)', ...
                  data_variable, file, strjoin(names',', '));
        end
        signal = S.(data_variable);
    elseif numel(names) == 1
        signal = S.(names{1});
    elseif isfield(S,'data')
        signal = S.data;
    else
        % Use the largest numeric matrix in the file
        sizes = cellfun(@(n) isnumeric(S.(n))*numel(S.(n)), names);
        [~,ind] = max(sizes);
        signal = S.(names{ind});
    end
    % Time series inside a struct
    if isstruct(signal)
        if ~isempty(data_variable) && isfield(signal,data_variable)
            signal = signal.(data_variable);
        elseif isfield(signal,'data')
            signal = signal.data;
        else
            error('LEiDA:variable','Could not find the time series in %s; set Data_variable', file);
        end
    end
else
    % Text files (.txt, .1D, ...)
    signal = importdata(file);
    if isstruct(signal)
        signal = signal.data;
    end
end
signal = double(signal);
end


function signal = Orient_signal(signal,orientation,n_areas,n_parcels,verbose)
%
% Return the signal in the format areas x time
%

switch lower(orientation)
    case 'areasxtime'
        % nothing to do
    case 'timexareas'
        signal = signal';
    case 'auto'
        % The dimension matching the number of areas/parcels is the areas
        n_expected = unique([n_parcels n_areas]);
        rows_match = any(size(signal,1) == n_expected);
        cols_match = any(size(signal,2) == n_expected);
        if cols_match && ~rows_match
            signal = signal';
        elseif ~rows_match && verbose
            disp(['      - Could not detect the orientation of the data (size ' mat2str(size(signal)) ...
                  '); assuming AreasxTime. Set Data_orientation if this is not correct.'])
        end
    otherwise
        error('LEiDA:orientation','Data_orientation must be ''AreasxTime'', ''TimexAreas'' or ''auto''');
end
if verbose
    disp(['      - Data has ' num2str(size(signal,1)) ' areas and ' num2str(size(signal,2)) ' time points'])
end
end


function [signal,bad_areas] = RepairSignal(signal)
%
% Repair areas with NaN, Inf or no signal so that the scan can be kept.
%
% INPUT:
% signal        (n_areas x T) parcellated fMRI signal
%
% OUTPUT:
% signal        repaired signal
% bad_areas     indices of the areas that were repaired
%

signal(isinf(signal)) = NaN;

% Areas with no usable signal: all NaN, or constant (e.g. all 0s)
nan_rows   = any(isnan(signal),2);
empty_rows = all(isnan(signal),2) | max(signal,[],2) == min(signal,[],2);
% Areas with only some missing time points
partial_rows = nan_rows & ~empty_rows;

bad_areas = find(nan_rows | empty_rows)';
if isempty(bad_areas)
    return
end

% Interpolate isolated missing time points from neighbouring time points
for n = find(partial_rows)'
    signal(n,:) = fillmissing(signal(n,:),'linear','EndValues','nearest');
end

% Amplitude of the noise: a small fraction of the typical signal amplitude.
% The Hilbert phase does not depend on the amplitude, so this value only
% needs to be small and non-zero to avoid numerical problems.
good_std = std(signal(~empty_rows,:),0,2);
if isempty(good_std)
    noise_std = 1e-3;
else
    noise_std = 1e-3*median(good_std);
end
signal(empty_rows,:) = noise_std*randn(nnz(empty_rows),size(signal,2));
end
