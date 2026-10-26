function [V1_all,Time_sessions,Data_info] = LEiDA_data(data_dir,parcellation,Subcortical,save_dir,n_areas,tmax,filter,flp,fhi,tr)
%
% For each subject compute the leading eigenvetor of the phase coherence
% matrix calculated at each recording frame.
%
% INPUT:
% data_dir      directory where the parcellated fMRI data are saved;
%               fMRI data should be one file per subject; data can be a
%               matrix containing the time series in the formats .mat,
%               .1D and .txt; if data is as a struct then fMRI signal
%               should correspond to a field called data
% save_dir      directory where the leading eigenvectors will be saved
% n_areas       number of brain areas to consider for analysis
% tmax          maximum number of volumes across fMRI sessions
% filter        0, temporal filtering (default); 1, no temporal filtering
% flp           lowpass frequency of filter
% fhi           highpass frequency of filter
% tr            TR of fMRI data
%
% OUTPUT:
% V1_all        (n_scans*(tmax-2) x n_areas) leading eigenvectors of all
%                subjects at each time point
% Time_sessions (1 x n_scans*(tmax-2)) scan number of each leading
%                eigenvector
% Data_info      parcellated data
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt
%         Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LEADING EIGENVECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ')


% Get the files containg the data (.mat or .1D or .txt)
Data_info = dir(fullfile(data_dir, '**', ['*' parcellation '.mat']));
% Total number of scans that will be read
n_scans = size(Data_info,1);
disp(['Total number of scans in folder: ' num2str(n_scans)])
disp(' ')

% Matrix to store the leading eigenvectors of all subjects at each TR
if Subcortical==0
    V1_all = zeros(n_scans*(tmax-2),n_areas);
elseif Subcortical==1
    V1_all = zeros(n_scans*(tmax-2),n_areas+19);
end
% Row vector with the scan number of each leading eigenvector
Time_sessions = zeros(1,n_scans*(tmax-2));  

t_all = 0;
discarded = [];
idx_data = [];
% Scans in which some areas had NaN/Inf or 0s and were repaired
% (columns: scan index, list of repaired areas)
repaired = {};
% Fixed seed so that the noise added to empty areas is reproducible
rng(0,'twister');
for s = 1:n_scans
    
    disp(['Computing the leading eigenvectors for scan ' num2str(s) ' ' Data_info(s).name]);
    
    % Handling the data differently depending on the file extension:
    % Regardless of the extension we import the data
    signal = importdata(fullfile(Data_info(s).folder,Data_info(s).name));
    if strcmp(parcellation,'Schaefer100')
        % The data in Schaefer is TxN
        signal = signal';
    end
    % Selecting only the areas specified by the user
    signal = signal(1:n_areas,:);

    if Subcortical
        subcort_file=dir(fullfile(Data_info(s).folder,'*subcortical.mat'));
        subcort_signal = importdata(fullfile(Data_info(s).folder,subcort_file.name));
        signal=cat(1,signal,subcort_signal);
    end
    
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
    else
        % Add index of participant to idx_data
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
        for t = 2:size(signal,2)-1 % exclude 1st and last TR of each fMRI signal

            % Save the leading eigenvector for time t
            [v1,~] = eigs(cos(signal(:,t)-signal(:,t)'),1);

            t_all = t_all + 1; % time point in V1_all

            % row t_all correponds to the computed eigenvector at time t for subject s
            V1_all(t_all,:) = v1;

            % to which subject the computed leading eigenvector belongs
            Time_sessions(t_all) = s;
        end
    end
end
% Reduce size in case some scans have less TRs than tmax
% In this case, these lines of code will not result in changes
V1_all(t_all+1:end,:) = [];
Time_sessions(:,t_all+1:end) = [];

disp(' ')
disp(['Total number of scans used to compute the leading eigenvectors: ' num2str(length(idx_data))]);

for d = 1:length(discarded)
    disp(['Attention: ' Data_info(discarded(d)).name ' was discarded'])
end
for d = 1:size(repaired,1)
    disp(['Attention: ' Data_info(repaired{d,1}).name ' was repaired in areas ' mat2str(repaired{d,2})])
end
disp(' ')

% Name of the file to save output
save_file = 'LEiDA_EigenVectors.mat';

save([save_dir save_file], 'V1_all','Time_sessions','Data_info','idx_data','repaired')

disp(['fMRI Phase Leading Eigenvectors saved successfully as ' save_file])
disp(' ')
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