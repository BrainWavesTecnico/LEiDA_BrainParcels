function SelectSample345(save_dir)

load([save_dir 'LEiDA_EigenVectors.mat'],'Data_info','V1_all','idx_data','repaired','Time_sessions')   % V1_all, Time_sessions, Data_info, idx_data, repaired

% Subjects to keep
keep_subjs = strtrim(readlines('SC-FC_common_subjs.txt'));
keep_subjs = keep_subjs(keep_subjs ~= "");

% Subject ID of each scan, e.g. 'sub-002S0413_Desikan.mat' -> 'sub-002S0413'
scan_subjs = extractBefore(string({Data_info.name}), '_');

keep_scan = ismember(scan_subjs, keep_subjs);          % 1 x n_scans logical
disp(['Removing: ' strjoin(scan_subjs(~keep_scan), ', ')])

% Keep only eigenvectors from those scans
keep_t        = keep_scan(Time_sessions);
V1_all        = V1_all(keep_t,:);
Time_sessions = Time_sessions(keep_t);

% Renumber scans 1..345 so Time_sessions matches the new Data_info
new_idx       = cumsum(keep_scan);                      % old index -> new index
Time_sessions = new_idx(Time_sessions);
idx_data      = new_idx(idx_data(keep_scan(idx_data)));
if ~isempty(repaired)
    repaired = repaired(keep_scan([repaired{:,1}]),:);
    repaired(:,1) = num2cell(new_idx([repaired{:,1}]))';
end
Data_info = Data_info(keep_scan);

% Checks
assert(numel(Data_info) == numel(keep_subjs), 'Some SC-FC subjects are missing')
assert(isequal(unique(Time_sessions), 1:numel(Data_info)))

save([save_dir 'LEiDA_EigenVectors.mat'], 'V1_all','Time_sessions','Data_info','idx_data','repaired')