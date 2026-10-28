function [Index_Conditions,Condition_tags,Scan_info] = LEiDA_conditions(data_dir,source,Condition_tags,metadata_file,subject_column,session_column,condition_column)
%
% Assign each scan to a condition, either from tags in the file names or
% from a metadata table (.csv, .txt, .xlsx, .xls).
%
% INPUT:
% data_dir          directory where the LEiDA results are saved
% source            'filename' or 'metadata'
% Condition_tags    'filename': tags to search in the file path of each
%                   scan, e.g. {'sub-control','sub-patient'}; a condition
%                   can group several tags, e.g.
%                   {{'ses-interictal','ses-preictal','ses-postictal'},'ses-ictal'}
%                   'metadata': values of condition_column to compare;
%                   {} to use all values found
% metadata_file     table with one row per subject (or subject and session)
% subject_column    column of the table with the subject IDs (as in the
%                   file paths, e.g. sub-control019)
% session_column    column of the table with the session IDs (as in the
%                   file paths, e.g. ses-midcycle); '' if one scan per
%                   subject
% condition_column  column of the table with the condition; in
%                   spreadsheets with several sheets, use 'sheet' for the
%                   name of the sheet
%
% OUTPUT:
% Index_Conditions  (1 x n_scans) condition of each scan (index in
%                   Condition_tags); 0 for scans not in any condition
% Condition_tags    (1 x n_conditions) name of each condition (tags of a
%                   condition with several tags joined by '+')
% Scan_info         table with the information of each scan, including
%                   all the columns of the metadata table
%
% The outputs are saved in LEiDA_Conditions.mat
%
% Author: Joana Cabral, University of Minho, joanacabral@med.uminho.pt

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ')

load(fullfile(data_dir,'LEiDA_EigenVectors.mat'),'Scan_info','idx_data');
n_scans = height(Scan_info);
Index_Conditions = zeros(1,n_scans);
if iscell(Condition_tags)
    Condition_tags = Condition_tags(:)';
else
    Condition_tags = cellstr(Condition_tags);
end

switch lower(source)

    case 'filename'
        for c = 1:numel(Condition_tags)
            % A scan belongs to condition c if its path contains any of its tags
            has_tag = contains(Scan_info.file, cellstr(Condition_tags{c}));
            if any(has_tag & Index_Conditions' > 0)
                error('LEiDA:conditions','Some files contain more than one condition tag, e.g. %s', ...
                      Scan_info.file{find(has_tag & Index_Conditions' > 0,1)});
            end
            Index_Conditions(has_tag) = c;
        end
        % Name of each condition (tags joined by '+')
        Condition_tags = cellfun(@(t) strjoin(cellstr(t),'+'), Condition_tags, 'UniformOutput', false);

    case 'metadata'
        [meta_names,meta_cells] = Read_metadata(metadata_file,subject_column);
        columns = {subject_column,session_column,condition_column};
        for col = columns(~cellfun(@isempty,columns))
            if ~ismember(col{1},meta_names)
                error('LEiDA:conditions','Column ''%s'' not found in %s. Columns: %s', col{1}, ...
                      metadata_file, strjoin(meta_names,', '));
            end
        end
        Metadata = @(col) meta_cells(:,strcmp(meta_names,col));
        Condition_tags = cellstr(Condition_tags);
        meta_subject = Cell2string(Metadata(subject_column));
        if ~isempty(session_column)
            meta_session = Cell2string(Metadata(session_column));
        end

        % Row of the metadata table corresponding to each scan
        meta_row = zeros(n_scans,1);
        for s = 1:n_scans
            rows = meta_subject == Scan_info.subject(s);
            if ~isempty(session_column)
                rows = rows & meta_session == Scan_info.session(s);
            end
            rows = find(rows);
            if numel(rows) == 1
                meta_row(s) = rows;
            elseif numel(rows) > 1
                disp(['Attention: ' Scan_info.file{s} ' matches ' num2str(numel(rows)) ...
                      ' rows of the metadata table; the first one is used'])
                meta_row(s) = rows(1);
            else
                disp(['Attention: ' Scan_info.file{s} ' was not found in the metadata table'])
            end
        end

        % Condition of each scan
        meta_condition = Cell2string(Metadata(condition_column));
        scan_condition = strings(n_scans,1);
        scan_condition(:) = missing;
        scan_condition(meta_row > 0) = meta_condition(meta_row(meta_row > 0));
        if isempty(Condition_tags)
            Condition_tags = cellstr(unique(scan_condition(~ismissing(scan_condition))))';
        end
        for c = 1:numel(Condition_tags)
            Index_Conditions(scan_condition == Condition_tags{c}) = c;
        end

        % Add the metadata of each scan to Scan_info
        meta_columns = ~ismember(meta_names,{subject_column,session_column});
        scan_cells = repmat({missing},n_scans,nnz(meta_columns));
        scan_cells(meta_row > 0,:) = meta_cells(meta_row(meta_row > 0),meta_columns);
        Scan_metadata = cell2table(scan_cells,'VariableNames',matlab.lang.makeUniqueStrings( ...
            meta_names(meta_columns),Scan_info.Properties.VariableNames));
        Scan_info = [Scan_info Scan_metadata];

    otherwise
        error('LEiDA:conditions','Condition_source must be ''filename'', ''metadata'' or ''none''');
end

if isempty(Condition_tags)
    error('LEiDA:conditions','No conditions defined');
end

% Summary of the scans in each condition
disp('Number of scans in each condition (excluding discarded scans):')
for c = 1:numel(Condition_tags)
    disp(['- ' Condition_tags{c} ': ' num2str(nnz(Index_Conditions(idx_data) == c))])
end
n_none = nnz(Index_Conditions(idx_data) == 0);
if n_none
    disp(['- not in any condition (excluded from the statistics): ' num2str(n_none)])
end
for c = 1:numel(Condition_tags)
    if nnz(Index_Conditions(idx_data) == c) < 2
        error('LEiDA:conditions','Condition ''%s'' has less than 2 scans; check Condition_tags', Condition_tags{c});
    end
end
disp(' ')

save(fullfile(data_dir,'LEiDA_Conditions.mat'),'Index_Conditions','Condition_tags','Scan_info')
end


function [all_names,merged] = Read_metadata(file,subject_column)
%
% Read a metadata table (.csv, .txt, .xlsx, .xls). The header is the first
% row containing subject_column, so rows above the header (e.g. units) are
% ignored. All sheets of a spreadsheet are read and merged, and the name of
% each sheet is added in a column named 'sheet'. Returns the column names
% and a cell array with the values as read (one row per subject/session).
%

[~,~,ext] = fileparts(file);
if any(strcmpi(ext,{'.xlsx','.xls','.xlsm','.ods'}))
    sheets = cellstr(sheetnames(file));
else
    sheets = {''};
end

all_names = {};
blocks = {};
for i = 1:numel(sheets)
    if isempty(sheets{i})
        C = readcell(file);
    else
        C = readcell(file,'Sheet',sheets{i});
    end

    % Header row: first row containing subject_column
    is_header = cellfun(@(x) (ischar(x) || isstring(x)) && strcmp(strtrim(string(x)),subject_column), C);
    [r,~] = find(is_header);
    if isempty(r)
        disp(['Attention: column ''' subject_column ''' not found in sheet ''' sheets{i} ''' of ' file ' -> sheet ignored'])
        continue
    end
    r = min(r);

    % Columns with a name
    has_name = cellfun(@(x) (ischar(x) || isstring(x)) && strlength(strtrim(string(x))) > 0, C(r,:));
    names = cellfun(@(x) char(strtrim(string(x))), C(r,has_name), 'UniformOutput', false);
    names = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(names));
    data = C(r+1:end,has_name);

    % Rows with a subject ID
    sub_col = strcmp(names,subject_column);
    has_subject = cellfun(@(x) (ischar(x) || isstring(x)) && strlength(strtrim(string(x))) > 0, data(:,sub_col));
    data = data(has_subject,:);

    if ~isempty(sheets{i})
        names = [names {'sheet'}]; %#ok<AGROW>
        data(:,end+1) = {sheets{i}}; %#ok<AGROW>
    end
    all_names = [all_names setdiff(names,all_names,'stable')]; %#ok<AGROW>
    blocks{end+1} = {names,data}; %#ok<AGROW>
end
if isempty(blocks)
    error('LEiDA:metadata','Column ''%s'' not found in %s', subject_column, file);
end

% Merge all sheets, filling the columns missing in a sheet with missing
merged = cell(0,numel(all_names));
for i = 1:numel(blocks)
    [names,data] = blocks{i}{:};
    block = repmat({missing},size(data,1),numel(all_names));
    [~,cols] = ismember(names,all_names);
    block(:,cols) = data;
    merged = [merged; block]; %#ok<AGROW>
end
end


function str = Cell2string(values)
%
% Convert a cell array column with mixed values to a string array, with
% leading/trailing spaces removed
%

str = strings(numel(values),1);
for i = 1:numel(values)
    v = values{i};
    if (isnumeric(v) || islogical(v)) && isscalar(v)
        str(i) = string(v);
    elseif ischar(v) || isstring(v)
        str(i) = strtrim(string(v));
    elseif isdatetime(v)
        str(i) = string(v);
    else
        str(i) = missing;
    end
end
end
