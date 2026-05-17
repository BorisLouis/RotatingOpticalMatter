function trackRes = load_tracking_data(folder, exposure_time_s)
% LOAD_TRACKING_DATA  Load and validate TrackRes.mat from the given folder.
%
%   TrackRes = load_tracking_data(folder)
%
%   INPUT
%   -----
%   folder   - full path to directory containing TrackRes.mat
%
%   OUTPUT
%   ------
%   TrackRes - struct loaded from file, guaranteed to have:
%                .traces   {N x 2} cell array
%                  column 1: [nFrames x 4] matrix  (x, y, z, t)
%                  column 2: scalar movie index
%
%   ERRORS
%   ------
%   Throws OT:FileNotFound if the .mat file is missing.
%   Throws OT:BadFormat   if expected fields / sizes are wrong.

    filepath = fullfile(folder, 'trackRes.mat');

    % --- Existence check -------------------------------------------------
    if ~isfile(filepath)
        error('OT:FileNotFound', ...
            'trackRes.mat not found in:\n  %s', folder);
    end

    % --- Load ------------------------------------------------------------
    fprintf('Loading %s ...\n', filepath);
    loaded = load(filepath, 'trackRes');

    if ~isfield(loaded, 'trackRes')
        error('OT:BadFormat', ...
            'TrackRes.mat does not contain a variable named ''TrackRes''.');
    end
    trackRes = loaded.trackRes;
    for k = 1:size(trackRes.traces, 1)
        trackRes.traces{k, 1}.t = trackRes.traces{k, 1}.t * exposure_time_s;
    end

    % --- Field check -----------------------------------------------------
    if ~isfield(trackRes, 'traces')
        error('OT:BadFormat', ...
            'TrackRes struct is missing the ''traces'' field.');
    end

    traces = trackRes.traces;

    if ~iscell(traces)
        error('OT:BadFormat', ...
            'TrackRes.traces must be a cell array.');
    end

    if size(traces, 2) < 2
        error('OT:BadFormat', ...
            'TrackRes.traces must have at least 2 columns (trace | movie_idx).');
    end

    n_rows = size(traces, 1);
    fprintf('  Found %d trace row(s) in TrackRes.traces.\n', n_rows);

  % --- Row-level sanity check (warn, do not error) --------------------
bad_rows = [];
for k = 1:n_rows
    tr = traces{k, 1};
    if isempty(tr)
        bad_rows(end+1) = k; %#ok<AGROW>
        continue;
    end
    if istable(tr)
        % Table format: check that required columns exist
        required_cols = {'col','row','z','t'};
        missing = required_cols(~ismember(required_cols, tr.Properties.VariableNames));
        if ~isempty(missing)
            warning('OT:MissingColumns', ...
                'Row %d: table is missing column(s): %s.', ...
                k, strjoin(missing, ', '));
        end
    elseif isnumeric(tr)
        % Numeric matrix fallback: check column count
        if size(tr, 2) < 4
            warning('OT:ShortTrace', ...
                'Row %d: matrix has only %d column(s); expected 4 (x,y,z,t).', ...
                k, size(tr, 2));
        end
    else
        bad_rows(end+1) = k; %#ok<AGROW>
    end
end

if ~isempty(bad_rows)
    warning('OT:EmptyTrace', ...
        '%d trace(s) are empty or unrecognised format (rows: %s). They will be skipped.', ...
        numel(bad_rows), num2str(bad_rows));
end

    fprintf('  Data loaded successfully.\n');
end
