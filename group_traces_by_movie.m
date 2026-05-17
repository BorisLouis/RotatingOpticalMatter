function movies = group_traces_by_movie(TrackRes)
% GROUP_TRACES_BY_MOVIE  Split TrackRes.traces into per-movie structs.
%
%   movies = group_traces_by_movie(TrackRes)
%
%   INPUT
%   -----
%   TrackRes  - struct with field .traces  {N x 2} cell array
%                 column 1: [nFrames x 4] trace matrix
%                 column 2: scalar movie index
%
%   OUTPUT
%   ------
%   movies   - {M x 1} cell array.  Each element is a struct:
%                .movie_idx   scalar, original movie index
%                .traces      {nParticles x 1} cell array of [nFrames x 4]
%                             matrices. Only valid (non-empty, numeric) rows
%                             are included.
%
%   NOTES
%   -----
%   - Traces belonging to the same movie index are grouped in order of
%     appearance in TrackRes.traces.
%   - Empty or non-numeric trace entries are silently skipped here
%     (they were already warned about in load_tracking_data).
%   - The function does NOT assume movies are sorted by index.

    traces    = TrackRes.traces;
    n_rows    = size(traces, 1);

    % --- Collect valid movie indices ------------------------------------
    movie_indices = NaN(n_rows, 1);
    for k = 1:n_rows
        idx = traces{k, 2};
        if isnumeric(idx) && isscalar(idx) && ~isnan(idx)
            movie_indices(k) = idx;
        end
    end

    unique_movies = unique(movie_indices(~isnan(movie_indices)));
    n_movies = numel(unique_movies);

    % --- Build per-movie structs ----------------------------------------
    movies = cell(n_movies, 1);

    for m = 1:n_movies
        midx = unique_movies(m);
        rows = find(movie_indices == midx);

        % Collect only non-empty numeric traces for this movie
        valid_traces = {};
        for k = 1:numel(rows)
            r  = rows(k);
            tr = traces{r, 1};
            if ~isempty(tr) && (istable(tr) || (isnumeric(tr) && size(tr, 2) >= 4))
                valid_traces{end+1} = tr; %#ok<AGROW>
            end
        end

        if isempty(valid_traces)
            warning('OT:EmptyMovie', ...
                'Movie index %d has no valid traces. Skipping.', midx);
            continue;
        end

        movies{m} = struct( ...
            'movie_idx', midx, ...
            'traces',    {valid_traces} );
    end

    % Remove empty slots (movies that had no valid traces)
    movies = movies(~cellfun(@isempty, movies));

    fprintf('  Grouped into %d movie(s).\n', numel(movies));
end
