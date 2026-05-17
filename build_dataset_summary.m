function summary = build_dataset_summary(movie_results)
% BUILD_DATASET_SUMMARY  Aggregate statistics across all movies.
%
%   summary = build_dataset_summary(movie_results)
%
%   INPUT
%   -----
%   movie_results  - {nMovies x 1} cell of per-movie structs
%                    (output of the main per-movie analysis loop)
%
%   OUTPUT
%   ------
%   summary  - struct with fields:
%
%     .n_movies           total number of movies analysed
%     .particle_counts    vector of n_particles per movie
%
%     .mean_std_x         dataset mean of per-movie (mean particle) std_x
%     .mean_std_y         dataset mean of per-movie (mean particle) std_y
%     .mean_std_r         dataset mean of per-movie (mean particle) std_r
%     .sem_std_x/y/r      standard error of the mean across movies
%
%     .mean_speed         dataset mean of mean_speed_total (per particle)
%     .sem_speed
%
%     .multi_particle     (only if any multi-particle movies exist)
%       .mean_dist            dataset mean of mean interparticle distance
%       .sem_dist
%
%       .mean_corr_x          dataset mean of full-trace x displacement corr
%       .mean_corr_y
%       .mean_corr_r
%       .mean_corr_theta      dataset mean of angular displacement corr
%                             (dtheta_i vs dtheta_j around CoM)
%       .sem_corr_x/y/r/theta
%
%       .mean_std_phi_circ    dataset mean of circular std of bond angle
%                             Primary orientational stability descriptor:
%                             0 = perfectly locked, sqrt(2) = fully random
%       .sem_std_phi_circ
%       .mean_std_phi_unwrap  dataset mean of std of unwrapped bond angle
%                             Sensitive to slow rotational drift
%       .sem_std_phi_unwrap
%
%   AGGREGATION STRATEGY
%   ---------------------
%   All means and SEMs are computed across pair observations pooled from
%   all multi-particle movies (one entry per pair per movie), consistent
%   with the existing approach for dist and corr_x/y/r.  This weights
%   each pair equally regardless of the number of particles in the movie.
%
%   std_phi_circ and std_phi_unwrap are averaged as scalars (not
%   mean_phi itself, whose circular average across movies is less
%   interpretable without a well-defined reference frame).

    n_movies = numel(movie_results);
    summary.n_movies = n_movies;

    % =====================================================================
    %  Pre-allocate per-movie accumulators
    % =====================================================================
    part_counts       = zeros(n_movies, 1);
    per_movie_std_x   = NaN(n_movies, 1);
    per_movie_std_y   = NaN(n_movies, 1);
    per_movie_std_r   = NaN(n_movies, 1);
    per_movie_speed   = NaN(n_movies, 1);

    % Pairwise accumulators (one entry per pair per movie, pooled)
    has_multi        = false;
    multi_dist       = [];
    multi_corr_x     = [];
    multi_corr_y     = [];
    multi_corr_r     = [];
    multi_corr_theta = [];      % NEW: angular displacement correlation
    multi_sphi_circ  = [];      % NEW: circular std of bond angle
    multi_sphi_unwrap = [];     % NEW: std of unwrapped bond angle

    % =====================================================================
    %  Per-movie loop
    % =====================================================================
    for m = 1:n_movies
        mov = movie_results{m};
        nP  = mov.n_particles;
        part_counts(m) = nP;
        pm  = mov.particle_metrics;

        % Average particle statistics across particles within this movie
        per_movie_std_x(m) = nanmean(cellfun(@(p) p.std_x,            pm));
        per_movie_std_y(m) = nanmean(cellfun(@(p) p.std_y,            pm));
        per_movie_std_r(m) = nanmean(cellfun(@(p) p.std_r,            pm));
        per_movie_speed(m) = nanmean(cellfun(@(p) p.mean_speed_total,  pm));

        % Pairwise metrics (only multi-particle movies)
        if nP < 2 || isempty(fieldnames(mov.pairwise_metrics))
            continue;
        end
        has_multi = true;
        pw = mov.pairwise_metrics;

        for k = 1:numel(pw.pairs)
            pr = pw.pairs{k};
            multi_dist(end+1)        = pr.mean_dist;        %#ok<AGROW>
            multi_corr_x(end+1)      = pr.corr_x;           %#ok<AGROW>
            multi_corr_y(end+1)      = pr.corr_y;           %#ok<AGROW>
            multi_corr_r(end+1)      = pr.corr_r;           %#ok<AGROW>
            multi_corr_theta(end+1)  = pr.corr_theta;       %#ok<AGROW>
            multi_sphi_circ(end+1)   = pr.std_phi_circ;     %#ok<AGROW>
            multi_sphi_unwrap(end+1) = pr.std_phi_unwrap;   %#ok<AGROW>
        end
    end

    % =====================================================================
    %  Dataset-level single-particle stats
    % =====================================================================
    summary.particle_counts = part_counts;

    summary.mean_std_x = nanmean(per_movie_std_x);
    summary.mean_std_y = nanmean(per_movie_std_y);
    summary.mean_std_r = nanmean(per_movie_std_r);
    summary.sem_std_x  = sem(per_movie_std_x);
    summary.sem_std_y  = sem(per_movie_std_y);
    summary.sem_std_r  = sem(per_movie_std_r);

    summary.mean_speed = nanmean(per_movie_speed);
    summary.sem_speed  = sem(per_movie_speed);

    % =====================================================================
    %  Multi-particle summary (appended only if relevant)
    % =====================================================================
    if has_multi
        mp.mean_dist   = nanmean(multi_dist);
        mp.sem_dist    = sem(multi_dist);

        mp.mean_corr_x = nanmean(multi_corr_x);
        mp.sem_corr_x  = sem(multi_corr_x);
        mp.mean_corr_y = nanmean(multi_corr_y);
        mp.sem_corr_y  = sem(multi_corr_y);
        mp.mean_corr_r = nanmean(multi_corr_r);
        mp.sem_corr_r  = sem(multi_corr_r);

        % Angular displacement correlation
        mp.mean_corr_theta = nanmean(multi_corr_theta);
        mp.sem_corr_theta  = sem(multi_corr_theta);

        % Bond angle stability
        mp.mean_std_phi_circ   = nanmean(multi_sphi_circ);
        mp.sem_std_phi_circ    = sem(multi_sphi_circ);
        mp.mean_std_phi_unwrap = nanmean(multi_sphi_unwrap);
        mp.sem_std_phi_unwrap  = sem(multi_sphi_unwrap);

        summary.multi_particle = mp;
    end
end


% =========================================================================
%  LOCAL UTILITY
% =========================================================================

function s = sem(x)
% Standard error of the mean, NaN-aware.
    x   = x(~isnan(x));
    n   = numel(x);
    if n < 2
        s = NaN;
    else
        s = std(x) / sqrt(n);
    end
end
