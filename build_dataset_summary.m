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
%       .mean_dist        dataset mean of mean interparticle distance
%       .sem_dist
%       .mean_corr_x      dataset mean of full-trace correlation x
%       .mean_corr_y
%       .mean_corr_r
%       .sem_corr_x/y/r
%
%   All means and SEMs are computed across movies (not across particles),
%   so they reflect movie-to-movie variability under the same condition.
%   Within a movie, per-particle values are first averaged across particles
%   before entering the cross-movie average.

    n_movies = numel(movie_results);
    summary.n_movies = n_movies;

    % Pre-allocate accumulators
    part_counts = zeros(n_movies, 1);
    per_movie_std_x   = NaN(n_movies, 1);
    per_movie_std_y   = NaN(n_movies, 1);
    per_movie_std_r   = NaN(n_movies, 1);
    per_movie_speed   = NaN(n_movies, 1);

    has_multi  = false;
    multi_dist   = [];
    multi_corr_x = [];
    multi_corr_y = [];
    multi_corr_r = [];

    for m = 1:n_movies
        mov = movie_results{m};
        nP  = mov.n_particles;
        part_counts(m) = nP;
        pm  = mov.particle_metrics;

        % Average particle statistics across particles within this movie
        std_x_vals = cellfun(@(p) p.std_x, pm);
        std_y_vals = cellfun(@(p) p.std_y, pm);
        std_r_vals = cellfun(@(p) p.std_r, pm);
        spd_vals   = cellfun(@(p) p.mean_speed_total, pm);

        per_movie_std_x(m) = nanmean(std_x_vals);
        per_movie_std_y(m) = nanmean(std_y_vals);
        per_movie_std_r(m) = nanmean(std_r_vals);
        per_movie_speed(m) = nanmean(spd_vals);

        % Pairwise metrics (only multi-particle movies)
        if nP > 1 && ~isempty(mov.pairwise_metrics)
            has_multi = true;
            pw = mov.pairwise_metrics;
            n_pairs = numel(pw.pairs);
            for k = 1:n_pairs
                multi_dist(end+1)   = pw.pairs{k}.mean_dist;    %#ok<AGROW>
                multi_corr_x(end+1) = pw.pairs{k}.corr_x;      %#ok<AGROW>
                multi_corr_y(end+1) = pw.pairs{k}.corr_y;      %#ok<AGROW>
                multi_corr_r(end+1) = pw.pairs{k}.corr_r;      %#ok<AGROW>
            end
        end
    end

    summary.particle_counts = part_counts;

    % Dataset-level single-particle stats
    summary.mean_std_x = nanmean(per_movie_std_x);
    summary.mean_std_y = nanmean(per_movie_std_y);
    summary.mean_std_r = nanmean(per_movie_std_r);
    summary.sem_std_x  = nanstd(per_movie_std_x) / sqrt(sum(~isnan(per_movie_std_x)));
    summary.sem_std_y  = nanstd(per_movie_std_y) / sqrt(sum(~isnan(per_movie_std_y)));
    summary.sem_std_r  = nanstd(per_movie_std_r) / sqrt(sum(~isnan(per_movie_std_r)));

    summary.mean_speed = nanmean(per_movie_speed);
    summary.sem_speed  = nanstd(per_movie_speed) / sqrt(sum(~isnan(per_movie_speed)));

    % Multi-particle summary (appended only if relevant)
    if has_multi
        mp.mean_dist  = nanmean(multi_dist);
        mp.sem_dist   = nanstd(multi_dist)   / sqrt(sum(~isnan(multi_dist)));
        mp.mean_corr_x = nanmean(multi_corr_x);
        mp.mean_corr_y = nanmean(multi_corr_y);
        mp.mean_corr_r = nanmean(multi_corr_r);
        mp.sem_corr_x  = nanstd(multi_corr_x) / sqrt(sum(~isnan(multi_corr_x)));
        mp.sem_corr_y  = nanstd(multi_corr_y) / sqrt(sum(~isnan(multi_corr_y)));
        mp.sem_corr_r  = nanstd(multi_corr_r) / sqrt(sum(~isnan(multi_corr_r)));
        summary.multi_particle = mp;
    end
end
