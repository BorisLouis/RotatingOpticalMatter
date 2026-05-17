function tables = build_summary_tables(movie_results, summary)
% BUILD_SUMMARY_TABLES  Assemble MATLAB tables for easy inspection and export.
%
%   tables = build_summary_tables(movie_results, summary)
%
%   OUTPUT
%   ------
%   tables  - struct with fields:
%
%     .movie_summary   - one row per movie, averaged across particles
%         MovieIndex | nParticles | MeanX | MeanY | StdX | StdY | StdR |
%         MeanSpeedTotal
%
%     .particle_detail - one row per (movie, particle) pair
%         MovieIndex | ParticleID | nFrames | MeanX | MeanY | StdX | StdY |
%         StdR | MeanSpeedTotal
%
%     .pairwise_summary - one row per (movie, pair) combination
%         MovieIndex | PairLabel | MeanDist | StdDist |
%         CorrX | CorrY | CorrR |
%         MeanRollCorrX | MeanRollCorrY | MeanRollCorrR
%         (only present if multi-particle movies exist)
%
%     .dataset_summary - one-row table of the dataset-level averages

    n_movies = numel(movie_results);

    % =====================================================================
    %  1.  Movie-level summary table
    % =====================================================================
    MovieIndex     = zeros(n_movies, 1);
    nParticles     = zeros(n_movies, 1);
    Mean_X         = NaN(n_movies, 1);
    Mean_Y         = NaN(n_movies, 1);
    Std_X          = NaN(n_movies, 1);
    Std_Y          = NaN(n_movies, 1);
    Std_R          = NaN(n_movies, 1);
    MeanSpeedTotal = NaN(n_movies, 1);
    MeanOmega      = NaN(n_movies, 1);
    MeanOmegaAbs   = NaN(n_movies, 1);
    RotFreq        = NaN(n_movies, 1);
    RotPeriod      = NaN(n_movies, 1);
    DomFreq        = NaN(n_movies, 1);
    
    for m = 1:n_movies
        mov = movie_results{m};
        pm  = mov.particle_metrics;
        MovieIndex(m)     = mov.movie_idx;
        nParticles(m)     = mov.n_particles;
        Mean_X(m)         = nanmean(cellfun(@(p) p.mean_x,           pm));
        Mean_Y(m)         = nanmean(cellfun(@(p) p.mean_y,           pm));
        Std_X(m)          = nanmean(cellfun(@(p) p.std_x,            pm));
        Std_Y(m)          = nanmean(cellfun(@(p) p.std_y,            pm));
        Std_R(m)          = nanmean(cellfun(@(p) p.std_r,            pm));
        MeanSpeedTotal(m) = nanmean(cellfun(@(p) p.mean_speed_total, pm));
        MeanOmega(m)      = nanmean(cellfun(@(p) p.mean_omega,       pm));
        MeanOmegaAbs(m)   = nanmean(cellfun(@(p) p.mean_omega_abs,   pm));
        RotFreq(m)        = nanmean(cellfun(@(p) p.rot_freq_mean,    pm));
        RotPeriod(m)      = nanmean(cellfun(@(p) p.rot_period_s,     pm));
        DomFreq(m)        = nanmean(cellfun(@(p) p.dominant_freq_hz, pm));
    end
    
    tables.movie_summary = table( ...
        MovieIndex, nParticles, Mean_X, Mean_Y, ...
        Std_X, Std_Y, Std_R, MeanSpeedTotal, ...
        MeanOmega, MeanOmegaAbs, RotFreq, RotPeriod, DomFreq, ...
        'VariableNames', {'MovieIndex','nParticles', ...
        'Mean_X','Mean_Y','Std_X','Std_Y','Std_R','MeanSpeedTotal', ...
        'MeanOmega','MeanOmegaAbs','RotFreq_Hz','RotPeriod_s','DomFreq_Hz'});
    % =====================================================================
    %  2.  Per-particle detail table
    % =====================================================================
    rows_mid  = [];
    rows_pid  = [];
    rows_nf   = [];
    rows_mx   = [];
    rows_my   = [];
    rows_sx   = [];
    rows_sy   = [];
    rows_sr   = [];
    rows_spd  = [];

    for m = 1:n_movies
        mov = movie_results{m};
        pm  = mov.particle_metrics;
        for p = 1:numel(pm)
            rows_mid(end+1) = mov.movie_idx;  %#ok<AGROW>
            rows_pid(end+1) = p;              %#ok<AGROW>
            rows_nf(end+1)  = numel(pm{p}.x);%#ok<AGROW>
            rows_mx(end+1)  = pm{p}.mean_x;  %#ok<AGROW>
            rows_my(end+1)  = pm{p}.mean_y;  %#ok<AGROW>
            rows_sx(end+1)  = pm{p}.std_x;   %#ok<AGROW>
            rows_sy(end+1)  = pm{p}.std_y;   %#ok<AGROW>
            rows_sr(end+1)  = pm{p}.std_r;   %#ok<AGROW>
            rows_spd(end+1) = pm{p}.mean_speed_total; %#ok<AGROW>
        end
    end

    tables.particle_detail = table( ...
        rows_mid(:), rows_pid(:), rows_nf(:), ...
        rows_mx(:), rows_my(:), rows_sx(:), rows_sy(:), rows_sr(:), rows_spd(:), ...
        'VariableNames', {'MovieIndex','ParticleID','nFrames', ...
        'Mean_X','Mean_Y','Std_X','Std_Y','Std_R','MeanSpeedTotal'});

    % =====================================================================
    %  3.  Pairwise summary table  (only if multi-particle movies exist)
    % =====================================================================
    has_pairs = false;
    pw_mid   = [];
    pw_lbl   = {};
    pw_md    = [];
    pw_sd    = [];
    pw_cx    = [];
    pw_cy    = [];
    pw_cr    = [];
    pw_rcx   = [];
    pw_rcy   = [];
    pw_rcr   = [];

    for m = 1:n_movies
        mov = movie_results{m};
        pm = mov.pairwise_metrics;
        if or(mov(1).n_particles < 2,isempty(fieldnames(pm)))
            continue;
        end
        has_pairs = true;
        pw = mov.pairwise_metrics;
        for k = 1:numel(pw.pairs)
            pr = pw.pairs{k};
            pw_mid(end+1)  = mov.movie_idx;           %#ok<AGROW>
            pw_lbl{end+1}  = pr.label;                %#ok<AGROW>
            pw_md(end+1)   = pr.mean_dist;            %#ok<AGROW>
            pw_sd(end+1)   = pr.std_dist;             %#ok<AGROW>
            pw_cx(end+1)   = pr.corr_x;               %#ok<AGROW>
            pw_cy(end+1)   = pr.corr_y;               %#ok<AGROW>
            pw_cr(end+1)   = pr.corr_r;               %#ok<AGROW>
            pw_rcx(end+1)  = pr.mean_roll_corr_x;     %#ok<AGROW>
            pw_rcy(end+1)  = pr.mean_roll_corr_y;     %#ok<AGROW>
            pw_rcr(end+1)  = pr.mean_roll_corr_r;     %#ok<AGROW>
        end
    end

    if has_pairs
        tables.pairwise_summary = table( ...
            pw_mid(:), pw_lbl(:), pw_md(:), pw_sd(:), ...
            pw_cx(:), pw_cy(:), pw_cr(:), ...
            pw_rcx(:), pw_rcy(:), pw_rcr(:), ...
            'VariableNames', { ...
            'MovieIndex','PairLabel','MeanDist','StdDist', ...
            'CorrX','CorrY','CorrR', ...
            'MeanRollCorrX','MeanRollCorrY','MeanRollCorrR'});
    end

    % =====================================================================
    %  4.  Dataset-level one-row summary table
    % =====================================================================
    ds_vars  = {'nMovies','MeanStdX','SEM_StdX','MeanStdY','SEM_StdY', ...
                'MeanStdR','SEM_StdR','MeanSpeed','SEM_Speed'};
    ds_vals  = {summary.n_movies, summary.mean_std_x, summary.sem_std_x, ...
                summary.mean_std_y, summary.sem_std_y, ...
                summary.mean_std_r, summary.sem_std_r, ...
                summary.mean_speed, summary.sem_speed};

    if isfield(summary, 'multi_particle')
        mp = summary.multi_particle;
        ds_vars = [ds_vars, {'MeanDist','SEM_Dist', ...
            'MeanCorrX','SEM_CorrX','MeanCorrY','SEM_CorrY', ...
            'MeanCorrR','SEM_CorrR'}];
        ds_vals = [ds_vals, {mp.mean_dist, mp.sem_dist, ...
            mp.mean_corr_x, mp.sem_corr_x, ...
            mp.mean_corr_y, mp.sem_corr_y, ...
            mp.mean_corr_r, mp.sem_corr_r}];
    end

    tables.dataset_summary = cell2table(ds_vals, 'VariableNames', ds_vars);
end
