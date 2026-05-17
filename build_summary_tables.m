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
%         MovieIndex | nParticles | Mean_X | Mean_Y | Std_X | Std_Y | Std_R |
%         MeanSpeedTotal | MeanOmega | MeanOmegaAbs | RotFreq_Hz |
%         RotPeriod_s | DomFreq_Hz
%
%     .particle_detail - one row per (movie, particle) pair
%         MovieIndex | ParticleID | nFrames |
%         Mean_X | Mean_Y | Std_X | Std_Y | Std_R | MeanSpeedTotal |
%         MeanOmega | MeanOmegaAbs | RotFreq_Hz | RotPeriod_s | DomFreq_Hz
%
%     .pairwise_summary - one row per (movie, pair) combination
%         MovieIndex | PairLabel |
%         MeanDist | StdDist |
%         MeanPhi | StdPhi_Circ | StdPhi_Unwrap |
%         CorrX | CorrY | CorrR | CorrTheta |
%         MeanRollCorrX | MeanRollCorrY | MeanRollCorrR | MeanRollCorrTheta
%         (only present if multi-particle movies exist)
%
%     .dataset_summary - one-row table of dataset-level averages
%
%   BOND ANGLE COLUMNS (pairwise_summary)
%   --------------------------------------
%   MeanPhi       : circular mean of the bond orientation angle [rad]
%                   angle(mean(exp(i*phi))) — the preferred direction of
%                   the i->j bond in the lab frame
%   StdPhi_Circ   : circular std [rad], sqrt(-2*ln(R_bar))
%                   0 = perfectly locked orientation; sqrt(2) = fully random
%                   Primary descriptor for assembly orientational stability
%   StdPhi_Unwrap : std of the unwrapped bond angle trace [rad]
%                   Sensitive to slow rotational drift of the whole assembly
%
%   THETA CORRELATION COLUMNS (pairwise_summary)
%   ---------------------------------------------
%   CorrTheta         : full-trace Pearson corr of angular displacements
%                       d(theta_i) vs d(theta_j) around the CoM
%                       +1 = co-rotation, -1 = counter-rotation
%   MeanRollCorrTheta : time-averaged rolling theta correlation

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
    %      Includes all per-particle scalar metrics: positional, speed,
    %      and angular (merged into particle_metrics by compute_angular_metrics)
    % =====================================================================
    rows_mid   = [];
    rows_pid   = [];
    rows_nf    = [];
    rows_mx    = [];
    rows_my    = [];
    rows_sx    = [];
    rows_sy    = [];
    rows_sr    = [];
    rows_spd   = [];
    rows_om    = [];
    rows_omabs = [];
    rows_rf    = [];
    rows_rp    = [];
    rows_df    = [];

    for m = 1:n_movies
        mov = movie_results{m};
        pm  = mov.particle_metrics;
        for p = 1:numel(pm)
            rows_mid(end+1)   = mov.movie_idx;              %#ok<AGROW>
            rows_pid(end+1)   = p;                          %#ok<AGROW>
            rows_nf(end+1)    = numel(pm{p}.x);            %#ok<AGROW>
            rows_mx(end+1)    = pm{p}.mean_x;              %#ok<AGROW>
            rows_my(end+1)    = pm{p}.mean_y;              %#ok<AGROW>
            rows_sx(end+1)    = pm{p}.std_x;               %#ok<AGROW>
            rows_sy(end+1)    = pm{p}.std_y;               %#ok<AGROW>
            rows_sr(end+1)    = pm{p}.std_r;               %#ok<AGROW>
            rows_spd(end+1)   = pm{p}.mean_speed_total;    %#ok<AGROW>
            rows_om(end+1)    = pm{p}.mean_omega;          %#ok<AGROW>
            rows_omabs(end+1) = pm{p}.mean_omega_abs;      %#ok<AGROW>
            rows_rf(end+1)    = pm{p}.rot_freq_mean;       %#ok<AGROW>
            rows_rp(end+1)    = pm{p}.rot_period_s;        %#ok<AGROW>
            rows_df(end+1)    = pm{p}.dominant_freq_hz;    %#ok<AGROW>
        end
    end

    tables.particle_detail = table( ...
        rows_mid(:), rows_pid(:), rows_nf(:), ...
        rows_mx(:), rows_my(:), rows_sx(:), rows_sy(:), rows_sr(:), ...
        rows_spd(:), rows_om(:), rows_omabs(:), ...
        rows_rf(:), rows_rp(:), rows_df(:), ...
        'VariableNames', { ...
        'MovieIndex','ParticleID','nFrames', ...
        'Mean_X','Mean_Y','Std_X','Std_Y','Std_R', ...
        'MeanSpeedTotal','MeanOmega','MeanOmegaAbs', ...
        'RotFreq_Hz','RotPeriod_s','DomFreq_Hz'});

    % =====================================================================
    %  3.  Pairwise summary table  (only if multi-particle movies exist)
    %
    %  Columns added vs previous version:
    %    MeanPhi          circular mean bond orientation angle [rad]
    %    StdPhi_Circ      circular std of bond angle [rad] — stability metric
    %    StdPhi_Unwrap    std of unwrapped bond angle [rad] — drift metric
    %    CorrTheta        full-trace Pearson corr of angular displacements
    %    MeanRollCorrTheta time-averaged rolling theta correlation
    % =====================================================================
    has_pairs = false;

    pw_mid    = [];
    pw_lbl    = {};
    pw_md     = [];
    pw_sd     = [];
    pw_mphi   = [];
    pw_sphic  = [];
    pw_sphiuw = [];
    pw_cx     = [];
    pw_cy     = [];
    pw_cr     = [];
    pw_cth    = [];
    pw_rcx    = [];
    pw_rcy    = [];
    pw_rcr    = [];
    pw_rcth   = [];

    for m = 1:n_movies
        mov = movie_results{m};
        pw  = mov.pairwise_metrics;
        if mov.n_particles < 2 || isempty(fieldnames(pw))
            continue;
        end
        has_pairs = true;
        for k = 1:numel(pw.pairs)
            pr = pw.pairs{k};
            pw_mid(end+1)    = mov.movie_idx;           %#ok<AGROW>
            pw_lbl{end+1}    = pr.label;                %#ok<AGROW>
            pw_md(end+1)     = pr.mean_dist;            %#ok<AGROW>
            pw_sd(end+1)     = pr.std_dist;             %#ok<AGROW>
            pw_mphi(end+1)   = pr.mean_phi;             %#ok<AGROW>
            pw_sphic(end+1)  = pr.std_phi_circ;         %#ok<AGROW>
            pw_sphiuw(end+1) = pr.std_phi_unwrap;       %#ok<AGROW>
            pw_cx(end+1)     = pr.corr_x;               %#ok<AGROW>
            pw_cy(end+1)     = pr.corr_y;               %#ok<AGROW>
            pw_cr(end+1)     = pr.corr_r;               %#ok<AGROW>
            pw_cth(end+1)    = pr.corr_theta;           %#ok<AGROW>
            pw_rcx(end+1)    = pr.mean_roll_corr_x;     %#ok<AGROW>
            pw_rcy(end+1)    = pr.mean_roll_corr_y;     %#ok<AGROW>
            pw_rcr(end+1)    = pr.mean_roll_corr_r;     %#ok<AGROW>
            pw_rcth(end+1)   = pr.mean_roll_corr_theta; %#ok<AGROW>
        end
    end

    if has_pairs
        tables.pairwise_summary = table( ...
            pw_mid(:), pw_lbl(:), ...
            pw_md(:), pw_sd(:), ...
            pw_mphi(:), pw_sphic(:), pw_sphiuw(:), ...
            pw_cx(:), pw_cy(:), pw_cr(:), pw_cth(:), ...
            pw_rcx(:), pw_rcy(:), pw_rcr(:), pw_rcth(:), ...
            'VariableNames', { ...
            'MovieIndex','PairLabel', ...
            'MeanDist','StdDist', ...
            'MeanPhi','StdPhi_Circ','StdPhi_Unwrap', ...
            'CorrX','CorrY','CorrR','CorrTheta', ...
            'MeanRollCorrX','MeanRollCorrY','MeanRollCorrR','MeanRollCorrTheta'});
    end

    % =====================================================================
    %  4.  Dataset-level one-row summary table
    %
    %  The base fields come directly from the summary struct (computed
    %  upstream, e.g. in compute_dataset_summary).  Multi-particle fields
    %  are added conditionally; new bond-angle and theta-correlation fields
    %  are added if present in summary.multi_particle.
    % =====================================================================
    ds_vars = {'nMovies', ...
               'MeanStdX','SEM_StdX', ...
               'MeanStdY','SEM_StdY', ...
               'MeanStdR','SEM_StdR', ...
               'MeanSpeed','SEM_Speed'};
    ds_vals = {summary.n_movies, ...
               summary.mean_std_x, summary.sem_std_x, ...
               summary.mean_std_y, summary.sem_std_y, ...
               summary.mean_std_r, summary.sem_std_r, ...
               summary.mean_speed, summary.sem_speed};

    if isfield(summary, 'multi_particle')
        mp = summary.multi_particle;

        % --- Core pairwise fields (always present in multi_particle) -----
        ds_vars = [ds_vars, { ...
            'MeanDist','SEM_Dist', ...
            'MeanCorrX','SEM_CorrX', ...
            'MeanCorrY','SEM_CorrY', ...
            'MeanCorrR','SEM_CorrR'}];
        ds_vals = [ds_vals, { ...
            mp.mean_dist, mp.sem_dist, ...
            mp.mean_corr_x, mp.sem_corr_x, ...
            mp.mean_corr_y, mp.sem_corr_y, ...
            mp.mean_corr_r, mp.sem_corr_r}];

        % --- Bond angle fields (added by updated compute_dataset_summary) -
        if isfield(mp, 'mean_std_phi_circ')
            ds_vars = [ds_vars, { ...
                'MeanStdPhi_Circ','SEM_StdPhi_Circ', ...
                'MeanStdPhi_Unwrap','SEM_StdPhi_Unwrap'}];
            ds_vals = [ds_vals, { ...
                mp.mean_std_phi_circ,   mp.sem_std_phi_circ, ...
                mp.mean_std_phi_unwrap, mp.sem_std_phi_unwrap}];
        end

        % --- Theta correlation fields -------------------------------------
        if isfield(mp, 'mean_corr_theta')
            ds_vars = [ds_vars, {'MeanCorrTheta','SEM_CorrTheta'}];
            ds_vals = [ds_vals, {mp.mean_corr_theta, mp.sem_corr_theta}];
        end
    end

    tables.dataset_summary = cell2table(ds_vals, 'VariableNames', ds_vars);
end
