function print_assembly_summary(Results)
% PRINT_ASSEMBLY_SUMMARY  Print mean values for each movie to the command window.
%
%   print_assembly_summary(Results)
%
%   Prints a formatted table for every movie, followed by a dataset average.
%   Covers all extracted metrics: position spread, speed, interparticle
%   distance, motion correlations, angular speed, rotation frequency,
%   and synchronisation index.

    movies = Results.movies;
    meta   = Results.meta;
    n_movies = numel(movies);

    sep  = repmat('=', 1, 72);
    sep2 = repmat('-', 1, 72);

    fprintf('\n%s\n', sep);
    fprintf('  ASSEMBLY SUMMARY  |  %d movie(s)  |  lambda=%d nm  |  exp=%.4f s\n', ...
        n_movies, meta.wavelength_nm, meta.exposure_time_s);
    fprintf('%s\n\n', sep);

    for m = 1:n_movies
        mov  = movies{m};
        pm   = mov.particle_metrics;
        nP   = mov.n_particles;
        pw   = mov.pairwise_metrics;
        has_pw   = nP > 1 && ~isempty(fieldnames(pw));
        has_ang  = isfield(mov, 'angular_metrics');
        has_sync = isfield(mov, 'sync_metrics') && nP > 1;

        fprintf('%s\n', sep2);
        fprintf('  Movie %d   |   %d particle(s)\n', mov.movie_idx, nP);
        fprintf('%s\n', sep2);

        % --- Per-particle positional stats --------------------------------
        fprintf('  %-28s  %s\n', 'Metric', ...
            strjoin(arrayfun(@(p) sprintf('   P%-4d', p), 1:nP, 'UniformOutput',false), ''));

        prt_row(pm, nP, 'mean_x',           'Mean col         (px)');
        prt_row(pm, nP, 'mean_y',           'Mean row         (px)');
        prt_row(pm, nP, 'std_x',            'Std  col         (px)');
        prt_row(pm, nP, 'std_y',            'Std  row         (px)');
        prt_row(pm, nP, 'std_r',            'Std  r           (px)');
        prt_row(pm, nP, 'mean_speed_total', 'Mean speed       (px/s)');

        % --- Angular metrics (per particle) -------------------------------
        if has_ang
            % New:
            prt_row(pm, nP, 'mean_omega',       'Mean omega       (rad/s)');
            prt_row(pm, nP, 'mean_omega_abs',   'Mean |omega|     (rad/s)');
            prt_row(pm, nP, 'rot_freq_mean',    'Rot. freq (fit)  (Hz)');
            prt_row(pm, nP, 'rot_period_s',     'Rot. period      (s)');
            prt_row(pm, nP, 'dominant_freq_hz', 'FFT dom. freq    (Hz)');
        end

        % --- Pairwise metrics ---------------------------------------------
        if has_pw
            fprintf('\n  Pairwise\n');
            for k = 1:numel(pw.pairs)
                pr = pw.pairs{k};
                fprintf('    %-20s  mean dist = %8.3f px   std dist = %8.3f px\n', ...
                    pr.label, pr.mean_dist, pr.std_dist);
                fprintf('    %-20s  corr x = %+.3f   corr y = %+.3f   corr r = %+.3f\n', ...
                    '', pr.corr_x, pr.corr_y, pr.corr_r);
                fprintf('    %-20s  <roll corr x> = %+.3f   <roll corr y> = %+.3f\n', ...
                    '', pr.mean_roll_corr_x, pr.mean_roll_corr_y);
            end
        end

        % --- Synchronisation ----------------------------------------------
        if has_sync
            sm = mov.sync_metrics;
            fprintf('\n  Synchronisation\n');
            fprintf('    Kuramoto R           mean = %.4f   std = %.4f\n', ...
                sm.mean_R, sm.std_R);
            for k = 1:numel(sm.pairwise_labels)
                fprintf('    %-12s  PLI = %.4f   mean dphi = %+.4f rad   std dphi = %.4f rad\n', ...
                    sm.pairwise_labels{k}, sm.phase_locking_index(k), ...
                    sm.pairwise_mean_dphi(k), sm.pairwise_std_dphi(k));
            end
        end

        fprintf('\n');
    end

    % =====================================================================
    %  Dataset-level averages
    % =====================================================================
    s = Results.summary;
    fprintf('%s\n', sep);
    fprintf('  DATASET AVERAGES  (mean +/- SEM across %d movies)\n', n_movies);
    fprintf('%s\n', sep);
    fprintf('  Std col          %8.4f +/- %.4f px\n',  s.mean_std_x, s.sem_std_x);
    fprintf('  Std row          %8.4f +/- %.4f px\n',  s.mean_std_y, s.sem_std_y);
    fprintf('  Std r            %8.4f +/- %.4f px\n',  s.mean_std_r, s.sem_std_r);
    fprintf('  Mean speed       %8.4f +/- %.4f px/s\n',s.mean_speed, s.sem_speed);
    if isfield(s, 'multi_particle')
        mp = s.multi_particle;
        fprintf('  Mean dist        %8.4f +/- %.4f px\n',  mp.mean_dist,   mp.sem_dist);
        fprintf('  Corr x           %+8.4f +/- %.4f\n',    mp.mean_corr_x, mp.sem_corr_x);
        fprintf('  Corr y           %+8.4f +/- %.4f\n',    mp.mean_corr_y, mp.sem_corr_y);
        fprintf('  Corr r           %+8.4f +/- %.4f\n',    mp.mean_corr_r, mp.sem_corr_r);
    end
    fprintf('%s\n\n', sep);
end


% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function prt_row(pm, nP, field, label)
% Print one metric row for all particles.
    vals = cellfun(@(p) p.(field), pm(1:nP));
    val_str = sprintf('%8.4f ', vals);
    fprintf('  %-28s  %s\n', label, val_str);
end

