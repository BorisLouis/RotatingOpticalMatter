function plot_angular_results(Results, movie_number)
% PLOT_ANGULAR_RESULTS  Diagnostic figure for angular and synchronisation metrics.
%
%   plot_angular_results(Results, movie_number)
%
%   Requires that Results.movies{m} contains .angular_metrics and
%   .sync_metrics fields (added after calling compute_angular_metrics and
%   compute_synchronisation_index in the main workflow).
%
%   Panels
%   ------
%   1. Unwrapped theta(t) for all particles  → reveals net rotation
%   2. Instantaneous angular speed omega(t)  (signed rad/s)
%   3. Rolling mean omega(t)                 → local rotation trend
%   4. FFT power spectrum of omega           → oscillation frequencies
%   5. Kuramoto R(t) + rolling mean          → instantaneous synchrony
%   6. Pairwise phase difference dphi(t)     → phase locking

    mov = Results.movies{movie_number};

    if ~isfield(mov, 'angular_metrics')
        error('OT:MissingField', ...
            'angular_metrics not found. Run compute_angular_metrics first.');
    end

    am  = mov.angular_metrics;
    nP  = mov.n_particles;
    pm  = mov.particle_metrics;
    colors = lines(nP);

    has_sync = isfield(mov, 'sync_metrics') && nP > 1;
    n_panels = 4 + 2*has_sync;

    figure('Name', sprintf('Angular analysis — Movie %d (%d particle(s))', ...
        mov.movie_idx, nP), ...
        'NumberTitle', 'off', 'Position', [50 50 1400 650]);

    % --- Panel 1: unwrapped theta ----------------------------------------
    ax1 = subplot(2, ceil(n_panels/2), 1);  hold(ax1,'on');
    for p = 1:nP
        t = pm(p).t;
        plot(ax1, t, am(p).theta_unwrap, '-', 'Color', colors(p,:), ...
            'LineWidth', 1, 'DisplayName', sprintf('P%d', p));
    end
    xlabel(ax1,'t (s)');  ylabel(ax1,'\theta_{unwrap} (rad)');
    title(ax1,'Cumulative angle');  legend(ax1,'show');  grid(ax1,'on');

    % --- Panel 2: instantaneous omega ------------------------------------
    ax2 = subplot(2, ceil(n_panels/2), 2);  hold(ax2,'on');
    for p = 1:nP
        t  = pm(p).t;
        t_mid = t(1:numel(am(p).omega));   % omega is nFrames-1 long
        plot(ax2, t_mid, am(p).omega, '-', 'Color', colors(p,:), ...
            'LineWidth', 0.8, 'DisplayName', sprintf('P%d', p));
    end
    yline(ax2, 0, 'k--');
    xlabel(ax2,'t (s)');  ylabel(ax2,'\omega (rad/s)');
    title(ax2,'Instantaneous angular speed');  legend(ax2,'show');  grid(ax2,'on');

    % --- Panel 3: rolling mean omega -------------------------------------
    ax3 = subplot(2, ceil(n_panels/2), 3);  hold(ax3,'on');
    for p = 1:nP
        t     = pm(p).t;
        t_mid = t(1:numel(am(p).roll_mean_omega));
        plot(ax3, t_mid, am(p).roll_mean_omega, '-', 'Color', colors(p,:), ...
            'LineWidth', 1.5, 'DisplayName', sprintf('P%d', p));
    end
    yline(ax3, 0, 'k--');
    xlabel(ax3,'t (s)');  ylabel(ax3,'\omega_{roll} (rad/s)');
    title(ax3,'Rolling mean \omega');  legend(ax3,'show');  grid(ax3,'on');

    % --- Panel 4: FFT power spectrum of omega ----------------------------
    ax4 = subplot(2, ceil(n_panels/2), 4);  hold(ax4,'on');
    for p = 1:nP
        f = am(p).fft_freq;
        P = am(p).fft_power;
        if ~isscalar(f) && ~isscalar(P)
            plot(ax4, f, P, '-', 'Color', colors(p,:), 'LineWidth', 1, ...
                'DisplayName', sprintf('P%d (f_{dom}=%.3f Hz)', p, am(p).dominant_freq_hz));
        end
    end
    xlabel(ax4,'Frequency (Hz)');  ylabel(ax4,'Power (rad^2/s^2/Hz)');
    title(ax4,'FFT of \omega(t)');  legend(ax4,'show');  grid(ax4,'on');
    set(ax4,'XScale','log','YScale','log');

    if has_sync
        sm = mov.sync_metrics;
        t_ref = pm(1).t(1:numel(sm.R_instantaneous));

        % --- Panel 5: Kuramoto R(t) --------------------------------------
        ax5 = subplot(2, ceil(n_panels/2), 5);  hold(ax5,'on');
        plot(ax5, t_ref, sm.R_instantaneous, 'Color', [0.6 0.6 0.6], ...
            'LineWidth', 0.8, 'DisplayName', 'R(t)');
        plot(ax5, t_ref, sm.roll_R, 'k-', 'LineWidth', 2, ...
            'DisplayName', sprintf('Rolling mean  <R>=%.3f', sm.mean_R));
        yline(ax5, sm.mean_R, 'r--', ...
            sprintf('Mean R = %.3f', sm.mean_R), 'LabelHorizontalAlignment','left');
        ylim(ax5, [0 1.05]);
        xlabel(ax5,'t (s)');  ylabel(ax5,'R (Kuramoto)');
        title(ax5,'Synchronisation order parameter');
        legend(ax5,'show');  grid(ax5,'on');

        % --- Panel 6: pairwise phase difference --------------------------
        ax6 = subplot(2, ceil(n_panels/2), 6);  hold(ax6,'on');
        pair_colors = lines(numel(sm.pairwise_labels));
        for k = 1:numel(sm.pairwise_labels)
            dphi = sm.pairwise_phase_diff{k};
            t_k  = t_ref(1:numel(dphi));
            plot(ax6, t_k, dphi, '-', 'Color', pair_colors(k,:), 'LineWidth', 1, ...
                'DisplayName', sprintf('%s  PLI=%.3f', ...
                sm.pairwise_labels{k}, sm.phase_locking_index(k)));
        end
        xlabel(ax6,'t (s)');  ylabel(ax6,'\Delta\phi (rad)');
        title(ax6,'Pairwise phase difference (unwrapped)');
        legend(ax6,'show');  grid(ax6,'on');
    end

    sgtitle(sprintf('Angular analysis — Movie %d  |  %d particle(s)', ...
        mov.movie_idx, nP));
end
