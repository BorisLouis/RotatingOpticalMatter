function plot_full_movie_report(mov, meta)
% PLOT_FULL_MOVIE_REPORT  All extracted parameters for one movie in one figure set.
%
%   plot_full_movie_report(mov, meta)
%
%   INPUT
%   -----
%   mov   - one element of Results.movies{m}  (per-movie struct)
%   meta  - Results.meta  (for units / labels)
%
%   Produces TWO figures:
%     Figure A — Positional / spatial metrics
%     Figure B — Angular / synchronisation metrics  (if angular data present)

    pm   = mov.particle_metrics;
    nP   = mov.n_particles;
    midx = mov.movie_idx;
    pw   = mov.pairwise_metrics;
    has_pw   = nP > 1 && ~isempty(fieldnames(pw));
    has_ang  = isfield(mov, 'angular_metrics');
    has_sync = isfield(mov, 'sync_metrics') && nP > 1;

    colors     = lines(nP);
    pair_colors = lines(max(1, nP*(nP-1)/2));

    % =====================================================================
    %  FIGURE A — Positional metrics  (3 x 3 grid)
    % =====================================================================
    fA = figure('Name', sprintf('Movie %d — Positional metrics', midx), ...
        'NumberTitle','off','Position',[30 30 1500 900]);

    % ---- A1: 2D localisation scatter ------------------------------------
    ax = subplot(3,3,1);  hold(ax,'on');  axis(ax,'equal');
    for p = 1:nP
        x = pm{p}.x;   y = pm{p}.y;
        scatter(ax, x, y, 4, colors(p,:), 'filled', 'MarkerFaceAlpha', 0.25, ...
            'DisplayName', sprintf('P%d', p));
        % Mark mean position with a cross
        plot(ax, pm{p}.mean_x, pm{p}.mean_y, '+', 'Color', colors(p,:), ...
            'MarkerSize', 12, 'LineWidth', 2, 'HandleVisibility','off');
    end
    xlabel(ax,'col (px)');  ylabel(ax,'row (px)');
    title(ax,'2D localisation');  legend(ax,'show');  grid(ax,'on');
    ax.YDir = 'reverse';   % image convention: row increases downward

    % ---- A2: x(t) -------------------------------------------------------
    ax = subplot(3,3,2);  hold(ax,'on');
    for p = 1:nP
        plot(ax, pm{p}.t, pm{p}.x, '-', 'Color', colors(p,:), 'LineWidth',0.8, ...
            'DisplayName', sprintf('P%d', p));
    end
    xlabel(ax,'t (s)');  ylabel(ax,'col (px)');
    title(ax,'Column position vs time');  legend(ax,'show');  grid(ax,'on');

    % ---- A3: y(t) -------------------------------------------------------
    ax = subplot(3,3,3);  hold(ax,'on');
    for p = 1:nP
        plot(ax, pm{p}.t, pm{p}.y, '-', 'Color', colors(p,:), 'LineWidth',0.8, ...
            'DisplayName', sprintf('P%d', p));
    end
    xlabel(ax,'t (s)');  ylabel(ax,'row (px)');
    title(ax,'Row position vs time');  legend(ax,'show');  grid(ax,'on');

    % ---- A4: r(t) -------------------------------------------------------
    ax = subplot(3,3,4);  hold(ax,'on');
    for p = 1:nP
        plot(ax, pm{p}.t, pm{p}.r, '-', 'Color', colors(p,:), 'LineWidth',0.8, ...
            'DisplayName', sprintf('P%d  <r>=%.2f', p, nanmean(pm{p}.r)));
    end
    xlabel(ax,'t (s)');  ylabel(ax,'r (px)');
    title(ax,'Radial distance from CoM');  legend(ax,'show');  grid(ax,'on');

    % ---- A5: rolling std x and y ----------------------------------------
    ax = subplot(3,3,5);  hold(ax,'on');
    for p = 1:nP
        plot(ax, pm{p}.t, pm{p}.roll_std_x, '-',  'Color', colors(p,:), 'LineWidth',1, ...
            'DisplayName', sprintf('P%d \sigma_x', p));
        plot(ax, pm{p}.t, pm{p}.roll_std_y, '--', 'Color', colors(p,:), 'LineWidth',1, ...
            'HandleVisibility','off');
    end
    xlabel(ax,'t (s)');  ylabel(ax,'\sigma (px)');
    title(ax,sprintf('Rolling \\sigma_x (solid) / \\sigma_y (dashed)   win=%d fr', ...
        meta.rolling_window_frames));
    legend(ax,'show');  grid(ax,'on');

    % ---- A6: total speed ------------------------------------------------
    ax = subplot(3,3,6);  hold(ax,'on');
    for p = 1:nP
        t_mid = pm{p}.t(1:numel(pm{p}.speed_total));
        plot(ax, t_mid, pm{p}.speed_total, '-', 'Color', colors(p,:), 'LineWidth',0.8, ...
            'DisplayName', sprintf('P%d  <v>=%.3f', p, pm{p}.mean_speed_total));
    end
    xlabel(ax,'t (s)');  ylabel(ax,'Speed (px/s)');
    title(ax,'Total instantaneous speed');  legend(ax,'show');  grid(ax,'on');

    % ---- A7: interparticle distance  (multi only) -----------------------
    ax = subplot(3,3,7);  hold(ax,'on');
    if has_pw
        for k = 1:numel(pw.pairs)
            pr  = pw.pairs{k};
            t_k = pm{1}.t(1:numel(pr.dist));
            plot(ax, t_k, pr.dist,      '-',  'Color', pair_colors(k,:), 'LineWidth',0.8, ...
                'DisplayName', sprintf('%s  <d>=%.2f', pr.label, pr.mean_dist));
            plot(ax, t_k, pr.roll_dist, '--', 'Color', pair_colors(k,:), 'LineWidth',1.5, ...
                'HandleVisibility','off');
        end
    else
        text(ax,0.5,0.5,'Single particle — N/A','Units','normalized', ...
            'HorizontalAlignment','center','Color',[0.5 0.5 0.5]);
    end
    xlabel(ax,'t (s)');  ylabel(ax,'Distance (px)');
    title(ax,'Interparticle distance (solid) + rolling mean (dashed)');
    legend(ax,'show');  grid(ax,'on');

    % ---- A8: motion correlation rolling  (multi only) -------------------
    ax = subplot(3,3,8);  hold(ax,'on');
    if has_pw
        for k = 1:numel(pw.pairs)
            pr  = pw.pairs{k};
            t_k = pm{1}.t(1:numel(pr.roll_corr_x));
            plot(ax, t_k, pr.roll_corr_x, '-',  'Color', pair_colors(k,:), 'LineWidth',1, ...
                'DisplayName', sprintf('%s x', pr.label));
            plot(ax, t_k, pr.roll_corr_y, '--', 'Color', pair_colors(k,:), 'LineWidth',1, ...
                'HandleVisibility','off');
        end
        yline(ax, 0, 'k--');
        ylim(ax, [-1.1 1.1]);
    else
        text(ax,0.5,0.5,'Single particle — N/A','Units','normalized', ...
            'HorizontalAlignment','center','Color',[0.5 0.5 0.5]);
    end
    xlabel(ax,'t (s)');  ylabel(ax,'Pearson r');
    title(ax,'Rolling motion correlation  x (solid) / y (dashed)');
    legend(ax,'show');  grid(ax,'on');

    % ---- A9: CoM drift  (multi only) ------------------------------------
    ax = subplot(3,3,9);  hold(ax,'on');
    if has_pw
        t_com = pm{1}.t(1:numel(pw.com_x));
        plot(ax, t_com, pw.com_x, 'b-', 'LineWidth',1, 'DisplayName','CoM col');
        plot(ax, t_com, pw.com_y, 'r-', 'LineWidth',1, 'DisplayName','CoM row');
    else
        plot(ax, pm{1}.t, pm{1}.x, 'b-', 'LineWidth',1, 'DisplayName','col');
        plot(ax, pm{1}.t, pm{1}.y, 'r-', 'LineWidth',1, 'DisplayName','row');
    end
    xlabel(ax,'t (s)');  ylabel(ax,'Position (px)');
    title(ax, 'Centre-of-mass drift');  legend(ax,'show');  grid(ax,'on');

    sgtitle(fA, sprintf('Movie %d  |  %d particle(s)  |  \\lambda=%d nm  |  exp=%.4f s', ...
        midx, nP, meta.wavelength_nm, meta.exposure_time_s), 'FontSize',13);

    % =====================================================================
    %  FIGURE B — Angular / synchronisation  (only if computed)
    % =====================================================================
    if ~has_ang
        return;
    end

    am = mov.angular_metrics;

    fB = figure('Name', sprintf('Movie %d — Angular metrics', midx), ...
        'NumberTitle','off','Position',[60 60 1500 750]);

    % ---- B1: unwrapped theta --------------------------------------------
    ax = subplot(2,3,1);  hold(ax,'on');
    for p = 1:nP
        plot(ax, pm{p}.t, am{p}.theta_unwrap, '-', 'Color', colors(p,:), 'LineWidth',1, ...
            'DisplayName', sprintf('P%d  f=%.4f Hz', p, am{p}.rot_freq_mean));
    end
    xlabel(ax,'t (s)');  ylabel(ax,'\theta_{unwrap} (rad)');
    title(ax,'Cumulative angle  (slope = rotation freq)');
    legend(ax,'show');  grid(ax,'on');

    % ---- B2: instantaneous omega ----------------------------------------
    ax = subplot(2,3,2);  hold(ax,'on');
    for p = 1:nP
        t_mid = pm{p}.t(1:numel(am{p}.omega));
        plot(ax, t_mid, am{p}.omega, '-', 'Color', colors(p,:)*0.7, 'LineWidth',0.6);
        plot(ax, t_mid, am{p}.roll_mean_omega, '-', 'Color', colors(p,:), 'LineWidth',2, ...
            'DisplayName', sprintf('P%d  <\\omega>=%.3f rad/s', p, am{p}.mean_omega));
    end
    yline(ax,0,'k--');
    xlabel(ax,'t (s)');  ylabel(ax,'\omega (rad/s)');
    title(ax,'Angular speed: raw (light) + rolling mean (bold)');
    legend(ax,'show');  grid(ax,'on');

    % ---- B3: FFT power spectrum -----------------------------------------
    ax = subplot(2,3,3);  hold(ax,'on');
    for p = 1:nP
        f = am{p}.fft_freq;
        P = am{p}.fft_power;
        if numel(f) > 1
            plot(ax, f(2:end), P(2:end), '-', 'Color', colors(p,:), 'LineWidth',1, ...
                'DisplayName', sprintf('P%d  f_{dom}=%.3f Hz', p, am{p}.dominant_freq_hz));
        end
    end
    xlabel(ax,'Frequency (Hz)');  ylabel(ax,'Power (rad^2 s^{-2} Hz^{-1})');
    title(ax,'Power spectrum of \omega  (DC removed)');
    legend(ax,'show');  grid(ax,'on');
    set(ax,'XScale','log','YScale','log');

    % ---- B4: rolling std of omega  (fluctuation amplitude) --------------
    ax = subplot(2,3,4);  hold(ax,'on');
    for p = 1:nP
        t_mid = pm{p}.t(1:numel(am{p}.roll_std_omega));
        plot(ax, t_mid, am{p}.roll_std_omega, '-', 'Color', colors(p,:), 'LineWidth',1, ...
            'DisplayName', sprintf('P%d', p));
    end
    xlabel(ax,'t (s)');  ylabel(ax,'\sigma_\omega (rad/s)');
    title(ax,'Rolling \sigma_\omega  (angular speed fluctuations)');
    legend(ax,'show');  grid(ax,'on');

    % ---- B5: Kuramoto R(t)  (multi only) --------------------------------
    ax = subplot(2,3,5);  hold(ax,'on');
    if has_sync
        sm    = mov.sync_metrics;
        t_ref = pm{1}.t(1:numel(sm.R_instantaneous));
        plot(ax, t_ref, sm.R_instantaneous, 'Color',[0.7 0.7 0.7], 'LineWidth',0.8);
        plot(ax, t_ref, sm.roll_R, 'k-', 'LineWidth',2, ...
            'DisplayName', sprintf('Rolling R  <R>=%.3f', sm.mean_R));
        yline(ax, sm.mean_R, 'r--', sprintf('Mean R=%.3f', sm.mean_R), ...
            'LabelHorizontalAlignment','left');
        ylim(ax,[0 1.05]);
    else
        text(ax,0.5,0.5,'Single particle — N/A','Units','normalized', ...
            'HorizontalAlignment','center','Color',[0.5 0.5 0.5]);
    end
    xlabel(ax,'t (s)');  ylabel(ax,'R (Kuramoto)');
    title(ax,'Synchronisation order parameter');  legend(ax,'show');  grid(ax,'on');

    % ---- B6: pairwise phase difference  (multi only) --------------------
    ax = subplot(2,3,6);  hold(ax,'on');
    if has_sync
        sm = mov.sync_metrics;
        t_ref = pm{1}.t(1:numel(sm.R_instantaneous));
        for k = 1:numel(sm.pairwise_labels)
            dphi = sm.pairwise_phase_diff{k};
            t_k  = t_ref(1:numel(dphi));
            plot(ax, t_k, dphi, '-', 'Color', pair_colors(k,:), 'LineWidth',1, ...
                'DisplayName', sprintf('%s  PLI=%.3f', ...
                sm.pairwise_labels{k}, sm.phase_locking_index(k)));
        end
    else
        text(ax,0.5,0.5,'Single particle — N/A','Units','normalized', ...
            'HorizontalAlignment','center','Color',[0.5 0.5 0.5]);
    end
    xlabel(ax,'t (s)');  ylabel(ax,'\Delta\phi (rad)');
    title(ax,'Pairwise phase difference (unwrapped)');
    legend(ax,'show');  grid(ax,'on');

    sgtitle(fB, sprintf('Movie %d  |  %d particle(s)  |  Angular analysis', ...
        midx, nP), 'FontSize',13);
end
