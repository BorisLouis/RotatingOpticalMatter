function plot_movie_results(Results, movie_number)
% PLOT_MOVIE_RESULTS  Quick diagnostic plots for a single movie.
%
%   plot_movie_results(Results)
%       Prompts for movie number interactively.
%
%   plot_movie_results(Results, movie_number)
%       Directly plots movie_number (1-based index into Results.movies).
%
%   Produces a figure with:
%     Panel 1: x(t) and y(t) for all particles
%     Panel 2: radial distance r(t) for all particles
%     Panel 3: rolling std of x (all particles)
%     Panel 4: rolling std of r (all particles)
%     Panel 5: (multi-particle) interparticle distance vs time
%     Panel 6: (multi-particle) rolling correlation in x
%
%   USAGE EXAMPLE
%   -------------
%   % After running main_analyze_optical_trapping:
%   plot_movie_results(Results, 1)

    n_movies = numel(Results.movies);

    if nargin < 2
        movie_number = input(sprintf('Enter movie number to plot [1-%d]: ', n_movies));
    end

    if movie_number < 1 || movie_number > n_movies
        error('OT:BadMovieIndex', 'movie_number must be between 1 and %d.', n_movies);
    end

    mov = Results.movies{movie_number};
    pm  = mov.particle_metrics;
    nP  = mov.n_particles;
    pw  = mov.pairwise_metrics;

    colors = lines(nP);   % one colour per particle

    n_panels = 4 + 2*(nP > 1);
    fig = figure('Name', sprintf('Movie %d  (%d particle(s))', ...
        mov.movie_idx, nP), 'NumberTitle', 'off', ...
        'Position', [50 50 1200 700]);

    % --- Panel 1: x(t) and y(t) -----------------------------------------
    ax1 = subplot(2, n_panels/2, 1);   hold(ax1, 'on');
    for p = 1:nP
        t = pm{p}.t;
        plot(ax1, t, pm{p}.x, '-', 'Color', colors(p,:), 'LineWidth', 1, ...
            'DisplayName', sprintf('P%d x', p));
        plot(ax1, t, pm{p}.y, '--', 'Color', colors(p,:), 'LineWidth', 1, ...
            'DisplayName', sprintf('P%d y', p));
    end
    xlabel(ax1, 't (s)');  ylabel(ax1, 'Position (px)');
    title(ax1, 'x(t), y(t)');  legend(ax1, 'show');  grid(ax1, 'on');

    % --- Panel 2: r(t) --------------------------------------------------
    ax2 = subplot(2, n_panels/2, 2);   hold(ax2, 'on');
    for p = 1:nP
        t = pm{p}.t;
        plot(ax2, t, pm{p}.r, '-', 'Color', colors(p,:), 'LineWidth', 1, ...
            'DisplayName', sprintf('P%d', p));
    end
    xlabel(ax2, 't (s)');  ylabel(ax2, 'r (px)');
    title(ax2, 'Radial distance from CoM');  legend(ax2, 'show');  grid(ax2, 'on');

    % --- Panel 3: rolling std of x --------------------------------------
    ax3 = subplot(2, n_panels/2, 3);   hold(ax3, 'on');
    for p = 1:nP
        t = pm{p}.t;
        plot(ax3, t, pm{p}.roll_std_x, '-', 'Color', colors(p,:), 'LineWidth', 1, ...
            'DisplayName', sprintf('P%d', p));
    end
    xlabel(ax3, 't (s)');  ylabel(ax3, '\sigma_x (px)');
    title(ax3, 'Rolling \sigma_x');  legend(ax3, 'show');  grid(ax3, 'on');

    % --- Panel 4: rolling std of r --------------------------------------
    ax4 = subplot(2, n_panels/2, 4);   hold(ax4, 'on');
    for p = 1:nP
        t = pm{p}.t;
        plot(ax4, t, pm{p}.roll_std_r, '-', 'Color', colors(p,:), 'LineWidth', 1, ...
            'DisplayName', sprintf('P%d', p));
    end
    xlabel(ax4, 't (s)');  ylabel(ax4, '\sigma_r (px)');
    title(ax4, 'Rolling \sigma_r');  legend(ax4, 'show');  grid(ax4, 'on');

    % --- Panels 5 & 6: multi-particle only ------------------------------
    if nP > 1 && ~isempty(fieldnames(pw))
        n_pairs = numel(pw.pairs);
        pair_colors = lines(n_pairs);
        t_ref = pm{1}.t(1:numel(pw.pairs{1}.dist));   % common time axis

        % Panel 5: interparticle distance
        ax5 = subplot(2, n_panels/2, 5);   hold(ax5, 'on');
        for k = 1:n_pairs
            pr = pw.pairs{k};
            t_k = t_ref(1:numel(pr.dist));
            plot(ax5, t_k, pr.dist, '-', 'Color', pair_colors(k,:), ...
                'LineWidth', 1, 'DisplayName', pr.label);
            plot(ax5, t_k, pr.roll_dist, '--', 'Color', pair_colors(k,:), ...
                'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
        xlabel(ax5, 't (s)');  ylabel(ax5, 'Distance (px)');
        title(ax5, 'Interparticle distance (solid) + rolling mean (dashed)');
        legend(ax5, 'show');  grid(ax5, 'on');

        % Panel 6: rolling correlation in x
        ax6 = subplot(2, n_panels/2, 6);   hold(ax6, 'on');
        for k = 1:n_pairs
            pr = pw.pairs{k};
            t_k = t_ref(1:numel(pr.roll_corr_x));
            plot(ax6, t_k, pr.roll_corr_x, '-', 'Color', pair_colors(k,:), ...
                'LineWidth', 1, 'DisplayName', [pr.label ' x']);
            plot(ax6, t_k, pr.roll_corr_y, '--', 'Color', pair_colors(k,:), ...
                'LineWidth', 1, 'DisplayName', [pr.label ' y']);
        end
        yline(ax6, 0, 'k--');
        xlabel(ax6, 't (s)');  ylabel(ax6, 'Pearson r');
        title(ax6, 'Rolling motion correlation (x solid, y dashed)');
        legend(ax6, 'show');  grid(ax6, 'on');
        ylim(ax6, [-1.1 1.1]);
    end

    sgtitle(fig, sprintf('Movie %d  |  %d particle(s)  |  %s', ...
        mov.movie_idx, nP, Results.meta.analysis_date));
end
