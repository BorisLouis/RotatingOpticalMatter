function pairwise_metrics = compute_pairwise_metrics(particle_metrics, win)
% COMPUTE_PAIRWISE_METRICS  Pairwise distances and motion correlations.
%
%   pairwise_metrics = compute_pairwise_metrics(particle_metrics, win)
%
%   Handles movies with 2, 3, or 4 particles.  For each ordered pair
%   (i < j) it computes interparticle distance and motion correlation.
%
%   INPUT
%   -----
%   particle_metrics  - {nP x 1} cell from compute_particle_metrics
%   win               - rolling window size in frames
%
%   OUTPUT
%   ------
%   pairwise_metrics  - struct with fields:
%
%     .pairs           - {nPairs x 1} cell of pair structs, each containing:
%         .label           string, e.g. 'P1-P2'
%         .i, .j           1-based particle indices
%
%         .dist            interparticle distance [nFrames_min x 1]
%         .mean_dist       scalar, time-averaged distance
%         .std_dist        scalar, std of distance
%
%         .roll_dist       rolling-window mean distance [nFrames_min x 1]
%
%         .corr_x          Pearson corr of dx displacements (full trace)
%         .corr_y          Pearson corr of dy displacements (full trace)
%         .corr_r          Pearson corr of dr (radial displacements)
%
%         .roll_corr_x     rolling-window Pearson corr of dx [nFrames-1 x 1]
%         .roll_corr_y     rolling-window Pearson corr of dy
%         .roll_corr_r     rolling-window Pearson corr of dr
%
%         .mean_roll_corr_x  scalar, average of roll_corr_x
%         .mean_roll_corr_y  scalar
%         .mean_roll_corr_r  scalar
%
%     .com_x, .com_y     instantaneous CoM across all particles [nMin x 1]
%     .mean_com_x/y      scalar, time-averaged CoM
%
%   PAIR NAMING
%   -----------
%   Pairs are named 'P{i}-P{j}' with i < j, matching particle order in
%   particle_metrics.  This naming is stable across all movies in the
%   dataset as long as the order of traces is consistent (which it is,
%   since group_traces_by_movie preserves appearance order).
%
%   TRACE LENGTH MISMATCH
%   ----------------------
%   All pairwise calculations are performed on the SHORTEST common trace
%   length across the two particles in the pair (trimming at the end).
%   The .dist and rolling windows respect this trimmed length.
%
%   CORRELATION DEFINITION
%   -----------------------
%   Pearson correlation of INSTANTANEOUS DISPLACEMENTS (not raw positions).
%   Displacement at frame k = position(k+1) - position(k).
%   This removes the effect of mean position offsets and tests whether
%   the two particles move in the same direction simultaneously.
%
%   ROLLING CORRELATION
%   --------------------
%   Within each window, Pearson correlation is computed between the two
%   displacement vectors.  Windows with fewer than 3 valid points return NaN.
%   The output vector has the same length as the displacement vector,
%   using a centred window (consistent with compute_particle_metrics).

    n_particles = numel(particle_metrics);

    % Build list of pairs
    pairs = {};
    pair_count = 0;
    for i = 1:n_particles
        for j = (i+1):n_particles
            pair_count = pair_count + 1;
            pairs{pair_count} = struct('i', i, 'j', j, ...
                'label', sprintf('P%d-P%d', i, j)); %#ok<AGROW>
        end
    end

    % =====================================================================
    %  Recompute instantaneous CoM using all particles (common length)
    % =====================================================================
    nF_all  = cellfun(@(pm) numel(pm.x), particle_metrics);
    nF_min  = min(nF_all);

    com_x = zeros(nF_min, 1);
    com_y = zeros(nF_min, 1);
    for p = 1:n_particles
        com_x = com_x + particle_metrics{p}.x(1:nF_min);
        com_y = com_y + particle_metrics{p}.y(1:nF_min);
    end
    com_x = com_x / n_particles;
    com_y = com_y / n_particles;

    % Overwrite the radial coordinate in particle_metrics to use true CoM
    % (This is the only place where particle_metrics is modified after
    %  compute_particle_metrics.  All downstream code uses this updated r.)
    for p = 1:n_particles
        x_p = particle_metrics{p}.x(1:nF_min);
        y_p = particle_metrics{p}.y(1:nF_min);
        dx  = x_p - com_x;
        dy  = y_p - com_y;
        r_new     = sqrt(dx.^2 + dy.^2);
        theta_new = atan2(dy, dx);

        % Extend with NaN if this particle's trace is longer than nF_min
        nF_p = numel(particle_metrics{p}.x);
        if nF_p > nF_min
            r_new     = [r_new;     NaN(nF_p - nF_min, 1)];
            theta_new = [theta_new; NaN(nF_p - nF_min, 1)];
        end
        particle_metrics{p}.r     = r_new;
        particle_metrics{p}.theta = theta_new;
        particle_metrics{p}.std_r = nanstd(r_new);
    end

    % =====================================================================
    %  Per-pair computation
    % =====================================================================
    for k = 1:pair_count
        i = pairs{k}.i;
        j = pairs{k}.j;

        % Trim to common length for this pair
        nF_i = numel(particle_metrics{i}.x);
        nF_j = numel(particle_metrics{j}.x);
        nF   = min([nF_i, nF_j, nF_min]);

        xi = particle_metrics{i}.x(1:nF);
        yi = particle_metrics{i}.y(1:nF);
        ri = particle_metrics{i}.r(1:nF);

        xj = particle_metrics{j}.x(1:nF);
        yj = particle_metrics{j}.y(1:nF);
        rj = particle_metrics{j}.r(1:nF);

        % --- Interparticle distance --------------------------------------
        dist = sqrt((xi - xj).^2 + (yi - yj).^2);
        mean_dist = nanmean(dist);
        std_dist  = nanstd(dist);

        % Rolling mean of distance
        roll_dist = rolling_nanmean(dist, win);

        % --- Displacements (for correlation) -----------------------------
        dxi = diff(xi);   dxi(isnan(dxi)) = 0;
        dyi = diff(yi);   dyi(isnan(dyi)) = 0;
        dri = diff(ri);   dri(isnan(dri)) = 0;

        dxj = diff(xj);   dxj(isnan(dxj)) = 0;
        dyj = diff(yj);   dyj(isnan(dyj)) = 0;
        drj = diff(rj);   drj(isnan(drj)) = 0;

        % Full-trace Pearson correlation of displacements
        corr_x = safe_corr(dxi, dxj);
        corr_y = safe_corr(dyi, dyj);
        corr_r = safe_corr(dri, drj);

        % Rolling Pearson correlation
        roll_corr_x = rolling_corr(dxi, dxj, win);
        roll_corr_y = rolling_corr(dyi, dyj, win);
        roll_corr_r = rolling_corr(dri, drj, win);

        mean_roll_corr_x = nanmean(roll_corr_x);
        mean_roll_corr_y = nanmean(roll_corr_y);
        mean_roll_corr_r = nanmean(roll_corr_r);

        % --- Store -------------------------------------------------------
        pairs{k}.dist             = dist;
        pairs{k}.mean_dist        = mean_dist;
        pairs{k}.std_dist         = std_dist;
        pairs{k}.roll_dist        = roll_dist;
        pairs{k}.corr_x           = corr_x;
        pairs{k}.corr_y           = corr_y;
        pairs{k}.corr_r           = corr_r;
        pairs{k}.roll_corr_x      = roll_corr_x;
        pairs{k}.roll_corr_y      = roll_corr_y;
        pairs{k}.roll_corr_r      = roll_corr_r;
        pairs{k}.mean_roll_corr_x = mean_roll_corr_x;
        pairs{k}.mean_roll_corr_y = mean_roll_corr_y;
        pairs{k}.mean_roll_corr_r = mean_roll_corr_r;
    end

    pairwise_metrics = struct( ...
        'pairs',      {pairs},   ...
        'com_x',      com_x,     ...
        'com_y',      com_y,     ...
        'mean_com_x', mean(com_x), ...
        'mean_com_y', mean(com_y)  );
end


% =========================================================================
%  LOCAL UTILITIES
% =========================================================================

function r = safe_corr(a, b)
% Pearson correlation with NaN-safety; returns NaN if insufficient data.
    valid = ~isnan(a) & ~isnan(b);
    if sum(valid) < 3
        r = NaN;
        return;
    end
    C = corrcoef(a(valid), b(valid));
    r = C(1, 2);
end


function out = rolling_corr(a, b, win)
% Rolling Pearson correlation of two equal-length column vectors.
% Output length equals input length; edges return NaN.
    n    = numel(a);
    out  = NaN(n, 1);
    half = floor(win / 2);

    for i = 1:n
        i0 = max(1, i - half);
        i1 = min(n, i + (win - half - 1));
        sa = a(i0:i1);
        sb = b(i0:i1);
        valid = ~isnan(sa) & ~isnan(sb);
        if sum(valid) >= 3
            C = corrcoef(sa(valid), sb(valid));
            out(i) = C(1, 2);
        end
    end
end


function out = rolling_nanmean(x, win)
% Rolling mean, centred window, NaN-aware.
    n    = numel(x);
    out  = NaN(n, 1);
    half = floor(win / 2);

    for i = 1:n
        i0   = max(1, i - half);
        i1   = min(n, i + (win - half - 1));
        seg  = x(i0:i1);
        if sum(~isnan(seg)) >= 1
            out(i) = nanmean(seg);
        end
    end
end
