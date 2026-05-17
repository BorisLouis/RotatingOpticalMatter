function pairwise_metrics = compute_pairwise_metrics(particle_metrics, win)
% COMPUTE_PAIRWISE_METRICS  Pairwise distances, angles, and motion correlations.
%
%   pairwise_metrics = compute_pairwise_metrics(particle_metrics, win)
%
%   Handles movies with 2, 3, or 4 particles.  For each ordered pair
%   (i < j) it computes interparticle distance, bond angle, and motion
%   correlations (x, y, r, theta).
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
%
%         .label           string, e.g. 'P1-P2'
%         .i, .j           1-based particle indices
%
%         --- Interparticle distance ---
%         .dist            instantaneous distance [nF x 1]
%         .mean_dist       scalar, time-averaged distance
%         .std_dist        scalar, std of distance
%         .roll_dist       rolling-window mean distance [nF x 1]
%
%         --- Bond angle (lab-frame orientation of the i->j vector) ---
%         .phi             atan2(yj-yi, xj-xi), wrapped to (-pi, pi] [nF x 1]
%         .phi_unwrapped   continuously unwrapped bond angle [nF x 1]
%         .mean_phi        circular mean of phi [scalar, radians]
%         .std_phi_circ    circular std of phi [scalar, radians]
%                            -> 0 = perfectly stable orientation
%                            -> sqrt(2) ~ 1.41 = maximally dispersed
%         .std_phi_unwrap  std of phi_unwrapped [scalar, radians]
%                            -> measures total angular wandering over time
%                            -> more sensitive to slow rotational drift
%         .roll_phi        rolling circular mean of phi [nF x 1]
%         .roll_std_phi    rolling circular std of phi [nF x 1]
%
%         --- Displacement correlations (x, y, r, theta) ---
%         .corr_x          Pearson corr of dx displacements (full trace)
%         .corr_y          Pearson corr of dy displacements
%         .corr_r          Pearson corr of dr (radial displacements)
%         .corr_theta      Pearson corr of dtheta (angular displacements,
%                            after unwrapping to remove 2pi jumps)
%
%         .roll_corr_x     rolling Pearson corr of dx [nF-1 x 1]
%         .roll_corr_y     rolling Pearson corr of dy
%         .roll_corr_r     rolling Pearson corr of dr
%         .roll_corr_theta rolling Pearson corr of dtheta
%
%         .mean_roll_corr_x     scalar averages of the rolling correlations
%         .mean_roll_corr_y
%         .mean_roll_corr_r
%         .mean_roll_corr_theta
%
%     .com_x, .com_y     instantaneous CoM across all particles [nMin x 1]
%     .mean_com_x/y      scalar, time-averaged CoM
%
%   BOND ANGLE DEFINITION
%   ----------------------
%   phi_ij(t) = atan2( yj(t) - yi(t),  xj(t) - xi(t) )
%
%   This is the orientation of the bond vector pointing from particle i
%   to particle j in the lab frame.  It characterises the absolute
%   orientation of each pair and is the primary descriptor for assembly
%   geometry classification (e.g. square vs rhombus for 4 particles).
%
%   Two complementary statistics are provided:
%     std_phi_circ   : circular std -> stability of the preferred orientation
%     std_phi_unwrap : std of the unwrapped trace -> sensitivity to slow
%                      rotational drift of the whole assembly
%
%   THETA CORRELATION DEFINITION
%   -----------------------------
%   theta_i is the polar angle of particle i relative to the CoM, updated
%   in the CoM block below.  Before differencing, each theta trace is
%   unwrapped (MATLAB unwrap) to remove discontinuities at ±pi.  The
%   displacement dtheta_i = diff(unwrap(theta_i)) is then correlated with
%   dtheta_j.  A high positive correlation means both particles rotate
%   around the CoM in the same direction simultaneously (co-rotation).
%
%   PAIR NAMING / TRACE LENGTH / ROLLING WINDOWS
%   ---------------------------------------------
%   See original header; these rules are unchanged.

    n_particles = numel(particle_metrics);

    % Build list of pairs -------------------------------------------------
    pairs      = {};
    pair_count = 0;
    for i = 1:n_particles
        for j = (i+1):n_particles
            pair_count = pair_count + 1;
            pairs{pair_count} = struct('i', i, 'j', j, ...
                'label', sprintf('P%d-P%d', i, j)); %#ok<AGROW>
        end
    end

    % =====================================================================
    %  Recompute instantaneous CoM and polar coordinates
    % =====================================================================
    nF_all = cellfun(@(pm) numel(pm.x), particle_metrics);
    nF_min = min(nF_all);

    com_x = zeros(nF_min, 1);
    com_y = zeros(nF_min, 1);
    for p = 1:n_particles
        com_x = com_x + particle_metrics{p}.x(1:nF_min);
        com_y = com_y + particle_metrics{p}.y(1:nF_min);
    end
    com_x = com_x / n_particles;
    com_y = com_y / n_particles;

    % Update r and theta for every particle using the true CoM.
    % theta is stored as the raw atan2 output in (-pi, pi]; unwrapping is
    % done locally below when needed for displacement computation.
    for p = 1:n_particles
        x_p = particle_metrics{p}.x(1:nF_min);
        y_p = particle_metrics{p}.y(1:nF_min);
        dx  = x_p - com_x;
        dy  = y_p - com_y;

        r_new     = sqrt(dx.^2 + dy.^2);
        theta_new = atan2(dy, dx);          % wrapped, (-pi, pi]

        nF_p = numel(particle_metrics{p}.x);
        if nF_p > nF_min
            r_new     = [r_new;     NaN(nF_p - nF_min, 1)];
            theta_new = [theta_new; NaN(nF_p - nF_min, 1)];
        end

        particle_metrics{p}.r     = r_new;
        particle_metrics{p}.theta = theta_new;          % wrapped
        particle_metrics{p}.std_r = nanstd(r_new);
    end

    % =====================================================================
    %  Per-pair computation
    % =====================================================================
    for k = 1:pair_count
        i = pairs{k}.i;
        j = pairs{k}.j;

        % Common trace length for this pair
        nF_i = numel(particle_metrics{i}.x);
        nF_j = numel(particle_metrics{j}.x);
        nF   = min([nF_i, nF_j, nF_min]);

        xi = particle_metrics{i}.x(1:nF);
        yi = particle_metrics{i}.y(1:nF);
        ri = particle_metrics{i}.r(1:nF);

        xj = particle_metrics{j}.x(1:nF);
        yj = particle_metrics{j}.y(1:nF);
        rj = particle_metrics{j}.r(1:nF);

        % Unwrap theta before use (particle-level, not pair-level)
        theta_i_raw = particle_metrics{i}.theta(1:nF);
        theta_j_raw = particle_metrics{j}.theta(1:nF);
        theta_i_uw  = safe_unwrap(theta_i_raw);
        theta_j_uw  = safe_unwrap(theta_j_raw);

        % -----------------------------------------------------------------
        %  Interparticle distance
        % -----------------------------------------------------------------
        dist      = sqrt((xi - xj).^2 + (yi - yj).^2);
        mean_dist = nanmean(dist);
        std_dist  = nanstd(dist);
        roll_dist = rolling_nanmean(dist, win);

        % -----------------------------------------------------------------
        %  Bond angle: orientation of the i->j vector in the lab frame
        % -----------------------------------------------------------------
        phi_raw   = atan2(yj - yi, xj - xi);       % wrapped, (-pi, pi]
        phi_uw    = safe_unwrap(phi_raw);            % continuously unwrapped

        % Circular mean: angle of the mean unit vector
        mean_phi  = angle(nanmean(exp(1i * phi_raw)));

        % Circular std: sqrt(-2*ln(R)) where R is the mean resultant length
        %   -> 0 for perfectly stable, sqrt(2) for fully dispersed
        R_mean         = abs(nanmean(exp(1i * phi_raw)));
        std_phi_circ   = sqrt(-2 * log(max(R_mean, eps)));   % guard log(0)

        % Std of unwrapped trace: captures slow rotational drift
        std_phi_unwrap = nanstd(phi_uw);

        % Rolling circular mean and std
        [roll_phi, roll_std_phi] = rolling_circ_stats(phi_raw, win);

        % -----------------------------------------------------------------
        %  Displacements (for translational and angular correlations)
        % -----------------------------------------------------------------
        dxi = diff(xi);      dxi(isnan(dxi)) = 0;
        dyi = diff(yi);      dyi(isnan(dyi)) = 0;
        dri = diff(ri);      dri(isnan(dri)) = 0;

        dxj = diff(xj);      dxj(isnan(dxj)) = 0;
        dyj = diff(yj);      dyj(isnan(dyj)) = 0;
        drj = diff(rj);      drj(isnan(drj)) = 0;

        % Angular displacements from unwrapped theta (avoids 2pi artefacts)
        dthi = diff(theta_i_uw);   dthi(isnan(dthi)) = 0;
        dthj = diff(theta_j_uw);   dthj(isnan(dthj)) = 0;

        % Full-trace Pearson correlations
        corr_x     = safe_corr(dxi,  dxj);
        corr_y     = safe_corr(dyi,  dyj);
        corr_r     = safe_corr(dri,  drj);
        corr_theta = safe_corr(dthi, dthj);

        % Rolling Pearson correlations
        roll_corr_x     = rolling_corr(dxi,  dxj,  win);
        roll_corr_y     = rolling_corr(dyi,  dyj,  win);
        roll_corr_r     = rolling_corr(dri,  drj,  win);
        roll_corr_theta = rolling_corr(dthi, dthj, win);

        mean_roll_corr_x     = nanmean(roll_corr_x);
        mean_roll_corr_y     = nanmean(roll_corr_y);
        mean_roll_corr_r     = nanmean(roll_corr_r);
        mean_roll_corr_theta = nanmean(roll_corr_theta);

        % -----------------------------------------------------------------
        %  Store all fields
        % -----------------------------------------------------------------
        % Distance
        pairs{k}.dist      = dist;
        pairs{k}.mean_dist = mean_dist;
        pairs{k}.std_dist  = std_dist;
        pairs{k}.roll_dist = roll_dist;

        % Bond angle
        pairs{k}.phi            = phi_raw;
        pairs{k}.phi_unwrapped  = phi_uw;
        pairs{k}.mean_phi       = mean_phi;
        pairs{k}.std_phi_circ   = std_phi_circ;
        pairs{k}.std_phi_unwrap = std_phi_unwrap;
        pairs{k}.roll_phi       = roll_phi;
        pairs{k}.roll_std_phi   = roll_std_phi;

        % Translational correlations
        pairs{k}.corr_x           = corr_x;
        pairs{k}.corr_y           = corr_y;
        pairs{k}.corr_r           = corr_r;
        pairs{k}.roll_corr_x      = roll_corr_x;
        pairs{k}.roll_corr_y      = roll_corr_y;
        pairs{k}.roll_corr_r      = roll_corr_r;
        pairs{k}.mean_roll_corr_x = mean_roll_corr_x;
        pairs{k}.mean_roll_corr_y = mean_roll_corr_y;
        pairs{k}.mean_roll_corr_r = mean_roll_corr_r;

        % Angular correlation (theta around CoM)
        pairs{k}.corr_theta           = corr_theta;
        pairs{k}.roll_corr_theta      = roll_corr_theta;
        pairs{k}.mean_roll_corr_theta = mean_roll_corr_theta;
    end

    pairwise_metrics = struct( ...
        'pairs',      {pairs},      ...
        'com_x',      com_x,        ...
        'com_y',      com_y,        ...
        'mean_com_x', mean(com_x),  ...
        'mean_com_y', mean(com_y)   );
end


% =========================================================================
%  LOCAL UTILITIES
% =========================================================================

function th_uw = safe_unwrap(th)
% SAFE_UNWRAP  Unwrap a phase vector, skipping NaN entries.
%
%   NaN values are left in place; unwrapping is applied only to the valid
%   segments so that isolated NaN gaps do not corrupt the cumulative phase.
    th_uw = th;
    valid = ~isnan(th);
    if sum(valid) < 2
        return
    end
    % Unwrap the valid samples, then put them back in their original positions
    th_uw(valid) = unwrap(th(valid));
end


function r = safe_corr(a, b)
% SAFE_CORR  Pearson correlation with NaN-safety.
    valid = ~isnan(a) & ~isnan(b);
    if sum(valid) < 3
        r = NaN;
        return
    end
    C = corrcoef(a(valid), b(valid));
    r = C(1, 2);
end


function out = rolling_nanmean(x, win)
% ROLLING_NANMEAN  Centred rolling mean, NaN-aware.
    n    = numel(x);
    out  = NaN(n, 1);
    half = floor(win / 2);
    for i = 1:n
        i0  = max(1, i - half);
        i1  = min(n, i + (win - half - 1));
        seg = x(i0:i1);
        if any(~isnan(seg))
            out(i) = nanmean(seg);
        end
    end
end


function out = rolling_corr(a, b, win)
% ROLLING_CORR  Centred rolling Pearson correlation of two equal-length vectors.
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


function [roll_mean, roll_std] = rolling_circ_stats(phi, win)
% ROLLING_CIRC_STATS  Centred rolling circular mean and std of a phase vector.
%
%   phi must be in radians.  Output length equals input length.
%
%   Circular mean : angle( mean( exp(i*phi) ) )   within each window
%   Circular std  : sqrt( -2 * ln( R ) )          where R = mean resultant length
%                     -> 0 for perfectly concentrated, sqrt(2) for uniform
    n         = numel(phi);
    roll_mean = NaN(n, 1);
    roll_std  = NaN(n, 1);
    half      = floor(win / 2);

    for i = 1:n
        i0  = max(1, i - half);
        i1  = min(n, i + (win - half - 1));
        seg = phi(i0:i1);
        valid = ~isnan(seg);
        if sum(valid) >= 2
            z            = exp(1i * seg(valid));
            R            = abs(mean(z));
            roll_mean(i) = angle(mean(z));
            roll_std(i)  = sqrt(-2 * log(max(R, eps)));
        end
    end
end
