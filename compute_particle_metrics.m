function particle_metrics = compute_particle_metrics(traces, win)
% COMPUTE_PARTICLE_METRICS  Compute per-particle position statistics and speed.
%
%   particle_metrics = compute_particle_metrics(traces, exposure_time_s, win)
%
%   INPUT
%   -----
%   traces          - {nParticles x 1} cell, each [nFrames x 4] (x,y,z,t)
%   exposure_time_s - scalar, camera exposure time in seconds
%   win             - scalar integer, rolling window size in frames
%
%   OUTPUT
%   ------
%   particle_metrics - {nParticles x 1} cell.  Each element is a struct:
%
%       .x, .y, .t          raw coordinate and time vectors [nFrames x 1]
%       .r                  radial distance from instantaneous CoM [nFrames x 1]
%       .theta              azimuthal angle from CoM [nFrames x 1]  (radians)
%
%       .mean_x, .mean_y    time-averaged position (scalar)
%       .std_x,  .std_y     standard deviation of x and y (scalar)
%       .std_r              std of radial distance r (scalar)
%
%       .speed_x            |dx/dt| frame-to-frame speed in x [nFrames-1 x 1]
%       .speed_y            |dy/dt| frame-to-frame speed in y [nFrames-1 x 1]
%       .speed_r            |dr/dt| frame-to-frame speed in r [nFrames-1 x 1]
%       .speed_total        sqrt(vx^2 + vy^2) total speed [nFrames-1 x 1]
%       .mean_speed_total   scalar
%
%       .roll_std_x         rolling std of x [nFrames x 1]  (NaN at edges)
%       .roll_std_y         rolling std of y [nFrames x 1]
%       .roll_std_r         rolling std of r [nFrames x 1]
%
%   RADIAL COORDINATE DEFINITION
%   -----------------------------
%   The origin for (r, theta) is the instantaneous centre-of-mass (CoM) of
%   ALL particles in the movie at EACH frame.  This is computed in
%   compute_pairwise_metrics AFTER all particle_metrics are computed, because
%   the CoM needs the positions of all particles simultaneously.
%   Here, a PRELIMINARY r is calculated relative to the time-averaged mean
%   position of each individual particle — adequate for single-particle
%   movies and as an initialisation step for multi-particle movies.
%   The pairwise function will OVERWRITE .r and .theta if it runs.
%
%   TIME VECTOR
%   -----------
%   The t column of each trace is used if it is monotonically increasing
%   (i.e. from a reliable clock).  If t is constant, all-zero, or
%   non-monotonic, the exposure time is used to construct a synthetic t.
%
%   TRACE LENGTH MISMATCH
%   ----------------------
%   Traces can have different lengths across particles in the same movie.
%   Each particle is analysed on its own length; paired analyses
%   (compute_pairwise_metrics) will trim to the shortest common length.
%
%   NaN HANDLING
%   ------------
%   NaN values in x or y are preserved in the vectors but excluded from
%   scalar summary statistics (mean, std) and speed calculations.

    n_particles = numel(traces);
    particle_metrics = cell(n_particles, 1);

    % =====================================================================
    %  Step 1: Extract raw coordinates and build time vectors
    % =====================================================================
    all_x = cell(n_particles, 1);
    all_y = cell(n_particles, 1);
    all_t = cell(n_particles, 1);

    for p = 1:n_particles
        
        tr = traces{p};
        if istable(tr)
            tr = [tr.col, tr.row, tr.z, tr.t];
        end
        x  = tr(:, 1);
        y  = tr(:, 2);
        % z = tr(:, 3);  -- discarded
        t_raw = tr(:, 4);

        % Decide whether the t column is trustworthy
        t_diff = diff(t_raw);
        t = t_raw;
        

        all_x{p} = x;
        all_y{p} = y;
        all_t{p} = t;
    end

    % =====================================================================
    %  Step 2: Compute instantaneous CoM across all particles per frame
    %          (clipped to shortest trace for CoM calculation only)
    % =====================================================================
    min_len = min(cellfun(@numel, all_x));

    com_x = zeros(min_len, 1);
    com_y = zeros(min_len, 1);
    for p = 1:n_particles
        com_x = com_x + all_x{p}(1:min_len);
        com_y = com_y + all_y{p}(1:min_len);
    end
    com_x = com_x / n_particles;
    com_y = com_y / n_particles;

    % =====================================================================
    %  Step 3: Per-particle derived quantities
    % =====================================================================
    for p = 1:n_particles
        x = all_x{p};
        y = all_y{p};
        t = all_t{p};
        nF = numel(x);

        % --- Radial / angular coordinate relative to instantaneous CoM --
        % For frames beyond min_len, fall back to per-particle mean position
        r     = NaN(nF, 1);
        theta = NaN(nF, 1);

        n_com = min(nF, min_len);   % frames with valid CoM
        dx_com = x(1:n_com) - com_x(1:n_com);
        dy_com = y(1:n_com) - com_y(1:n_com);
        r(1:n_com)     = sqrt(dx_com.^2 + dy_com.^2);
        theta(1:n_com) = atan2(dy_com, dx_com);

        if nF > min_len
            % Extra frames: use individual particle mean as fallback origin
            x0_ind = nanmean(x);
            y0_ind = nanmean(y);
            for f = min_len+1 : nF
                dx = x(f) - x0_ind;
                dy = y(f) - y0_ind;
                r(f)     = sqrt(dx^2 + dy^2);
                theta(f) = atan2(dy, dx);
            end
        end

        % --- Scalar position statistics ----------------------------------
        mean_x = nanmean(x);
        mean_y = nanmean(y);
        std_x  = nanstd(x);
        std_y  = nanstd(y);
        std_r  = nanstd(r);

        % --- Frame-to-frame speed ----------------------------------------
        dt       = diff(t);
        
        dx = diff(x);   dx(isnan(dx)) = 0;
        dy = diff(y);   dy(isnan(dy)) = 0;
        dr = diff(r);   dr(isnan(dr)) = 0;

        speed_x     = abs(dx) ./ dt;
        speed_y     = abs(dy) ./ dt;
        speed_r     = abs(dr) ./ dt;
        speed_total = sqrt((dx./dt).^2 + (dy./dt).^2);

        mean_speed_total = nanmean(speed_total);

        % --- Rolling-window standard deviation ---------------------------
        roll_std_x = rolling_nanstd(x, win);
        roll_std_y = rolling_nanstd(y, win);
        roll_std_r = rolling_nanstd(r, win);

        % --- Store -------------------------------------------------------
        particle_metrics{p} = struct( ...
            'x',                x,                  ...
            'y',                y,                  ...
            't',                t,                  ...
            'r',                r,                  ...
            'theta',            theta,              ...
            'mean_x',           mean_x,             ...
            'mean_y',           mean_y,             ...
            'std_x',            std_x,              ...
            'std_y',            std_y,              ...
            'std_r',            std_r,              ...
            'speed_x',          speed_x,            ...
            'speed_y',          speed_y,            ...
            'speed_r',          speed_r,            ...
            'speed_total',      speed_total,        ...
            'mean_speed_total', mean_speed_total,   ...
            'roll_std_x',       roll_std_x,         ...
            'roll_std_y',       roll_std_y,         ...
            'roll_std_r',       roll_std_r          );
    end
end


% =========================================================================
%  LOCAL UTILITY: rolling_nanstd
% =========================================================================
function out = rolling_nanstd(x, win)
% Compute rolling standard deviation with window size win.
% Output is same length as x; edges where the full window is not yet
% available are set to NaN.
%   - Centred window is used when possible.
%   - Requires Statistics and Machine Learning Toolbox for movstd;
%     falls back to a manual loop if not available.

    n   = numel(x);
    out = NaN(n, 1);
    half = floor(win / 2);

    for i = 1:n
        i_start = max(1, i - half);
        i_end   = min(n, i + (win - half - 1));
        seg     = x(i_start : i_end);
        if sum(~isnan(seg)) >= 2   % need at least 2 points for std
            out(i) = nanstd(seg);
        end
    end
end
