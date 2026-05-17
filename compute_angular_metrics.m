function particle_metrics = compute_angular_metrics(particle_metrics, exposure_time_s, win)
% COMPUTE_ANGULAR_METRICS  Angular speed and rotation frequency per particle.
%
%   angular_metrics = compute_angular_metrics(particle_metrics, exposure_time_s, win)
%
%   Operates on the .theta and .t fields already computed in
%   compute_particle_metrics (azimuthal angle relative to instantaneous CoM).
%   Should be called AFTER compute_pairwise_metrics, because pairwise
%   overwrites .theta with the CoM-corrected version.
%
%   INPUT
%   -----
%   particle_metrics  - struct array from compute_particle_metrics
%                       (after pairwise CoM correction)
%   exposure_time_s   - scalar, used as fallback dt
%   win               - rolling window size in frames
%
%   OUTPUT
%   ------
%   angular_metrics   - struct array (one element per particle) with fields:
%
%     .theta_unwrap     unwrapped azimuthal angle [nFrames x 1]  (rad)
%     .omega            instantaneous angular speed [nFrames-1 x 1] (rad/s)
%                       signed: positive = CCW, negative = CW
%     .omega_abs        |omega|  (rad/s)
%     .mean_omega       time-averaged signed angular velocity (rad/s)
%     .mean_omega_abs   time-averaged angular speed (rad/s)
%     .std_omega        std of instantaneous omega
%
%     .roll_mean_omega  rolling mean of signed omega [nFrames-1 x 1]
%     .roll_std_omega   rolling std of omega
%
%     .rot_freq_mean    rotation frequency from linear fit to theta_unwrap
%                       (Hz).  Positive = net CCW, negative = net CW.
%                       This is the most robust estimator of steady rotation.
%     .rot_period_s     1 / |rot_freq_mean|  in seconds  (Inf if ~0)
%
%     .fft_freq         frequency axis of FFT  (Hz)  [nFFT x 1]
%     .fft_power        single-sided power spectrum of omega  [nFFT x 1]
%     .dominant_freq_hz dominant frequency from FFT peak  (Hz)
%     .dominant_period_s 1 / dominant_freq_hz  (s)
%
%   DEFINITIONS
%   -----------
%   theta is atan2(dy_from_CoM, dx_from_CoM) in [-pi, pi].
%   Unwrapping removes 2*pi jumps to give a continuous cumulative angle.
%   Angular velocity:  omega(k) = d(theta_unwrap)/dt  at frame k.
%   Rotation frequency (Hz) from linear fit:  theta_unwrap ≈ 2*pi*f*t + phi0
%   so f = slope / (2*pi).
%
%   FFT IS COMPUTED ON omega (not theta) to highlight oscillatory rotation.
%   The DC component (mean omega) dominates the spectrum for steady rotation.

    n_particles = numel(particle_metrics);
   
    for p = 1:n_particles
        pm  = particle_metrics{p};
        t   = pm.t;
        th  = pm.theta;

        % Remove leading/trailing NaNs for unwrap (NaN breaks unwrap)
        valid = ~isnan(th);
        th_valid = th(valid);
        t_valid  = t(valid);

        if numel(th_valid) < 3
            warning('OT:TooFewAngles', ...
                'Particle %d: fewer than 3 valid theta values. Skipping.', p);
            angular_metrics(p) = empty_angular_struct(numel(t));
            continue;
        end

        % ------------------------------------------------------------------
        %  1. Unwrap theta
        % ------------------------------------------------------------------
        th_unwrap_valid = unwrap(th_valid);

        % Map back to full-length vector (NaN where original was NaN)
        th_unwrap = NaN(size(th));
        th_unwrap(valid) = th_unwrap_valid;

        % ------------------------------------------------------------------
        %  2. Instantaneous angular velocity  omega = d(theta)/dt
        % ------------------------------------------------------------------
        dt = diff(t_valid);
        dt(dt == 0) = exposure_time_s;   % safeguard

        dth    = diff(th_unwrap_valid);
        omega_valid = dth ./ dt;          % signed rad/s

        % Map back to full-length-minus-1 vector
        % (same indexing convention as speed vectors in particle_metrics)
        n_full = numel(t);
        omega     = NaN(n_full - 1, 1);
        omega_abs = NaN(n_full - 1, 1);
        valid_idx = find(valid);
        for k = 1:numel(omega_valid)
            % omega_valid(k) sits between valid_idx(k) and valid_idx(k+1)
            % Place it at valid_idx(k) in the output vector
            if valid_idx(k) <= n_full - 1
                omega(valid_idx(k))     = omega_valid(k);
                omega_abs(valid_idx(k)) = abs(omega_valid(k));
            end
        end

        mean_omega     = nanmean(omega);
        mean_omega_abs = nanmean(omega_abs);
        std_omega      = nanstd(omega);

        % ------------------------------------------------------------------
        %  3. Rolling statistics on omega
        % ------------------------------------------------------------------
        roll_mean_omega = rolling_nanmean_vec(omega, win);
        roll_std_omega  = rolling_nanstd_vec(omega, win);

        % ------------------------------------------------------------------
        %  4. Rotation frequency from linear fit  (robust, low-noise)
        % ------------------------------------------------------------------
        % Fit:  theta_unwrap = slope * t + intercept
        % Rotation frequency f = slope / (2*pi)
        coeff = polyfit(t_valid, th_unwrap_valid, 1);
        slope = coeff(1);                   % rad/s
        rot_freq_mean  = slope / (2 * pi);  % Hz
        rot_period_s   = freq_to_period(rot_freq_mean);

        % ------------------------------------------------------------------
        %  5. FFT of omega  (spectral content of angular fluctuations)
        % ------------------------------------------------------------------
        omega_clean = omega_valid;
        omega_clean(isnan(omega_clean)) = 0;   % zero-fill NaNs for FFT

        n_fft   = numel(omega_clean);
        dt_mean = mean(dt);
        fs      = 1 / dt_mean;               % sampling frequency (Hz)

        Y       = fft(omega_clean);
        P2      = abs(Y / n_fft).^2;
        n_half  = floor(n_fft / 2) + 1;
        P1      = P2(1:n_half);
        P1(2:end-1) = 2 * P1(2:end-1);      % single-sided spectrum

        fft_freq = fs * (0:(n_half-1))' / n_fft;

        % Dominant frequency: ignore DC bin (k=1)
        [~, peak_idx] = max(P1(2:end));
        dominant_freq_hz  = fft_freq(peak_idx + 1);
        dominant_period_s = freq_to_period(dominant_freq_hz);

        % ------------------------------------------------------------------
        %  6. Store
        % ------------------------------------------------------------------
        particle_metrics{p}.theta_unwrap      = th_unwrap;
        particle_metrics{p}.omega             = omega;
        particle_metrics{p}.omega_abs         = omega_abs;
        particle_metrics{p}.mean_omega        = mean_omega;
        particle_metrics{p}.mean_omega_abs    = mean_omega_abs;
        particle_metrics{p}.std_omega         = std_omega;
        particle_metrics{p}.roll_mean_omega   = roll_mean_omega;
        particle_metrics{p}.roll_std_omega    = roll_std_omega;
        particle_metrics{p}.rot_freq_mean     = rot_freq_mean;
        particle_metrics{p}.rot_period_s      = rot_period_s;
        particle_metrics{p}.fft_freq          = fft_freq;
        particle_metrics{p}.fft_power         = P1;
        particle_metrics{p}.dominant_freq_hz  = dominant_freq_hz;
        particle_metrics{p}.dominant_period_s = dominant_period_s;
    end
end


% =========================================================================
%  LOCAL UTILITIES
% =========================================================================

function T = freq_to_period(f)
    if abs(f) < 1e-12
        T = Inf;
    else
        T = 1 / abs(f);
    end
end

function s = empty_angular_struct(n)
% Return a struct of NaNs for particles with insufficient data.
    s.theta_unwrap      = NaN(n, 1);
    s.omega             = NaN(n-1, 1);
    s.omega_abs         = NaN(n-1, 1);
    s.mean_omega        = NaN;
    s.mean_omega_abs    = NaN;
    s.std_omega         = NaN;
    s.roll_mean_omega   = NaN(n-1, 1);
    s.roll_std_omega    = NaN(n-1, 1);
    s.rot_freq_mean     = NaN;
    s.rot_period_s      = NaN;
    s.fft_freq          = NaN;
    s.fft_power         = NaN;
    s.dominant_freq_hz  = NaN;
    s.dominant_period_s = NaN;
end

function out = rolling_nanmean_vec(x, win)
    n    = numel(x);
    out  = NaN(n, 1);
    half = floor(win / 2);
    for i = 1:n
        seg = x(max(1,i-half) : min(n, i+(win-half-1)));
        if sum(~isnan(seg)) >= 1
            out(i) = nanmean(seg);
        end
    end
end

function out = rolling_nanstd_vec(x, win)
    n    = numel(x);
    out  = NaN(n, 1);
    half = floor(win / 2);
    for i = 1:n
        seg = x(max(1,i-half) : min(n, i+(win-half-1)));
        if sum(~isnan(seg)) >= 2
            out(i) = nanstd(seg);
        end
    end
end
