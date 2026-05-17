function sync_metrics = compute_synchronisation_index(particle_metrics, win)
% COMPUTE_SYNCHRONISATION_INDEX  Kuramoto order parameter for multi-particle movies.
%
%   sync_metrics = compute_synchronisation_index(particle_metrics, angular_metrics, win)
%
%   Quantifies the degree of phase synchronisation across all particles
%   in a movie using the Kuramoto order parameter R(t).
%
%   INPUT
%   -----
%   particle_metrics  - struct array (output of compute_particle_metrics,
%                       after pairwise CoM correction)
%   angular_metrics   - struct array (output of compute_angular_metrics)
%   win               - rolling window size in frames
%
%   OUTPUT
%   ------
%   sync_metrics  - struct with fields:
%
%     .R_instantaneous   Kuramoto R at each frame [nFrames x 1], in [0,1]
%     .Psi_instantaneous mean phase angle at each frame [nFrames x 1] (rad)
%     .mean_R            time-averaged R  (scalar, the synchronisation index)
%     .std_R             std of R over time
%     .roll_R            rolling mean of R  [nFrames x 1]
%
%     .pairwise_phase_diff   {nPairs x 1} cell, each [nFrames x 1]
%                            instantaneous phase difference phi_i - phi_j (rad)
%                            unwrapped for continuity
%     .pairwise_labels       {nPairs x 1} cell of strings, e.g. 'P1-P2'
%     .pairwise_mean_dphi    [nPairs x 1]  mean phase difference per pair
%     .pairwise_std_dphi     [nPairs x 1]  std of phase difference per pair
%                            Small std_dphi → phase-locked pair.
%
%     .phase_locking_index   [nPairs x 1]  |mean(exp(i*(phi_i - phi_j)))|
%                            1 = perfectly locked, 0 = incoherent
%
%   DEFINITION OF KURAMOTO ORDER PARAMETER
%   ----------------------------------------
%   For N particles with instantaneous phases theta_1(t) ... theta_N(t):
%
%       R(t) * exp(i*Psi(t)) = (1/N) * sum_j( exp(i * theta_j(t)) )
%
%   R(t) in [0,1]:
%       R = 1  →  all particles at the same angular phase (full sync)
%       R = 0  →  phases uniformly distributed (incoherent)
%
%   The phases used here are the UNWRAPPED theta values from angular_metrics
%   (mod 2*pi to bring back to the unit circle for the complex exponential).
%
%   PAIRWISE PHASE LOCKING INDEX (PLI)
%   ------------------------------------
%   PLI_ij = | mean_t( exp(i * (theta_i(t) - theta_j(t))) ) |
%   Equivalent to the length of the mean resultant vector of the phase
%   difference distribution.  PLI = 1 means the phase difference is
%   constant (locked); PLI = 0 means the difference drifts uniformly.
%
%   NOTE
%   ----
%   For a single-particle movie this function returns R = 1 trivially
%   (a single phasor always has unit magnitude). It should only be
%   interpreted for N >= 2.

    n_particles = numel(particle_metrics);
   

    if n_particles < 2
        warning('OT:SingleParticle', ...
            'Synchronisation index is trivially 1 for a single particle.');
    end
    

    % =====================================================================
    %  1.  Collect unwrapped phases, map to [0, 2*pi) for Kuramoto
    % =====================================================================
    % Use the shortest common trace length
    nF = numel(particle_metrics{1}.theta);      % Build phase matrix [nF x nParticles]
    % Kuramoto needs phases on the unit circle → use mod(theta_unwrap, 2*pi)
    phase_mat = NaN(nF, n_particles);
    for p = 1:n_particles
        th_uw = particle_metrics{p}.theta_unwrap(1:nF);
        phase_mat(:, p) = mod(th_uw, 2*pi);
    end

    % =====================================================================
    %  2.  Instantaneous Kuramoto order parameter R(t) and mean phase Psi(t)
    % =====================================================================
    % Z(t) = mean over particles of exp(i * theta_p(t))
    Z = mean(exp(1i * phase_mat), 2, 'omitnan');   % [nF x 1] complex

    R_inst   = abs(Z);       % order parameter
    Psi_inst = angle(Z);     % mean phase angle

    mean_R = nanmean(R_inst);
    std_R  = nanstd(R_inst);
    roll_R = rolling_nanmean_vec(R_inst, win);

    % =====================================================================
    %  3.  Pairwise phase difference and phase locking index
    % =====================================================================
    pair_labels  = {};
    pair_dphi    = {};
    mean_dphi    = [];
    std_dphi     = [];
    pli          = [];

    for i = 1:n_particles
        for j = (i+1):n_particles
            label = sprintf('P%d-P%d', i, j);

            % Use unwrapped phases for difference (gives smoother signal)
            phi_i = particle_metrics{i}.theta_unwrap(1:nF);
            phi_j = particle_metrics{j}.theta_unwrap(1:nF);

            dphi        = phi_i - phi_j;
            dphi_unwrap = unwrap(dphi);   % remove 2*pi jumps in difference

            % Phase locking index: |<exp(i*dphi)>|  (use wrapped dphi)
            pli_val = abs(nanmean(exp(1i * dphi)));

            pair_labels{end+1}  = label;         %#ok<AGROW>
            pair_dphi{end+1}    = dphi_unwrap;   %#ok<AGROW>
            mean_dphi(end+1)    = nanmean(dphi_unwrap); %#ok<AGROW>
            std_dphi(end+1)     = nanstd(dphi_unwrap);  %#ok<AGROW>
            pli(end+1)          = pli_val;               %#ok<AGROW>
        end
    end

    % =====================================================================
    %  4.  Store
    % =====================================================================
    sync_metrics = struct( ...
        'R_instantaneous',      R_inst,          ...
        'Psi_instantaneous',    Psi_inst,        ...
        'mean_R',               mean_R,          ...
        'std_R',                std_R,           ...
        'roll_R',               roll_R,          ...
        'pairwise_phase_diff',  {pair_dphi},     ...
        'pairwise_labels',      {pair_labels},   ...
        'pairwise_mean_dphi',   mean_dphi(:),    ...
        'pairwise_std_dphi',    std_dphi(:),     ...
        'phase_locking_index',  pli(:)           );
end


% =========================================================================
%  LOCAL UTILITY
% =========================================================================
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
