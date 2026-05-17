function [Results, tables] = main_analyze_optical_trapping()
% MAIN_ANALYZE_OPTICAL_TRAPPING
%   Entry point for the optical trapping / optical matter analysis pipeline.
%
%   WORKFLOW
%   --------
%   1. User selects folder containing TrackRes.mat
%   2. User provides experimental metadata (exposure time, wavelength, window)
%   3. Traces are loaded, validated, and grouped by movie
%   4. Per-particle metrics are computed (position stats, speed)
%   5. Pairwise metrics are computed for multi-particle movies
%      (interparticle distance, motion correlations)
%   6. Rolling-window statistics are computed for all metrics
%   7. Movie-level and dataset-level averages are assembled into tables
%   8. Results are saved as .mat and .csv in the selected folder
%
%   OUTPUT
%   ------
%   Results  -  struct with fields:
%       .meta       experimental metadata + pipeline settings
%       .movies     per-movie analysis structs (cell array)
%       .summary    dataset-level aggregated statistics
%       .tables     MATLAB tables (movie summary, pairwise summary)
%
%   DESIGN NOTES
%   ------------
%   - z-coordinate is loaded but NOT used in any computation
%   - Radial coordinate r is defined relative to the center-of-mass of all
%     particles in a given movie at each frame (see compute_particle_metrics)
%   - Rolling-window correlation uses Pearson correlation of instantaneous
%     displacements (dx, dy, dr) within each window
%   - Speed is computed as |delta_position| / delta_t, using the t column
%     of each trace if available; falls back to exposure_time if t is
%     constant or missing
%   - Particle pairs are named "P1-P2", "P1-P3", etc. consistently
%
%   See also: load_tracking_data, group_traces_by_movie,
%             compute_particle_metrics, compute_pairwise_metrics,
%             compute_rolling_metrics, build_summary_tables, save_results

% =========================================================================
%  1. FOLDER SELECTION
% =========================================================================
folder = uigetdir(pwd, 'Select folder containing TrackRes.mat');
if isequal(folder, 0)
    error('OT:UserCancelled', 'No folder selected. Analysis aborted.');
end

% =========================================================================
%  2. USER INPUTS  (exposure time, wavelength, rolling window)
% =========================================================================
% =========================================================================
%  2. EXPERIMENTAL PARAMETERS  — edit these before running
% =========================================================================
meta.exposure_time_s       = 0.05;    % [s]   exposure time per frame
meta.wavelength_nm         = 975;     % [nm]  laser wavelength
meta.rolling_window_frames = 15;      % [frames] rolling window size
meta.folder                = folder;
meta.analysis_date         = datetime("now");

fprintf('\n=== Optical Trapping Analysis Pipeline ===\n');
fprintf('Folder      : %s\n', folder);
fprintf('Exposure    : %.4f s\n', meta.exposure_time_s);
fprintf('Wavelength  : %d nm\n', meta.wavelength_nm);
fprintf('Roll window : %d frames\n\n', meta.rolling_window_frames);

% =========================================================================
%  3. LOAD DATA
% =========================================================================

trackRes = load_tracking_data(folder, meta.exposure_time_s);
% Traces are in nm, time is in seconds


% =========================================================================
%  4. GROUP TRACES BY MOVIE
% =========================================================================
movies_raw = group_traces_by_movie(trackRes);
n_movies   = numel(movies_raw);
fprintf('Found %d movie(s).\n\n', n_movies);

% =========================================================================
%  5. PER-MOVIE ANALYSIS
% =========================================================================
movie_results = cell(n_movies, 1);

for m = 1:n_movies
    mov = movies_raw{m};
    n_particles = numel(mov.traces);    % number of particles in this movie

    fprintf('[Movie %d]  %d particle(s) ...\n', mov.movie_idx, n_particles);

    % --- Per-particle metrics (positions, std, speed) --------------------
    particle_metrics = compute_particle_metrics( ...
        mov.traces, meta.rolling_window_frames);
    
    particle_metrics = compute_angular_metrics(particle_metrics, meta.rolling_window_frames);
    
    % --- Pairwise metrics (distance, correlation) -----------------------
    if n_particles > 1
        pairwise_metrics = compute_pairwise_metrics( ...
            particle_metrics, meta.rolling_window_frames);
        sync_metrics = compute_synchronisation_index( ...
            particle_metrics, meta.rolling_window_frames);
    else
        pairwise_metrics = struct();   % empty for single-particle movies
        sync_metrics = strut(); % empty for single-particle movies
    end

    
    % --- Assemble movie struct -------------------------------------------
    movie_results{m} = struct( ...
        'movie_idx',        mov.movie_idx,      ...
        'n_particles',      n_particles,        ...
        'particle_metrics', {particle_metrics},   ...
        'sync_metrics', sync_metrics,   ...
        'pairwise_metrics', pairwise_metrics    );
end

% =========================================================================
%  6. DATASET-LEVEL SUMMARY
% =========================================================================
summary = build_dataset_summary(movie_results);

% =========================================================================
%  7. TABLES
% =========================================================================
tables = build_summary_tables(movie_results, summary);

% =========================================================================
%  8. ASSEMBLE & SAVE
% =========================================================================
Results = struct( ...
    'meta',    meta,          ...
    'movies',  {movie_results}, ...
    'summary', summary,       ...
    'tables',  tables         );

% =========================================================================
%  9. PLOTS  — full diagnostic report for the first movie
%             Change the index or wrap in a loop to inspect more movies
% =========================================================================
movie_to_plot = 1;   % <-- change this to plot a different movie
fprintf('Plotting full report for movie index %d ...\n', ...
    Results.movies{movie_to_plot}.movie_idx);
plot_full_movie_report(Results.movies{movie_to_plot}, Results.meta);

% =========================================================================
% 10. COMMAND-LINE SUMMARY  — mean values for every assembly
% =========================================================================
print_assembly_summary(Results);

% =========================================================================
% 11. SAVE
% =========================================================================
save_results(Results, folder);
fprintf('\nDone. Results saved to:\n  %s\n', folder);

T = Results.tables.movie_summary;

figure; hold on
for m = 1:height(T)
    th = T.MeanThetaUnwrap{m};
    t  = (0:numel(th)-1) * meta.exposure_time_s;
    plot(t, th / (2*pi),'DisplayName', sprintf('Mov %d', m))   % full turns
    
end
xlabel('Time (s)');  ylabel('Mean rotation (full turns)');
legen('show')

end

% =========================================================================
%  LOCAL HELPER:  prompt_user_inputs
% =========================================================================
function meta = prompt_user_inputs()
% Prompt the user for the three required experimental parameters.
% Defaults are provided so pressing Enter keeps a sensible value.

    meta = struct();

    % --- Exposure time ---------------------------------------------------
    raw = input('Enter exposure time per frame [s] (e.g. 0.02): ', 's');
    raw = strtrim(raw);
    if isempty(raw)
        meta.exposure_time_s = 0.02;
        fprintf('  -> Using default: 0.02 s\n');
    else
        val = str2double(raw);
        if isnan(val) || val <= 0
            error('OT:BadInput', 'Exposure time must be a positive number.');
        end
        meta.exposure_time_s = val;
    end

    % --- Laser wavelength ------------------------------------------------
    raw = input('Enter laser wavelength [nm] (e.g. 532): ', 's');
    raw = strtrim(raw);
    if isempty(raw)
        meta.wavelength_nm = 532;
        fprintf('  -> Using default: 532 nm\n');
    else
        val = str2double(raw);
        if isnan(val) || val <= 0
            error('OT:BadInput', 'Wavelength must be a positive number.');
        end
        meta.wavelength_nm = val;
    end

    % --- Rolling window --------------------------------------------------
    raw = input('Enter rolling window size [frames] (default 15): ', 's');
    raw = strtrim(raw);
    if isempty(raw)
        meta.rolling_window_frames = 15;
        fprintf('  -> Using default: 15 frames\n');
    else
        val = round(str2double(raw));
        if isnan(val) || val < 2
            error('OT:BadInput', 'Rolling window must be >= 2 frames.');
        end
        meta.rolling_window_frames = val;
    end
end
