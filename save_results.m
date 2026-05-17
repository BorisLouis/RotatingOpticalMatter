function save_results(Results, folder)
% SAVE_RESULTS  Persist the analysis output to disk.
%
%   save_results(Results, folder)
%
%   Saves:
%     AnalysisResults.mat          full Results struct
%     table_movies.csv             movie-level summary table
%     table_particles.csv          per-particle detail table
%     table_pairwise.csv           pairwise statistics (if present)
%     table_dataset.csv            dataset-level one-row summary
%
%   Existing files with the same name in the folder are overwritten.

    % --- .mat -----------------------------------------------------------
    mat_path = fullfile(folder, 'AnalysisResults.mat');
    save(mat_path, 'Results', '-v7.3');
    fprintf('  Saved: %s\n', mat_path);

    % --- CSV helpers ----------------------------------------------------
    function write_csv(tbl, name)
        if isempty(tbl), return; end
        p = fullfile(folder, name);
        try
            writetable(tbl, p);
            fprintf('  Saved: %s\n', p);
        catch ME
            warning('OT:CSVWriteFail', 'Could not write %s: %s', p, ME.message);
        end
    end

    t = Results.tables;
    write_csv(t.movie_summary,    'table_movies.csv');
    write_csv(t.particle_detail,  'table_particles.csv');
    write_csv(t.dataset_summary,  'table_dataset.csv');

    if isfield(t, 'pairwise_summary')
        write_csv(t.pairwise_summary, 'table_pairwise.csv');
    end
end
