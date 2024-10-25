function update_transducers_and_run(subject_id, parameters, ID, vary)

    % Loop through each transducer and apply the updates
    for i = 1:numel(parameters.transducers)
        % Get current position and focus for the transducer
        t = parameters.transducers(i).pos_t1_grid;
        f = parameters.transducers(i).focus_pos_t1_grid;
        
        vt = round(-3 + 6 * rand(1, 3));  % Random values for pos_t1_grid
        vf = round(-3 + 6 * rand(1, 3));  % Random values for focus_pos_t1_grid

        % Conditionally set vt and vf based on the 'vary' parameter
        switch vary
            case 'transducer'
                vf = [0, 0, 0];            % No variation in focus_pos_t1_grid
            case 'focus'
                vt = [0, 0, 0];            % No variation in pos_t1_grid
            case 'none'
                vf = [0, 0, 0];
                vt = [0, 0, 0];
        end
        
        % Update pos_t1_grid and focus_pos_t1_grid with the calculated values
        parameters.transducers(i).pos_t1_grid = t + vt;
        parameters.transducers(i).focus_pos_t1_grid = f + vf;
    end
    
    % Set the results filename
    parameters.results_filename_affix = ID;
    
    % Run the pipeline
    single_subject_pipeline_with_slurm(subject_id, parameters);
    % single_subject_pipeline(subject_id, parameters);
end