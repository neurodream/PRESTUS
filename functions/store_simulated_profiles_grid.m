function [] = store_simulated_profiles_grid(parameters_fname, fpath, fname, depths, widths)

    parameters = load_parameters(parameters_fname);

    % Initialize figure and progress bar
    figure;
    goal_intensity = 100;
    total_iterations = numel(depths) * numel(widths);  % Total iterations for progress bar
    progress = waitbar(0, 'Processing...');  % Initialize progress bar

    i = 1; % Initialize subplot index
    current_iteration = 0; % Initialize iteration counter

    for depth = depths
        for width = widths
            % Perform calculation
            [parameters, axial_position, i_axial_oneil] = calculate_transducer_phases(parameters, 1, depth, width, goal_intensity, false);
            desired_function = create_boxcar(depth, width, axial_position, 0, goal_intensity);
            % desired_function = create_mexican_hat(depth, width, axial_position, 0, goal_intensity);
            % desired_function = create_gaussian(depth, width, axial_position, 0, goal_intensity);
            
            % Adjust the source_amp of contralateral-focused transducer if needed
            peaks_list = findpeaks(i_axial_oneil); 
            peak = peaks_list(end);
            
            % TODO should not be needed anymore when implemented in
            % calculate_transducer_phases function
            if goal_intensity > peak
                i_axial_oneil = i_axial_oneil * (goal_intensity / peak);
            end

            % Plot in subplot
            subplot(numel(depths), numel(widths), i);
            plot(axial_position, i_axial_oneil);
            hold on;
            plot(axial_position, desired_function);
            ylim([0 120]);  % Apply ylim to each subplot

            % Update progress bar
            current_iteration = current_iteration + 1;
            waitbar(current_iteration / total_iterations, progress, ...
                    sprintf('Processing... %.0f%%', (current_iteration / total_iterations) * 100));

            % Update loop variables
            i = i + 1;
        end
    end

    % Close the progress bar after the loop is complete
    close(progress);

    % Save the figure
    savefig(fullfile(fpath, fname));

end
