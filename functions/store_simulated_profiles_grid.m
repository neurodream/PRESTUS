function [] = store_simulated_profiles_grid(parameters_fname, fpath, fname, depths, widths)

parameters = load_parameters(parameters_fname);

% Initialize figure and progress bar
figure;

i = 1; % Initialize subplot index
for depth = depths
    for width = widths
        % Perform calculation
        [parameters, axial_position, i_axial_oneil] = calculate_transducer_phases(parameters, 1, depth, width, 100, false);
        

        % Plot in subplot
        subplot(numel(depths), numel(widths), i);
        plot(axial_position, i_axial_oneil);
        ylim([0 120]); % Apply ylim to each subplot

        % Update loop variables
        i = i + 1;
    end
end

savefig(fullfile(fpath, fname));

end