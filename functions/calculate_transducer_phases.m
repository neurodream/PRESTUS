function [parameters, axial_position, i_axial_oneil] = calculate_transducer_phases(parameters, transducer_ind, expected_focal_distance_mm, ROI_width_mm, goal_intensity, sham)

%% optimize phases (active condition)

transducer = parameters.transducers(transducer_ind);

axial_position = (1:parameters.default_grid_dims(3))*(parameters.grid_step_mm);
axial_position = axial_position';
dist_to_exit_plane = transducer.curv_radius_mm - transducer.dist_to_plane_mm;
ax_pos = axial_position + dist_to_exit_plane; % TODO sanity check if distance to exit plane makes sense - deleted for now

opt_limits = [axial_position(2,1), axial_position(end,1)];

% % create the weights function
% [flhm_center, flhm_center_index] = get_flhm_center_position( ...
%     axial_position, ...
%     desired_function ...
%     );
% weights = normpdf_nostatstoolbox( ...
%     axial_position(:,1), ...
%     axial_position(flhm_center_index,1), ...
%     axial_position(flhm_center_index,1)/.5 ...
%     );

desired_function = create_boxcar(expected_focal_distance_mm, ROI_width_mm, axial_position, 0, goal_intensity);

% temporarily set a field "transducer" to be able to use the old function
parameters.transducer = transducer;

optimize_phases = @(phases_and_velocity) phase_optimization_annulus_full_curve( ...
    phases_and_velocity(1:transducer.n_elements-1), ...
    parameters, ...
    phases_and_velocity(transducer.n_elements),...
    axial_position, ...
    desired_function, ...
    0, ...
    opt_limits, ...
    ones(1, numel(axial_position)) ...
);

rng(100,'twister') % setting seed for consistency

velocity = transducer.source_amp(1)/(parameters.medium.water.density*parameters.medium.water.sound_speed);   % [m/s]
func = optimize_phases;
x0 = [randi(360, [1 transducer.n_elements - 1])/180*pi velocity];
lb = zeros(1,transducer.n_elements);
ub = [2*pi*ones(1,transducer.n_elements - 1) 0.2];
options = setoptimoptions( ...
    'popsize',1000, ...
    'FinDiffType', 'central', ...
    'MaxFunEvals', 1e6, ...
    'MaxIter', 1e7 ...
    );
[opt_phases_and_velocity, ~, ~, ~] = minimize(func, x0, [],[],[],[],lb, ub, [], options);

opt_phases = [0 opt_phases_and_velocity(1:9)];
opt_velocity = opt_phases_and_velocity(end);

parameters.transducers(transducer_ind).source_phase_rad = opt_phases;
% TODO set velocity in Pascal

% note: I adjusted the phase_optimization_annulus_full_curve function so
% that it returns the axial intensity as the last output parameter
p_axial_oneil = focusedAnnulusONeil( ...
    transducer.curv_radius_mm/1e3, ...
    [transducer.Elements_ID_mm; transducer.Elements_OD_mm]/1e3, ...
    repmat(opt_velocity, 1, transducer.n_elements), ...
    opt_phases, ...
    transducer.source_freq_hz, ...
    parameters.medium.water.sound_speed, ...
    parameters.medium.water.density, ...
    (axial_position-0.5)*1e-3 ... % before here was ax_pos
    );

i_axial_oneil = p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;

% figure;
% i_axial_oneil_active = i_axial_oneil;
% plot(axial_position, i_axial_oneil_active)
% hold on;

if sham

    % automatic detection of cutoff
    minima_indices = islocalmin(i_axial_oneil);
    min_dists = axial_position(minima_indices);
    % Find the closest value in min_dists that is lower than expected_focal_distance_mm
    lower_values = min_dists(min_dists < expected_focal_distance_mm);  % Filter values in l that are less than expected focal distance
    cutoff = max(lower_values);  % Find the maximum of the filtered values

    desired_function = i_axial_oneil .* single(axial_position < cutoff); % near field

    optimize_phases = @(phases_and_velocity) phase_optimization_annulus_full_curve( ...
        phases_and_velocity(1:transducer.n_elements-1), ...
        parameters, ...
        phases_and_velocity(transducer.n_elements),...
        axial_position, ...
        desired_function, ...
        0, ...
        opt_limits, ...
        ones(1, numel(axial_position)) ...
    );
    
    rng(100,'twister') % setting seed for consistency
    
    velocity = transducer.source_amp(1)/(parameters.medium.water.density*parameters.medium.water.sound_speed);   % [m/s]
    func = optimize_phases;
    x0 = [randi(360, [1 transducer.n_elements - 1])/180*pi velocity];
    lb = zeros(1,transducer.n_elements);
    ub = [2*pi*ones(1,transducer.n_elements - 1) 0.2];
    options = setoptimoptions( ...
        'popsize',1000, ...
        'FinDiffType', 'central', ...
        'MaxFunEvals', 1e6, ...
        'MaxIter', 1e7 ...
        );
    [opt_phases_and_velocity, ~, ~, ~] = minimize(func, x0, [],[],[],[],lb, ub, [], options);
    
    opt_phases = [0 opt_phases_and_velocity(1:9)];
    opt_velocity = opt_phases_and_velocity(end);
    
    parameters.transducers(transducer_ind).source_phase_rad = opt_phases;
    % TODO set velocity in Pascal
    
    % note: I adjusted the phase_optimization_annulus_full_curve function so
    % that it returns the axial intensity as the last output parameter
    p_axial_oneil = focusedAnnulusONeil( ...
        transducer.curv_radius_mm/1e3, ...
        [transducer.Elements_ID_mm; transducer.Elements_OD_mm]/1e3, ...
        repmat(opt_velocity, 1, transducer.n_elements), ...
        opt_phases, ...
        transducer.source_freq_hz, ...
        parameters.medium.water.sound_speed, ...
        parameters.medium.water.density, ...
        (axial_position-0.5)*1e-3 ... % before here was ax_pos
        );
    
    i_axial_oneil = p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;

end

% plot(axial_position, i_axial_oneil);

parameters = rmfield(parameters, 'transducer');



% [~, ~, ~, h, i_axial_oneil] = phase_optimization_annulus_full_curve( ...
%     opt_phases_and_velocity(1:parameters.transducer.n_elements-1), ...
%     parameters, ...
%     opt_phases_and_velocity(parameters.transducer.n_elements),...
%     axial_position, ...
%     desired_function, ...
%     1, ...
%     opt_limits, ...
%     weights ...
%     );  
% 
% close(h);

% TODO handle if first element should be 0 or not

%% postprocessing
% degs = rad2deg([0 opt_phases_and_velocity(1:9)]);
% degs_str = arrayfun(@(x) sprintf('%.3f', x), degs, 'UniformOutput', false);
% fprintf('[%s]\n', strjoin(degs_str, ', '));
% 
% disp(['delta focus: ' num2str(axial_position(find(i_axial_oneil == max(i_axial_oneil))) - expected_focal_distance_mm)]);
% disp(['delta intensity: ' num2str(max(i_axial_oneil) - goal_intensity)]);
% 
% % intermediate plot to check where (if) there should be a cutoff
% h = figure;
% plot(axial_position, i_axial_oneil);
% hold on;
% plot(axial_position, desired_function/2);

end