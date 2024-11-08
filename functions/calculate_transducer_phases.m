function [parameters, axial_position, i_axial_oneil] = calculate_transducer_phases(parameters, transducer_ind, expected_focal_distance_mm, ROI_width_mm, goal_intensity, sham, desired_function, use_all_phases)

% TODO adjust the velocity to match the intended strength

if nargin < 7
    desired_function = @create_boxcar; % possible functions: create_boxcar create_mexican_hat create_gaussian
    use_all_phases = false;
end

if nargin < 8
    use_all_phases = false;
end

%% optimize phases (active condition)

transducer = parameters.transducers(transducer_ind);

axial_position = (1:parameters.default_grid_dims(3))*(parameters.grid_step_mm);
axial_position = axial_position';
dist_to_exit_plane = transducer.curv_radius_mm - transducer.dist_to_plane_mm;
ax_pos = axial_position + dist_to_exit_plane; % TODO sanity check if distance to exit plane makes sense - deleted for now

opt_limits = [axial_position(2,1), axial_position(end,1)];

ymin = 0;

desired_profile = desired_function(expected_focal_distance_mm, ROI_width_mm, axial_position, ymin, goal_intensity);

% temporarily set a field "transducer" to be able to use the old function
parameters.transducer = transducer;

if use_all_phases
    stop_before = 0;
    multiply_velocity = 0;
    add_velocity = 0.15;
    optimization_function = @phase_optimization_annulus_full_curve_all_elements;
    ub = 2*pi*ones(1,transducer.n_elements);
    x0 = randi(360, [1 transducer.n_elements])/180*pi;
else
    stop_before = 1;
    multiply_velocity = 1;
    add_velocity = 0;
    optimization_function = @phase_optimization_annulus_full_curve;
    ub = [2*pi*ones(1,transducer.n_elements - 1) 0.2];
    x0 = [randi(360, [1 transducer.n_elements - 1])/180*pi 0.15];
end

optimize_phases = @(phases_and_velocity) optimization_function( ...
    phases_and_velocity(1:transducer.n_elements-stop_before), ...
    parameters, ...
    phases_and_velocity(transducer.n_elements)*multiply_velocity + add_velocity,...
    axial_position, ...
    desired_profile, ...
    0, ...
    opt_limits, ...
    ones(1, numel(axial_position)) ... % weights: keep equal
    );

% optimize_phases = @(phases_and_velocity) phase_optimization_annulus( ...
%     phases_and_velocity(1:transducer.n_elements-1), ...
%     parameters, ...
%     phases_and_velocity(transducer.n_elements),...
%     axial_position, ...
%     expected_focal_distance_mm ...
%     );

rng(100,'twister') % setting seed for consistency

% velocity = transducer.source_amp(1)/(parameters.medium.water.density*parameters.medium.water.sound_speed);   % [m/s]
func = optimize_phases;
lb = zeros(1,transducer.n_elements);
options = setoptimoptions( ...
    'popsize',1000, ...
    'FinDiffType', 'central', ...
    'MaxFunEvals', 1e7, ...
    'MaxIter', 1e6 ...
    );
[opt_phases_and_velocity, ~, ~, ~] = minimize(func, x0, [],[],[],[],lb, ub, [], options);

if use_all_phases
    opt_phases = opt_phases_and_velocity;
    opt_velocity = 0.15;
else
    opt_phases = [0 opt_phases_and_velocity(1:9)];
    opt_velocity = opt_phases_and_velocity(end);
end

parameters.transducers(transducer_ind).source_phase_rad = opt_phases;
% TODO set velocity in Pascal

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

% adjust intensity levels to desired intensity
peaks_list = findpeaks(i_axial_oneil); 
peak = peaks_list(end);
adjustment_ratio = max([goal_intensity peak])/min([goal_intensity peak]);
parameters.transducers(transducer_ind).source_amp = transducer.source_amp*sqrt(adjustment_ratio);

% also adjust the intensity so that the sham works correctly
i_axial_oneil = i_axial_oneil*adjustment_ratio;

if sham

    % automatic detection of cutoff
    minima_indices = islocalmin(i_axial_oneil);
    min_dists = axial_position(minima_indices);
    % Find the closest value in min_dists that is lower than expected_focal_distance_mm
    lower_values = min_dists(min_dists < expected_focal_distance_mm);  % Filter values in l that are less than expected focal distance
    cutoff = max(lower_values);  % Find the maximum of the filtered values

    desired_profile = i_axial_oneil .* single(axial_position < cutoff); % near field

    % store the global maximum of the desired sham profile for later
    % intensity adjustment
    peak_desired_sham_profile = max(desired_profile);

    optimize_phases = @(phases_and_velocity) optimization_function( ...
        phases_and_velocity(1:transducer.n_elements-stop_before), ...
        parameters, ...
        phases_and_velocity(transducer.n_elements)*multiply_velocity + add_velocity,...
        axial_position, ...
        desired_profile, ...
        0, ...
        opt_limits, ...
        ones(1, numel(axial_position)) ... % weights: keep equal
        );

    % refresh the optimization function
    func = optimize_phases;

    % optimize_phases = @(phases_and_velocity) phase_optimization_annulus( ...
    %     phases_and_velocity(1:transducer.n_elements-1), ...
    %     parameters, ...
    %     phases_and_velocity(transducer.n_elements),...
    %     axial_position, ...
    %     expected_focal_distance_mm ...
    %     );
    
    % rng(100,'twister') % setting seed for consistency
    
    % velocity = transducer.source_amp(1)/(parameters.medium.water.density*parameters.medium.water.sound_speed);   % [m/s]

    [opt_phases_and_velocity, ~, ~, ~] = minimize(func, x0, [],[],[],[],lb, ub, [], options);
    
    if use_all_phases
        opt_phases = opt_phases_and_velocity;
        opt_velocity = 0.15;
    else
        opt_phases = [0 opt_phases_and_velocity(1:9)];
        opt_velocity = opt_phases_and_velocity(end);
    end
    
    parameters.transducers(transducer_ind).source_phase_rad = opt_phases;
    % TODO set velocity in Pascal
    
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
    
    % adjust intensity levels to desired intensity
    sham_peak = max(i_axial_oneil);
    adjustment_ratio = max([peak_desired_sham_profile sham_peak])/min([peak_desired_sham_profile sham_peak]);
    parameters.transducers(transducer_ind).source_amp = transducer.source_amp*sqrt(adjustment_ratio);

    % for the output param
    i_axial_oneil = i_axial_oneil*adjustment_ratio;

end

parameters = rmfield(parameters, 'transducer');

end