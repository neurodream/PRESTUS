function [parameters, axial_position, i_axial_oneil] = calculate_transducer_phases(parameters, transducer_ind, expected_focal_distance_mm, ROI_width_mm, goal_intensity, sham, desired_function, use_all_phases)

% TODO adjust the velocity to match the intended strength
% TODO replace ROI_width_mm with config file
% TODO maybe just return transducer parameters

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

parameters.transducers(transducer_ind).source_amp = repmat(opt_velocity*parameters.medium.water.sound_speed*parameters.medium.water.density, 1, 10); % TODO remove the hardcoding of element numbers


% % TODO run this only when validated phases should be used;
% % should be read out if desired 
% 
%             if transducer_ind == 1
%                 parameters.transducers(1).source_phase_rad = deg2rad([163.85, 148.61, 133.58, 117.97, 101.91, 85.67, 68.80, 51.59, 34.31, 16.72]);
%                 parameters.transducers(1).source_amp = repmat(350000, 1, 10);
%             else
%                 parameters.transducers(2).source_phase_rad = deg2rad([232.81, 230.72, 228.65, 226.51, 224.31, 222.10, 219.80, 217.46, 215.11, 212.73]);
%                 parameters.transducers(2).source_amp = repmat(250000, 1, 10);
%             end
% 
%             parameters = get_simulated_axial_intensity(parameters);
% 
%             transducer = parameters.transducers(transducer_ind);
% 
%             i_axial_oneil = transducer.axial_intensity_sim_FW;
% 
%             lb = zeros(1,transducer.n_elements);
%             if use_all_phases
%                 stop_before = 0;
%                 multiply_velocity = 0;
%                 add_velocity = 0.15;
%                 optimization_function = @phase_optimization_annulus_full_curve_all_elements;
%                 ub = 2*pi*ones(1,transducer.n_elements);
%                 x0 = randi(360, [1 transducer.n_elements])/180*pi;
%             else
%                 stop_before = 1;
%                 multiply_velocity = 1;
%                 add_velocity = 0;
%                 optimization_function = @phase_optimization_annulus_full_curve;
%                 ub = [2*pi*ones(1,transducer.n_elements - 1) 0.2];
%                 x0 = [randi(360, [1 transducer.n_elements - 1])/180*pi 0.15];
%             end
%             options = setoptimoptions( ...
%                 'popsize',1000, ...
%                 'FinDiffType', 'central', ...
%                 'MaxFunEvals', 1e7, ...
%                 'MaxIter', 1e6 ...
%                 );


% adjust intensity levels to desired intensity

% % TODO remove debug loop (get some data if scaling effect of source_amp is linear)
% figure; hold on;
% for scale_factor = 1.82:0.001:1.83
%     source_amp = transducer.source_amp;
%     source_amp = source_amp*scale_factor;
%     p_axial_oneil = focusedAnnulusONeil( ...
%         transducer.curv_radius_mm/1e3, ...
%         [transducer.Elements_ID_mm; transducer.Elements_OD_mm]/1e3, ...
%         (parameters.transducers(transducer_ind).source_amp*scale_factor)/(parameters.medium.water.density*parameters.medium.water.sound_speed),...%repmat(opt_velocity*scale_factor, 1, 10), ...
%         opt_phases, ...
%         transducer.source_freq_hz, ...
%         parameters.medium.water.sound_speed, ...
%         parameters.medium.water.density, ...
%         (axial_position-0.5)*1e-3 ... % before here was ax_pos
%         );
%     i_axial_oneil = p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;
%     plot(i_axial_oneil);
% end
% % legend([0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2]);
% hold off;

% % new version
% 
% % define the desired intensity in pressure
% goal_pressure = sqrt(goal_intensity * 2 * parameters.medium.water.sound_speed * parameters.medium.water.density * 1e4);
% [peaks_list, peaks_inds] = findpeaks(p_axial_oneil); 
% peak = peaks_list(end);
% % attempt to robustly exclude "fake peaks" in the far field with super low intensity (TODO check if really robust)
% if numel(peaks_list) > 1
%     if peaks_list(end-1) > peaks_list(end)*10
%         peak = peaks_list(end-1);
%     end
% end
% 
% adjustment_ratio = goal_pressure/peak;
% 
% parameters.transducers(transducer_ind).source_amp = parameters.transducers(transducer_ind).source_amp*adjustment_ratio;
% 
% % check if correctly adjusted, and adjust intensity for sham to work correctly:
% p_axial_oneil = focusedAnnulusONeil( ...
%     transducer.curv_radius_mm/1e3, ...
%     [transducer.Elements_ID_mm; transducer.Elements_OD_mm]/1e3, ...
%     parameters.transducers(transducer_ind).source_amp/(parameters.medium.water.density*parameters.medium.water.sound_speed), ...%repmat(opt_velocity*adjustment_ratio, 1, 10), ...
%     opt_phases, ...
%     transducer.source_freq_hz, ...
%     parameters.medium.water.sound_speed, ...
%     parameters.medium.water.density, ...
%     (axial_position-0.5)*1e-3 ... % before here was ax_pos
%     );
% i_axial_oneil = p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;
% 
% disp('')

% old version

[peaks_list, peaks_inds] = findpeaks(i_axial_oneil); 
peak = peaks_list(end);
% attempt to robustly exclude "fake peaks" in the far field with super low
% intensity (TODO check if really robust)
if numel(peaks_list) > 1
    if peaks_list(end-1) > peaks_list(end)*4
        peak = peaks_list(end-1);
    end
end
adjustment_ratio = goal_intensity/peak; %max([goal_intensity peak])/min([goal_intensity peak]);
parameters.transducers(transducer_ind).source_amp = parameters.transducers(transducer_ind).source_amp*sqrt(adjustment_ratio); % TODO check if incorrectly adjusted!

% also adjust the intensity so that the sham works correctly
i_axial_oneil = i_axial_oneil*adjustment_ratio;
i_axial_oneil_active = i_axial_oneil; % keep for reference for debugging
desired_profile_active = desired_profile;

% % debug
% figure; hold on;
% plot(axial_position, i_axial_oneil_active);
% plot(axial_position, desired_profile);
% hold off;

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

    parameters.transducers(transducer_ind).source_amp = repmat(opt_velocity*parameters.medium.water.sound_speed*parameters.medium.water.density, 1, 10); % TODO remove the hardcoding of element numbers

    % % adjust intensity levels to desired intensity
    % sham_peak = max(i_axial_oneil);
    % adjustment_ratio = max([peak_desired_sham_profile sham_peak])/min([peak_desired_sham_profile sham_peak]);
    % parameters.transducers(transducer_ind).source_amp = transducer.source_amp*sqrt(adjustment_ratio);
    % 
    % % for the output param
    % i_axial_oneil = i_axial_oneil*adjustment_ratio;

end

disp(transducer_ind)
disp(sham)
disp(rad2deg(opt_phases));
disp(opt_velocity*parameters.medium.water.sound_speed*parameters.medium.water.density);

% % debug
% figure; hold on;
% plot(axial_position, i_axial_oneil_active);
% plot(axial_position, i_axial_oneil);
% plot(axial_position, desired_profile_active);
% hold off;

parameters = rmfield(parameters, 'transducer');

end