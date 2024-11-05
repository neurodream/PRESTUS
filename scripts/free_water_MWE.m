clear; close all; clc;

% load the functions of PRESTUS and kWave

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(fileparts(currentFile));

addpath(fullfile(rootpath, 'functions'));
addpath(fullfile(rootpath, 'toolboxes', 'k-wave-toolbox-version-1.4', 'k-Wave'));
addpath(fullfile(rootpath, 'toolboxes', 'FEX-minimize'));

%% params (TODO load them from config to integrate with PRESTUS)

% equipment params (values for 100mm Imasonics transducer)
source_freq_hz = 300e3;
Elements_ID_mm = [10.0, 22.1, 29.8, 36.0, 41.4, 46.3, 50.7, 54.9, 58.7, 62.4];
Elements_OD_mm = [21.1, 28.8, 35.0, 40.4, 45.3, 49.7, 53.9, 57.8, 61.5, 65.0];
curv_radius_mm = 100;
source_phase_rad = deg2rad([0.000, 359.988, 359.999, 14.804, 352.467, 355.567, 0.000, 341.045, 305.927, 0.001]); 
dist_to_exit_plane = 7;
source_amp = 200608;

% design params
expected_focal_distance_mm = 62; % 96.1332
ROI_width_mm = 20;
goal_intensity = 200;
voxels_from_edge = 11;
first_phase_free = false;

% space params
grid_dims = [144, 144, 300];
grid_step_m = 0.5/1000;

% water params
sound_speed = 1500;
density = 994;

% derived params
grid_size = grid_dims*grid_step_m;
axial_position = (1:grid_dims(3))*(grid_step_m*1000);
axial_position = axial_position';
axial_position = axial_position + dist_to_exit_plane;
trans_pos = [72 72 voxels_from_edge];
focus_pos = [72 72 round(voxels_from_edge + expected_focal_distance_mm/(grid_step_m*1000))];
stop_ind = numel(Elements_ID_mm);
if ~first_phase_free
    stop_ind = stop_ind - 1;
end

desired_function = create_boxcar(expected_focal_distance_mm, ROI_width_mm, axial_position, 0, goal_intensity);

%% plot the profile with current set of phases

% run script investigate_intensity_profiles (and adjust phases in radians as desired)

%% optimize phases (active condition)

opt_limits = [axial_position(2,1), axial_position(end,1)];

% create the weights function
[flhm_center, flhm_center_index] = get_flhm_center_position( ...
    axial_position, ...
    desired_function ...
    );
weights = normpdf_nostatstoolbox( ...
    axial_position(:,1), ...
    axial_position(flhm_center_index,1), ...
    axial_position(flhm_center_index,1)/.5 ...
    );

% bring parameters in PRESTUS-readable format
parameters = [];
parameters.transducer.curv_radius_mm = curv_radius_mm;
parameters.transducer.Elements_ID_mm = Elements_ID_mm;
parameters.transducer.Elements_OD_mm = Elements_OD_mm;
parameters.transducer.n_elements     = numel(Elements_ID_mm);
parameters.transducer.source_freq_hz = source_freq_hz;
parameters.medium.water.sound_speed  = sound_speed;
parameters.medium.water.density      = density;

optimize_phases = @(phases_and_velocity) phase_optimization_annulus_full_curve( ...
    phases_and_velocity(1:parameters.transducer.n_elements-1), ...
    parameters, ...
    phases_and_velocity(parameters.transducer.n_elements),...
    axial_position, ...
    desired_function, ...
    0, ...
    opt_limits, ...
    weights);

rng(100,'twister') % setting seed for consistency

velocity = source_amp(1)/(parameters.medium.water.density*parameters.medium.water.sound_speed);   % [m/s]
func = optimize_phases;
x0 = [randi(360, [1 parameters.transducer.n_elements - 1])/180*pi velocity];
lb = zeros(1,parameters.transducer.n_elements);
ub = [2*pi*ones(1,parameters.transducer.n_elements - 1) 0.2];
options = setoptimoptions( ...
    'popsize',1000, ...
    'FinDiffType', 'central', ...
    'MaxFunEvals', 1e6, ...
    'MaxIter', 1e7 ...
    );
[opt_phases_and_velocity, ~, ~, ~] = minimize(func, x0, [],[],[],[],lb, ub, [], options);

% note: I adjusted the phase_optimization_annulus_full_curve function so
% that it returns the axial intensity as the last output parameter
[~, ~, ~, h, i_axial_oneil] = phase_optimization_annulus_full_curve( ...
    opt_phases_and_velocity(1:parameters.transducer.n_elements-1), ...
    parameters, ...
    opt_phases_and_velocity(parameters.transducer.n_elements),...
    axial_position, ...
    desired_function, ...
    1, ...
    opt_limits, ...
    weights ...
    );  

close(h);

% TODO handle if first element should be 0 or not

% postprocessing
degs = rad2deg([0 opt_phases_and_velocity(1:9)]);
degs_str = arrayfun(@(x) sprintf('%.3f', x), degs, 'UniformOutput', false);
fprintf('[%s]\n', strjoin(degs_str, ', '));

disp(['delta focus: ' num2str(axial_position(find(i_axial_oneil == max(i_axial_oneil))) - expected_focal_distance_mm)]);
disp(['delta intensity: ' num2str(max(i_axial_oneil) - goal_intensity)]);

% intermediate plot to check where (if) there should be a cutoff
h = figure;
plot(axial_position, i_axial_oneil);
hold on;
plot(axial_position, desired_function/2);

%% optimize phases (sham condition)

% close(h); % close previous figure

% manual step: set the axial position in mm from which point on the
% intensity should ideally be zero
% TODO maybe just automatically set it to the last minimum before the peak
% closest to intended peak
intensity_cutoff = 51;

desired_function_sham = i_axial_oneil .* single(axial_position < intensity_cutoff); % near field

optimize_phases = @(phases_and_velocity) phase_optimization_annulus_full_curve( ...
    phases_and_velocity(1:parameters.transducer.n_elements-1), ...
    parameters, ...
    phases_and_velocity(parameters.transducer.n_elements),...
    axial_position, ...
    desired_function_sham, ...
    0, ...
    opt_limits, ...
    weights);

velocity = source_amp(1)/(parameters.medium.water.density*parameters.medium.water.sound_speed);   % [m/s]
func = optimize_phases;
x0 = [randi(360, [1 parameters.transducer.n_elements-1])/180*pi velocity];
lb = zeros(1,parameters.transducer.n_elements);
ub = [2*pi*ones(1,parameters.transducer.n_elements-1) 0.2];
options = setoptimoptions( ...
    'popsize',1000, ...
    'FinDiffType', 'central', ...
    'MaxFunEvals', 1e6, ...
    'MaxIter', 1e7 ...
    ); % number of initializations [, 'TolCon', 1e-8]
[opt_phases_and_velocity_sham, min_err, ~, output] = minimize(func, x0, [],[],[],[],lb, ub, [], options);

% note: I adjusted the phase_optimization_annulus_full_curve function so
% that it returns the axial intensity as the last output parameter
[~, ~, ~, h, i_axial_oneil_sham] = phase_optimization_annulus_full_curve( ...
    opt_phases_and_velocity_sham(1:parameters.transducer.n_elements-1), ...
    parameters, ...
    opt_phases_and_velocity_sham(parameters.transducer.n_elements),...
    axial_position, ...
    desired_function_sham, ...
    1, ...
    opt_limits, ...
    weights ...
    );  

close(h);

degs = rad2deg([0 opt_phases_and_velocity_sham(1:9)]);
degs_str = arrayfun(@(x) sprintf('%.3f', x), degs, 'UniformOutput', false);
fprintf('[%s]\n', strjoin(degs_str, ', '));

disp(['delta focus: ' num2str(axial_position(find(i_axial_oneil_sham == max(i_axial_oneil_sham))) - expected_focal_distance_mm)]);
disp(['delta intensity: ' num2str(max(i_axial_oneil_sham) - goal_intensity)]);

%% plot the optimization results

figure('Position', [100, 100, 1000, 400]);
plot(axial_position, i_axial_oneil, 'LineWidth', 2); % NOTE: 0 in axial position should be surface of transducer! therefore needs correcting term
hold on;
plot(axial_position, i_axial_oneil_sham, 'LineWidth', 2);
xticks(linspace(0, max(axial_position), 10));
ylim([0 max(desired_function)]);
xlabel('focal axis (mm)');
ylabel('intensity');
legend('active', 'sham');