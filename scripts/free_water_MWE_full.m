clear; close all; clc;

% goal output: interactive plots of intensity along focal axis

% load the functions of PRESTUS and kWave

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(fileparts(currentFile));

addpath(fullfile(rootpath, 'functions'));
addpath(fullfile(rootpath, 'toolboxes', 'k-wave-toolbox-version-1.4', 'k-Wave'));
addpath(fullfile(rootpath, 'toolboxes', 'FEX-minimize'));

%% params

% equipment params
source_freq_hz = 300e3;

% values for 100mm transducer:
Elements_ID_mm = [10.0, 22.1, 29.8, 36.0, 41.4, 46.3, 50.7, 54.9, 58.7, 62.4];
Elements_OD_mm = [21.1, 28.8, 35.0, 40.4, 45.3, 49.7, 53.9, 57.8, 61.5, 65.0];
curv_radius_mm = 100; %99.60;

% % values for 75mm transducer:
% Elements_ID_mm = [10.0, 22.3, 30.0, 36.3, 41.7, 46.5, 51.0, 55.1, 58.9, 62.5];
% Elements_OD_mm = [21.3, 29.1, 35.3, 40.7, 45.6, 50.0, 54.1, 58.0, 61.6, 65.0];
% curv_radius_mm = 75;

source_phase_rad = [0,4.76846291859276,4.69990987623293,5.01648467727667,5.49696858436295,0,0,6.28318530717959,6.28318530717959,1.28034037800725]; % zeros(1,10); 
dist_to_exit_plane = 7;

% space params
grid_dims = [144, 144, 256];
grid_step_m = 0.5/1000;
grid_size = grid_dims*grid_step_m;
axial_position = (1:grid_dims(3))*(grid_step_m*1000);
axial_position = axial_position';
axial_position = axial_position + dist_to_exit_plane;

% water params
sound_speed = 1500;
density = 994;

% source_phase_rad = 2 * pi * rand(1, 10);
% source_phase_rad = repmat(1, 1, 10);
source_amp = 200608;%200608;

% design params (center of grid in x and y, beginning of z)
expected_focal_distance_mm = 40;
voxels_from_edge = 11; %11;
trans_pos = [72 72 voxels_from_edge];
focus_pos = [72 72 round(voxels_from_edge + expected_focal_distance_mm/(grid_step_m*1000))];

%% plot the profile with current set of phases

parameters = [];
parameters.transducer.curv_radius_mm = curv_radius_mm;
parameters.transducer.Elements_ID_mm = Elements_ID_mm;
parameters.transducer.Elements_OD_mm = Elements_OD_mm;
parameters.transducer.n_elements = numel(Elements_ID_mm);
parameters.transducer.source_freq_hz = source_freq_hz;
parameters.medium.water.sound_speed = sound_speed;
parameters.medium.water.density = density;

opt_limits = [axial_position(2,1), axial_position(end,1)];
weights = ones(size(axial_position));

[~, ~, ~, h, i_axial_oneil1] = phase_optimization_annulus_full_curve( ...
    [6.27608181712397,0.734280960265038,0.440929510248336,6.18930404671481,0.194866010985167,0.145035813643178,5.91785298815213,0.294489404689003,0.0147502836540872], ...
    parameters, ...
    0.15,...
    axial_position, ...
    axial_position, ...
    1, ...
    opt_limits, ...
    weights ...
    );
close(h);
[~, ~, ~, h, i_axial_oneil2] = phase_optimization_annulus_full_curve( ...
    [4.75078622392856,4.67912300484168,5.00544721508706,5.48415612232406,0,0,6.28318530717959,6.28318530717959,1.28793779624118], ...
    parameters, ...
    0.15,...
    axial_position, ...
    axial_position, ...
    1, ...
    opt_limits, ...
    weights ...
    );
close(h);


h = figure('Position', [100, 100, 800, 200]);

plot(axial_position, i_axial_oneil1);
hold on;
plot(axial_position, i_axial_oneil2);
xlabel('Axial Position [mm]');
ylabel('Intensity [W/cm^2]');
title('Pressure along the beam axis');

%% set up medium

voxels = ones(grid_dims);

% see paper: https://doi.org/10.1121/1.4894790
alpha_coeff = fitPowerLawParamsMulti( ...
    voxels*0, ...
    voxels*1.2, ...
    voxels*sound_speed, ...
    source_freq_hz, ...
    2 ...
    );

kwave_medium = struct( ...
    'sound_speed',          voxels*sound_speed, ...
    'density',              voxels*density, ...
    'alpha_coeff',          alpha_coeff,...
    'alpha_power',          2 , ...
    'thermal_conductivity', voxels*0.6,...
    'specific_heat',        voxels*4178 ...
    );

%% set up grid

% Calculate the time step using an integer number of points per period
points_per_wavelength = sound_speed / (source_freq_hz * grid_step_m);      % points per wavelength
points_per_period = ceil(points_per_wavelength / 0.3);                     % points per period (CFL number is kwave default)
wave_period   = 1 / source_freq_hz;                  % period [s]
grid_time_step = (wave_period / points_per_period)/2;                      % time step [s]

% Calculate the number of time steps to reach steady state
t_end = sqrt(grid_size(1).^2 + grid_size(2).^2 + grid_size(3).^2) / sound_speed; % [s]
simulation_time_points = round(t_end / grid_time_step);

kgrid = kWaveGrid( ...
    grid_dims(1), grid_step_m, ...
    grid_dims(2), grid_step_m, ...
    grid_dims(3), grid_step_m);
kgrid.setTime(simulation_time_points, grid_time_step);

%% set up sensor
% Creates a sensor that records the the maximum and final pressure
% values in every point of the grid
sensor = struct();
sensor.mask = ones(grid_dims);
sensor.record = {'p_max_all','p_final'};

% Record the last 3 cycles in steady state
num_periods = 3;
time_points_to_record = round(num_periods * wave_period / kgrid.dt);
sensor.record_start_index = simulation_time_points - time_points_to_record + 1;

%% set up source
focus_pos = focus_pos';
trans_pos = trans_pos';

Elements_OD = 2*floor(Elements_OD_mm / (grid_step_m*1000) / 2) + 1; % [grid points]
Elements_ID = 2*floor(Elements_ID_mm / (grid_step_m*1000) / 2) + 1; % [grid points]
Elements_ID(Elements_ID_mm==0) = 0;
radius_grid = round(curv_radius_mm / (grid_step_m*1000));  % [grid points]

cw_signal = createCWSignals( ...
    kgrid.t_array, ...
    source_freq_hz, ...
    repmat(source_amp, 1, numel(Elements_OD_mm)), ...
    source_phase_rad ...
    );

source = struct();

% create empty kWaveArray
karray = kWaveArray('BLITolerance', 0.1, 'UpsamplingRate', 10, 'BLIType', 'sinc');

% add bowl shaped element
karray.addAnnularArray( ...
    [kgrid.x_vec(trans_pos(1)) kgrid.y_vec(trans_pos(2)) kgrid.z_vec(trans_pos(3))], ...
    curv_radius_mm*1e-3, ...
    [Elements_ID_mm; Elements_OD_mm]*1e-3, ...
    [kgrid.x_vec(focus_pos(1)) kgrid.y_vec(focus_pos(2)) kgrid.z_vec(focus_pos(3))] ...
    );

% somehow compute the weights of sources adjusted for grid coarseness (?); 
% see https://doi.org/10.1121/1.5116132
grid_weights_4d = zeros(numel(Elements_OD_mm), kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1));
for ind = 1:numel(Elements_OD_mm)
    fprintf('Computing weights for element %i...', ind); 
    grid_weights_4d(ind,:,:,:) = karray.getElementGridWeights(kgrid, ind);                
    fprintf(' done\n');
end

% get the binary mask
binary_mask = squeeze(sum(grid_weights_4d, 1)) ~= 0;
mask_ind = find(binary_mask);
num_source_points = sum(binary_mask(:));

% number of time points in the signal
Nt = size(cw_signal, 2);
 
% bits of code from kWaveArray.m
% initialise the source signal
distributed_source_signal = zeros(num_source_points, Nt);

% source_labels = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz); % max(kgrid.Nz, 1)

% loop through the elements
for ind = 1:numel(Elements_OD_mm)

    % get the offgrid source weights
    source_weights = squeeze(grid_weights_4d(ind,:,:,:));
    el_binary_mask = source_weights ~= 0;
    % get indices of the non-zero points 
    element_mask_ind = find(el_binary_mask);
    % convert these to indices in the distributed source
    local_ind = ismember(mask_ind, element_mask_ind);

    % add to distributed source
    distributed_source_signal(local_ind, :) = distributed_source_signal(local_ind, :) + source_weights(element_mask_ind) * cw_signal(ind, :);

    % source_labels = source_labels + ind * el_binary_mask;

end
        
% assign binary mask
source.p_mask = binary_mask;

% assign source signals
source.p = distributed_source_signal;


%% RUN ACOUSTIC SIMULATION

% Defines the edge of the simulation as the edge of the PML layer (see line 148)
kwave_input_args = struct( ...
    'PMLInside', true, ...
    'PMLSize', 10, ...
    'PlotPML', true ...d
    );

input_args.PlotSim = false;
input_args.PlotScale = [-1, 1] * source_amp;
input_args.DataCast = 'gpuArray-single';

input_args_cell = zip_fields(input_args);

try
    kwave_medium = rmfield(kwave_medium,'thermal_conductivity');
    kwave_medium = rmfield(kwave_medium,'specific_heat');
catch
end

sensor_data = kspaceFirstOrder3D(kgrid, kwave_medium, source, sensor, input_args_cell{:});

% post processing
p_final = gather(sensor_data.p_final);
p_max = gather(sensor_data.p_max_all);
intensity = p_max .^2/(2*sound_speed*density) .* 1e-4;
intensity_1D = squeeze(intensity(72,72,:));

[flhm_center, flhm_center_index] = get_flhm_center_position( ...
    axial_position, ...
    intensity_1D ...
    );
weights = normpdf_nostatstoolbox( ...
    axial_position(:,1), ...
    axial_position(flhm_center_index,1), ...
    axial_position(flhm_center_index,1)/.5 ...
    );

%% plotting

% Create a new figure for subplots
h = figure('Position', [10 10 1400 500]);

% Adjust the space between subplots
subplot(1, 2, 1);
imagesc((1:size(p_max,1))*(grid_step_m*1000), ...
    (1:size(p_max,3))*(grid_step_m*1000), ...
    squeeze(p_max(:,trans_pos(2),:))');
axis image;
colormap(getColorMap);
xlabel('Lateral Position [mm]');
ylabel('Axial Position [mm]');
colorbar;
title('Pressure for the focal plane');

% Second subplot: Focal axis pressure
subplot('Position', [0.51, 0.1, 0.48, 0.8]); % Adjusted position for the longer plot
plot(axial_position, intensity_1D);
xlabel('Axial Position [mm]');
ylabel('Intensity [W/cm^2]');
title('Pressure along the beam axis');
hold on;

plot(axial_position, weights/max(weights) - min(weights)/max(weights));

%% optimize phases
% TODO not quite clear why the last element is dropped and replaced with
% velocity (also, what is velocity in this case?)

opt_limits = [axial_position(2,1), axial_position(end,1)];

i_axial_oneil(axial_position > 42.5) = 0; desired_function = i_axial_oneil;
% desired_function = create_boxcar(67.4, 20, axial_position, 0, 200); % 850

% recreate the weights function
[flhm_center, flhm_center_index] = get_flhm_center_position( ...
    axial_position, ...
    desired_function ...sign
    );
weights = normpdf_nostatstoolbox( ...
    axial_position(:,1), ...
    axial_position(flhm_center_index,1), ...
    axial_position(flhm_center_index,1)/.5 ...
    );
% weights = ones(size(axial_position));

parameters = [];
parameters.transducer.curv_radius_mm = curv_radius_mm;
parameters.transducer.Elements_ID_mm = Elements_ID_mm;
parameters.transducer.Elements_OD_mm = Elements_OD_mm;
parameters.transducer.n_elements = numel(Elements_ID_mm);
parameters.transducer.source_freq_hz = source_freq_hz;
parameters.medium.water.sound_speed = sound_speed;
parameters.medium.water.density = density;

cd(rootpath);
% load('scripts/tests/parameter_comparison_MWE.mat');
% load('scripts/tests/parameter_comparison_MWE_2.mat');
% parameters = param_cmp.parameters;

% % TODO remove debug movie
% video = VideoWriter('real_time_animation.avi');
% video.FrameRate = 60; % Set the frame rate
% open(video); % Open the video file for writing
% f = figure('Position', [100, 100, 1200, 400]);
% 
% optimize_phases = @(phases_and_velocity) phase_optimization_annulus_full_curve_video( ...
%     phases_and_velocity(1:parameters.transducer.n_elements-1), ...
%     parameters, ...
%     phases_and_velocity(parameters.transducer.n_elements),...
%     axial_position, ...
%     desired_function, ...
%     1, ...
%     opt_limits, ...
%     weights, ...
%     video, ...
%     f);


% optimize_phases = @(phases_and_velocity) phase_optimization_annulus_full_curve_all_elements( ...
%     phases_and_velocity(1:parameters.transducer.n_elements), ...
%     parameters, ...
%     phases_and_velocity(parameters.transducer.n_elements + 1),...
%     axial_position, ...
%     desired_function, ...
%     0, ...
%     opt_limits, ...
%     weights);

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
% x0 = [randi(360, [1 parameters.transducer.n_elements])/180*pi velocity];
% lb = zeros(1,parameters.transducer.n_elements + 1);
% ub = [2*pi*ones(1,parameters.transducer.n_elements) 0.2];
x0 = [randi(360, [1 parameters.transducer.n_elements - 1])/180*pi velocity];
lb = zeros(1,parameters.transducer.n_elements);
ub = [2*pi*ones(1,parameters.transducer.n_elements - 1) 0.2];
options = setoptimoptions( ...
    'popsize',1000, ...
    'FinDiffType', 'central', ...
    'MaxFunEvals', 1e6, ...
    'MaxIter', 1e7 ...
    ); % number of initializations [, 'TolCon', 1e-8]
[opt_phases_and_velocity, min_err, ~, output] = minimize(func, x0, [],[],[],[],lb, ub, [], options);

% close(video);

opt_velocity = opt_phases_and_velocity(end);
% opt_velocity = velocity;
opt_phases = [0 opt_phases_and_velocity(1:numel(source_phase_rad))];
% opt_phases = opt_phases_and_velocity;

% plot optimization results

% [~, ~, ~, h, i_axial_oneil] = phase_optimization_annulus_full_curve_all_elements( ...
%     opt_phases_and_velocity(1:parameters.transducer.n_elements), ...
%     parameters, ...
%     opt_phases_and_velocity(parameters.transducer.n_elements + 1),...
%     axial_position, ...
%     desired_function, ...
%     1, ...
%     opt_limits, ...
%     weights ...
%     );  

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

% +(voxels_from_edge*grid_step_m*1000)
figure('Position', [100, 100, 1000, 400]);
plot(axial_position, i_axial_oneil); % NOTE: 0 in axial position should be surface of transducer! therefore needs correcting term
hold on;
plot(axial_position, desired_function);
plot(axial_position, weights/max(weights) - min(weights)/max(weights));
xticks(linspace(0, max(axial_position), 10));
ylim([0 max(desired_function)]);

% % store phases
% save('../../scans/phases_lowerint', 'opt_phases_and_velocity');

% % try how it looks with double the intensity...
% 
% % plot optimization results
% [~, ~, ~, h, i_axial_oneil_double] = phase_optimization_annulus_full_curve_all_elements( ...
%     opt_phases_and_velocity(1:parameters.transducer.n_elements), ...
%     parameters, ...
%     opt_phases_and_velocity(parameters.transducer.n_elements + 1) * 2,...
%     axial_position, ...
%     desired_function, ...
%     1, ...
%     opt_limits, ...
%     weights ...
%     );  
% 
% [~, ~, ~, h, i_axial_oneil_double] = phase_optimization_annulus_full_curve( ...
%     opt_phases_and_velocity(1:parameters.transducer.n_elements-1), ...
%     parameters, ...
%     opt_phases_and_velocity(parameters.transducer.n_elements + 1) * 2,...
%     axial_position, ...
%     desired_function, ...
%     1, ...
%     opt_limits, ...
%     weights ...
%     );  
% 
% close(h);
% 
% plot(axial_position, i_axial_oneil_double); % NOTE: 0 in axial position should be surface of transducer! therefore needs correcting term
% 
% ylim([0 max([max(desired_function) max(i_axial_oneil_double)])]);

%% store the near field

desired_function = i_axial_oneil .* single(axial_position<49); % near field

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
x0 = [randi(360, [1 parameters.transducer.n_elements-1])/180*pi velocity];
lb = zeros(1,parameters.transducer.n_elements);
ub = [2*pi*ones(1,parameters.transducer.n_elements-1) 0.2];
options = setoptimoptions( ...
    'popsize',1000, ...
    'FinDiffType', 'central', ...
    'MaxFunEvals', 1e6, ...
    'MaxIter', 1e7 ...
    ); % number of initializations [, 'TolCon', 1e-8]
[opt_phases_and_velocity, min_err, ~, output] = minimize(func, x0, [],[],[],[],lb, ub, [], options);

% save('../../scans/phases_near_field_lowerint', 'opt_phases_and_velocity');

%% test if phases are sufficient for replication

opt_phases_target = load('../../scans/phases_lowerint');
opt_phases_target = opt_phases_target.opt_phases_and_velocity;
opt_phases_control = load('../../scans/phases_near_field_lowerint');
opt_phases_control = opt_phases_control.opt_phases_and_velocity;

[~, ~, ~, h, i_axial_oneil_target] = phase_optimization_annulus_full_curve( ...
    opt_phases_target(1:parameters.transducer.n_elements-1), ...
    parameters, ...
    opt_phases_target(parameters.transducer.n_elements),...
    axial_position, ...
    desired_function, ...
    1, ...
    opt_limits, ...
    weights ...
    );
close(h);

[~, ~, ~, h, i_axial_oneil_control] = phase_optimization_annulus_full_curve( ...
    opt_phases_control(1:parameters.transducer.n_elements-1), ...
    parameters, ...
    opt_phases_control(parameters.transducer.n_elements),...
    axial_position, ...
    desired_function, ...
    1, ...
    opt_limits, ...
    weights ...
    );
close(h);

% +(voxels_from_edge*grid_step_m*1000)
figure('Position', [100, 100, 1000, 400]);
plot(axial_position, i_axial_oneil_target, 'LineWidth', 2); % NOTE: 0 in axial position should be surface of transducer! therefore needs correcting term
hold on;
plot(axial_position, i_axial_oneil_control, 'LineWidth', 2);
% plot(axial_position, weights/max(weights) - min(weights)/max(weights));
xlabel('focal axis');
ylabel('free water intensity');
legend('active', 'control');
xticks(linspace(0, max(axial_position), 10));
ylim([0 1.3*max(i_axial_oneil_target)]);

% near_field = i_axial_oneil .* single(axial_position<49);
% 
% figure;
% plot(axial_position, i_axial_oneil);
% hold on;
% plot(axial_position, i_axial_oneil);