%% identify if visual inspection is correct (results: yes, it is correct)

% Assuming focal_planeT is a 512x400x22 double matrix
focal_planeT_CPU = gather(focal_planeT);
[nX, nY, nTime] = size(focal_planeT_CPU);

% Initialize a logical array to track unchanged time steps
unchanged_steps = false(nTime-1, 1);

% Compare each time step with the previous one
for t = 2:nTime
    if isequal(focal_planeT_CPU(:,:,t), focal_planeT_CPU(:,:,t-1))
        unchanged_steps(t-1) = true;
    end
end

% Identify every third time step
every_third_step = 3:3:nTime;

% Check if the values do not change at every third time step
faulty_steps = unchanged_steps(every_third_step - 1);

% Display the results
disp('Indices of unchanged time steps:');
disp(find(unchanged_steps) + 1);

disp('Indices of unchanged every third time steps:');
disp(every_third_step(faulty_steps));


%% run_heating_simulations debug - step 1: load data

% function: 
% run_heating_simulations(sensor_data, kgrid, kwave_medium, sensor, source, parameters, trans_pos)
% i.e. needs: sensor_data, kgrid, kwave_medium, sensor, source, parameters, trans_pos
% might be loaded from 

load('M:\Documents\repos\PRESTUS_forked\data\sims\sbj_3_duty_23_ntrl_4_amp_105000_targ_128-142-160_tpos_209-146-180\sub-003_layered_kwave_source.mat');
% parameters = load_parameters('config_CTX500-024-010_77.0mm.yaml'); % not sure if needed
% still needs: sensor_data, kgrid, kwave_medium - from file below
load('M:\Documents\repos\PRESTUS_forked\data\sims\sbj_3_duty_23_ntrl_4_amp_105000_targ_128-142-160_tpos_209-146-180\sub-003_water_results');

%% run_heating_simulations debug - step 2: call function

trans_pos = [140 140 250];

parameters.data_path = '/project/3015999.02/andche_sandbox/TUS_sims/tusim/data';
parameters.simnibs_bin_path = 'C:/Users/nicade/SimNIBS-4.0/bin';
parameters.ld_library_path = 'H:/sleep/nicade/lib64';
parameters.segmentation_software = 'charm.cmd';

parameters.thermal.n_trials = 4;
parameters.transducer.source_amp = 10*200608;

% Sets up a window to perform simulations in
heating_window_dims(2,3) = parameters.grid_dims(3);
sensor.mask(heating_window_dims(1,1):heating_window_dims(2,1), heating_window_dims(1,2):heating_window_dims(2,2), :) = 1;

% convert to CPU
kwave_medium.density = gather(kwave_medium.density);
kwave_medium.sound_speed = gather(kwave_medium.sound_speed);
kwave_medium.alpha_power = gather(kwave_medium.alpha_power);

% Convert the absorption coefficient to nepers/m
alpha_np = db2neper(kwave_medium.alpha_coeff, kwave_medium.alpha_power) * (2 * pi * parameters.transducer.source_freq_hz).^kwave_medium.alpha_power;

% % Extract the pressure amplitude at each position
% p = extractAmpPhase(sensor_data.p_max_all, 1/kgrid.dt, parameters.transducer.source_freq_hz);
% 
% % Reshape the data, and calculate the volume rate of heat deposition
% p = reshape(p, kgrid.Nx, kgrid.Ny);

p = sensor_data.p_max_all;
p = gather(p);
source.Q = alpha_np .* p.^2 ./ (kwave_medium.density .* kwave_medium.sound_speed); % Heat delivered to the system
source.T0 = parameters.thermal.temp_0; % Initial temperature distribution

% create kWaveDiffusion object
if isfield(parameters.thermal,'record_t_at_every_step') && ~parameters.thermal.record_t_at_every_step 
    sensor = [];
end
thermal_diff_obj = kWaveDiffusion(kgrid, kwave_medium, source, sensor, 'PlotSim', boolean(parameters.interactive));

% New field temperature
% if contains(parameters.code_type, 'gpu') % logical(parameters.using_donders_hpc)
%     thermal_diff_obj.T = gpuArray(thermal_diff_obj.T);
% end
maxT = thermal_diff_obj.T;
maxCEM = thermal_diff_obj.cem43;

[~, on_off_repetitions, on_steps_n,  on_steps_dur, off_steps_n, off_steps_dur, post_stim_steps_n, post_stim_time_step_dur] = check_thermal_parameters(parameters);

time_status_seq = struct('status', {'off'}, 'time', {0}, 'step', {0}, 'recorded', {1});

% Set up some last parameters for simulation
total_timepoints = 1+parameters.thermal.n_trials*(on_off_repetitions*(1+double(off_steps_n>0))+1);
cur_timepoint = 1;
focal_planeT = zeros([size(maxT,[1,3]) total_timepoints]);
focal_planeT(:,:,cur_timepoint) = squeeze(maxT(:,trans_pos(2),:));
focal_planeCEM43 = zeros([size(maxT,[1,3]) total_timepoints]);
% focal_planeCEM43(:,:,cur_timepoint) = squeeze(maxCEM(:,trans_pos(2),:));

% Loop over trials, so that the heat can accumulate
for trial_i = 1:parameters.thermal.n_trials
  fprintf('Trial %i\n', trial_i)
  % Calculate the heat accumulation within a trial
  for pulse_i = 1:on_off_repetitions
      fprintf('Pulse %i\n', pulse_i)
      
      % Pulse on
      thermal_diff_obj.Q = source.Q;
      thermal_diff_obj.takeTimeStep(on_steps_n, on_steps_dur);
      time_status_seq = [time_status_seq struct('time', num2cell(max([time_status_seq(:).time]) + (1:on_steps_n)*on_steps_dur), ...
                                                'step', num2cell(max([time_status_seq(:).step]) + (1:on_steps_n)), 'status', "on",'recorded',0)];
      curT = thermal_diff_obj.T;
      curCEM43 = thermal_diff_obj.cem43;
      cur_timepoint = cur_timepoint+1;
      time_status_seq(end).recorded = 1;
      focal_planeT(:,:,cur_timepoint) = cat(3, squeeze(curT(:,trans_pos(2),:)));
      focal_planeCEM43(:,:,cur_timepoint) = cat(3, squeeze(curCEM43(:,trans_pos(2),:)));
      if any(maxT>curT,'all')
          maxT = curT;
      end
      
      % Pulse off 
      if off_steps_n>0
          thermal_diff_obj.Q = 0;
          thermal_diff_obj.takeTimeStep(off_steps_n, off_steps_dur);
          time_status_seq = [time_status_seq struct('time', num2cell(max([time_status_seq(:).time]) + (1:off_steps_n)*off_steps_dur), 'step', num2cell(max([time_status_seq(:).step]) + (1:off_steps_n)), 'status', "off",'recorded',0)];
          cur_timepoint = cur_timepoint+1;
          time_status_seq(end).recorded = 1;
          curT = thermal_diff_obj.T;
          curCEM43 = thermal_diff_obj.cem43;

          focal_planeT(:,:,cur_timepoint) = cat(3, squeeze(curT(:,trans_pos(2),:)));
          focal_planeCEM43(:,:,cur_timepoint) = cat(3, squeeze(curCEM43(:,trans_pos(2),:)));

      end

  end

end

%% run plotting function
plot_heating_sims(focal_planeT, time_status_seq, parameters, trans_pos, []);