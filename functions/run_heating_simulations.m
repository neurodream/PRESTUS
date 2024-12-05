function [thermal_diff_obj, time_status_seq, maxT,focal_planeT, maxCEM43, CEM43, tissue_specific_heat, tissue_specific_CEM43] = run_heating_simulations(sensor_data, kgrid, kwave_medium, sensor, source, parameters, trans_pos)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                    Runs the heating simulations                   %
%                                                                   %
% Instead of modelling each oscillation for measuring heating       %
% effects, k-wave divides each duty cycle up into segments called   %
% 'sim_time_steps' with the total of one duty cycle being           %
% represented as the 'on_off_step_duration'.                        %
% The stable state of the acoustic simulations (what is seen in the %
% figures and nifti files) is used as the input for the temperature %
% simulations.                                                      %
%                                                                   %
% For a detailed explanation on how to correctly configure your     %
% thermal parameters, see thermal simulations getting started.      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



% Convert the absorption coefficients to nepers/m (!)
% for the following, see also fitPowerLawParamsMulti.m
% define frequency in rad/s
w = 2*pi*parameters.transducers(1).source_freq_hz; % TODO note assumes same source freq for every transducer!! enable different source freqs between transducers
% convert absorption to Nepers/((rad/s)^y m)
a0_np = db2neper(kwave_medium.alpha_coeff, kwave_medium.alpha_power);
alpha_np = a0_np.*w.^kwave_medium.alpha_power;
clear w a0_np;

% alternative simplified conversion dB to Nepers
% alpha_np = (100 * kwave_medium.alpha_coeff .* (parameters.transducer.source_freq_hz/10^6)^kwave_medium.alpha_power)/8.686;

% Get the maximum pressure (in Pa) and calculate Q, the volume rate of heat deposition
p = gather(abs(sensor_data.p_max_all));
source.Q = (alpha_np .* p.^2) ./ (kwave_medium.density .* kwave_medium.sound_speed); % Heat delivered to the system (W/m3)
source.T0 = kwave_medium.temp_0; %parameters.thermal.temp_0; % Initial temperature distribution
clear alpha_np p;

% split temp_0 off from kWave_medium
kwave_medium = rmfield(kwave_medium, 'temp_0');

% create kWaveDiffusion object
if isfield(parameters.thermal,'record_t_at_every_step') && ~parameters.thermal.record_t_at_every_step 
    sensor = [];
end
thermal_diff_obj = kWaveDiffusion(kgrid, kwave_medium, source, sensor, 'PlotSim', boolean(parameters.interactive));

% New field temperature
thermal_diff_obj.T = gpuArray(thermal_diff_obj.T);
maxT = thermal_diff_obj.T;

% new field for cem43 (replicate above)
thermal_diff_obj.cem43 = gpuArray(zeros(size(thermal_diff_obj.T)));
maxCEM43 = thermal_diff_obj.cem43;

[~, on_off_repetitions, on_steps_n,  on_steps_dur, off_steps_n, off_steps_dur, post_stim_steps_n, post_stim_time_step_dur] = ...
    check_thermal_parameters(parameters);







% Calculate total time steps during trials
if off_steps_n > 0
    steps_per_pulse = on_steps_n + off_steps_n;
else
    steps_per_pulse = on_steps_n;
end

steps_per_trial = on_off_repetitions * steps_per_pulse;
total_steps_during_trials = parameters.thermal.n_trials * steps_per_trial;

% Calculate total trials time
total_trials_time = parameters.thermal.n_trials * on_off_repetitions * ...
    (on_steps_n * on_steps_dur + off_steps_n * off_steps_dur);

% Set cooling step duration and number of steps
cooling_step_duration = on_steps_dur + off_steps_dur; % or choose a suitable value
cooling_steps_n = ceil(total_trials_time / cooling_step_duration);

% Adjust cooling_steps_n if necessary
if cooling_steps_n ~= round(cooling_steps_n)
    cooling_steps_n = round(cooling_steps_n);
    cooling_step_duration = total_trials_time / cooling_steps_n;
    warning('Cooling steps adjusted to integer number of steps.');
end

% Calculate total time points
total_timepoints = 1 + total_steps_during_trials + post_stim_steps_n + cooling_steps_n;

% Initialize focal_planeT and CEM43 with the total number of time points
focal_planeT = gpuArray(NaN([size(maxT,[1,3]) total_timepoints]));
focal_planeT(:,:,1) = squeeze(maxT(:,trans_pos(2),:));

CEM43 = gpuArray(NaN([size(maxT,[1,3]) total_timepoints]));
CEM43(:,:,1) = squeeze(maxCEM43(:,trans_pos(2),:));

cur_timepoint = 1;

% Initialize time_status_seq
time_status_seq = struct('status', {'off'}, 'time', {0}, 'step', {0}, 'recorded', {1});

% initialize an output that for every tissue and time point shows the
% maximum value
tissues = fieldnames(parameters.medium);
tissue_specific_heat  = zeros(numel(tissues), total_timepoints);
tissue_specific_CEM43 = zeros(numel(tissues), total_timepoints);

% initialize a logical 3D matrix the shape of the whole medium
tissue_labels = zeros(size(kwave_medium.sound_speed));
for tissue_ID = 1:length(tissues)
    tissue = tissues{tissue_ID};
    tissue_soundspeed = parameters.medium.(tissue).sound_speed;
    tissue_labels(kwave_medium.sound_speed == tissue_soundspeed) = tissue_ID;
end

% % old approach to have cool-off period
% time_status_seq = struct('status', {'off'}, 'time', {0}, 'step', {0}, 'recorded', {1});
% 
% % Set final parameters for simulation
% total_timepoints = 1+parameters.thermal.n_trials*(on_off_repetitions*(1+double(off_steps_n>0))+1);
% total_timepoints = total_timepoints*2;
% cur_timepoint = 1;
% focal_planeT = gpuArray(NaN([size(maxT,[1,3]) total_timepoints]));
% focal_planeT(:,:,cur_timepoint) = squeeze(maxT(:,trans_pos(2),:));
% 
% CEM43 = gpuArray(NaN([size(maxT,[1,3]),total_timepoints]));
% CEM43(:,:,cur_timepoint) = squeeze(maxCEM43(:,trans_pos(2),:));



% Loop over trials, so that the heat can accumulate
for trial_i = 1:parameters.thermal.n_trials
  fprintf('Trial %i\n', trial_i)

    is_break_period = 0;
    if isfield(parameters, 'start_break_trials')
        for i_break = 1:length(parameters.start_break_trials)
            if trial_i >= parameters.start_break_trials(i_break) && trial_i <= parameters.stop_break_trials(i_break)
                is_break_period = 1;
            end
        end
    end
    
  % Calculate the heat accumulation within a trial
  for pulse_i = 1:on_off_repetitions
      fprintf('Pulse %i\n', pulse_i)

      if is_break_period == 0
          thermal_diff_obj.Q = source.Q;
          tmp_status = 'on';
      else
          thermal_diff_obj.Q(:,:,:) = 0;
          tmp_status = 'off';
          disp('Break active');
      end

      thermal_diff_obj.takeTimeStep(on_steps_n, on_steps_dur);
      
      time_status_seq = [time_status_seq, ...
          struct('time', num2cell(max([time_status_seq(:).time]) + (1:on_steps_n)*on_steps_dur), ...
                 'step', num2cell(max([time_status_seq(:).step]) + (1:on_steps_n)), ...
                 'status', tmp_status,...
                 'recorded',1)];
      curT = thermal_diff_obj.T;
      curCEM43 = thermal_diff_obj.cem43;

      cur_timepoint = cur_timepoint+1;
      focal_planeT(:,:,cur_timepoint) = squeeze(curT(:,trans_pos(2),:));
      CEM43(:,:,cur_timepoint) = squeeze(curCEM43(:,trans_pos(2),:));
      % fill in the maximum tissue thermal properties
      for tissue_ID = 1:length(tissues)
          cur_max_T = max(curT(tissue_labels == tissue_ID));
          cur_max_CEM43 = max(curCEM43(tissue_labels == tissue_ID));
          if ~isempty(cur_max_T)
              tissue_specific_heat(tissue_ID, cur_timepoint) = cur_max_T;
          end
          if ~isempty(cur_max_CEM43)
              tissue_specific_CEM43(tissue_ID, cur_timepoint) = cur_max_CEM43;
          end
      end

      % update maxT and maxCEM43 where applicable
      maxT = max(maxT, curT);
      maxCEM43 = max(maxCEM43, curCEM43);

      % Pulse off 
      if off_steps_n>0
          thermal_diff_obj.Q(:,:,:) = 0;
          thermal_diff_obj.takeTimeStep(off_steps_n, off_steps_dur);
          time_status_seq = [time_status_seq,...
              struct('time', num2cell(max([time_status_seq(:).time]) + (1:off_steps_n)*off_steps_dur), ...
              'step', num2cell(max([time_status_seq(:).step]) + (1:off_steps_n)),...
              'status', "off",...
              'recorded',1)];
          cur_timepoint = cur_timepoint+1;
          curT = thermal_diff_obj.T;
          curCEM43 = thermal_diff_obj.cem43;
          focal_planeT(:,:,cur_timepoint) = ...
              squeeze(curT(:,trans_pos(2),:));
          CEM43(:,:,cur_timepoint) = ...
              squeeze(curCEM43(:,trans_pos(2),:));
          % fill in the maximum tissue thermal properties
          for tissue_ID = 1:length(tissues)
              cur_max_T = max(curT(tissue_labels == tissue_ID));
              cur_max_CEM43 = max(curCEM43(tissue_labels == tissue_ID));
              if ~isempty(cur_max_T)
                  tissue_specific_heat(tissue_ID, cur_timepoint) = cur_max_T;
              end
              if ~isempty(cur_max_CEM43)
                  tissue_specific_CEM43(tissue_ID, cur_timepoint) = cur_max_CEM43;
              end
          end
      end
  end
end

% Post-stimulation period
if post_stim_steps_n > 0 
  thermal_diff_obj.Q(:,:,:) = 0;
  thermal_diff_obj.takeTimeStep(post_stim_steps_n, post_stim_time_step_dur);
  time_status_seq = [time_status_seq struct(...
      'time', num2cell(max([time_status_seq(:).time]) + (1:off_steps_n)*off_steps_dur), ...
      'step', num2cell(max([time_status_seq(:).step]) + (1:off_steps_n)), ...
      'status', "off",...
      'recorded',1)];
  cur_timepoint = cur_timepoint+1;
  curT = thermal_diff_obj.T;
  curCEM43 = thermal_diff_obj.cem43;
  focal_planeT(:,:,cur_timepoint) = squeeze(curT(:,trans_pos(2),:));
  CEM43(:,:,cur_timepoint) = squeeze(curCEM43(:,trans_pos(2),:));
  % fill in the maximum tissue thermal properties
  for tissue_ID = 1:length(tissues)
      cur_max_T = max(curT(tissue_labels == tissue_ID));
      cur_max_CEM43 = max(curCEM43(tissue_labels == tissue_ID));
      if ~isempty(cur_max_T)
        tissue_specific_heat(tissue_ID, cur_timepoint) = cur_max_T;
      end
      if ~isempty(cur_max_CEM43)
          tissue_specific_CEM43(tissue_ID, cur_timepoint) = cur_max_CEM43;
      end
  end
end

% Cool-off period
fprintf('Starting cool-off period...\n');

thermal_diff_obj.Q(:,:,:) = 0; % Turn off the sound source

for step_i = 1:cooling_steps_n
    thermal_diff_obj.takeTimeStep(1, cooling_step_duration);
    
    % Update time_status_seq
    time_status_seq = [time_status_seq, ...
        struct('time', max([time_status_seq(:).time]) + cooling_step_duration, ...
               'step', max([time_status_seq(:).step]) + 1, ...
               'status', "off", ...
               'recorded', 1)];
    
    % Update temperatures and CEM43
    cur_timepoint = cur_timepoint + 1;
    if cur_timepoint > total_timepoints
        error('cur_timepoint exceeds total_timepoints during cool-off period');
    end
    curT = thermal_diff_obj.T;
    curCEM43 = thermal_diff_obj.cem43;
    focal_planeT(:,:,cur_timepoint) = squeeze(curT(:,trans_pos(2),:));
    CEM43(:,:,cur_timepoint) = squeeze(curCEM43(:,trans_pos(2),:));
    % fill in the maximum tissue thermal properties
    for tissue_ID = 1:length(tissues)
      cur_max_T = max(curT(tissue_labels == tissue_ID));
      cur_max_CEM43 = max(curCEM43(tissue_labels == tissue_ID));
      if ~isempty(cur_max_T)
        tissue_specific_heat(tissue_ID, cur_timepoint) = cur_max_T;
      end
      if ~isempty(cur_max_CEM43)
          tissue_specific_CEM43(tissue_ID, cur_timepoint) = cur_max_CEM43;
      end
    end
    
    % Update maxT and maxCEM43 where applicable
    maxT = max(maxT, curT);
    maxCEM43 = max(maxCEM43, curCEM43);
end

end