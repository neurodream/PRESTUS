function [kgrid, source, sensor, source_labels] = setup_grid_source_sensor(parameters, max_sound_speed, trans_pos_final, focus_pos_final)
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                  Set up the transducer and sensor                 %
    %                                                                   %
    % This function sets up the transducer in the grid, the timeperiod  %
    % during which simulations will take place and the sensor that      %
    % records the pressure-levels in the grid.                          %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % TODO trans_pos_final, focus_pos_final not referenced anymore

    % Creates a simulation grid in 3 or 2 dimensions
    if parameters.n_sim_dims == 3
        kgrid = kWaveGrid(parameters.grid_dims(1), parameters.grid_step_m, ...
                      parameters.grid_dims(2), parameters.grid_step_m, ...
                      parameters.grid_dims(3), parameters.grid_step_m);
    else
        kgrid = kWaveGrid(parameters.grid_dims(1), parameters.grid_step_m, ...
                      parameters.grid_dims(2), parameters.grid_step_m);
    end
    
    % TODO figure out how to handle multiple source frequencies
    source_freq_hz = parameters.transducers(1).source_freq_hz;
    
    wave_period   = 1 / source_freq_hz;                                                       % period [s]
    
    % Check the number of input arguments
    if nargin < 5
        % Calculate the time step using an integer number of points per period
        points_per_wavelength = max_sound_speed / (source_freq_hz * parameters.grid_step_m);      % points per wavelength
        cfl = 0.3;                                                                                                      % CFL number (kwave default)
        points_per_period = ceil(points_per_wavelength / cfl);                                                          % points per period
        grid_time_step = (wave_period / points_per_period)/2;                                                           % time step [s]  
     end

    % Calculate the number of time steps to reach steady state
    t_end = sqrt(kgrid.x_size.^2 + kgrid.z_size.^2 + kgrid.y_size.^2) / max_sound_speed;                           % [s]
    simulation_time_points = round(t_end / grid_time_step);

    % Create the time array
    kgrid.setTime(simulation_time_points, grid_time_step);

    % Create source (transducer with focuspoint)
    % TODO notice that changed - make it a parameter if source file should
    % be stored
    parameters.kwave_source_filename  = fullfile(parameters.output_dir, sprintf('sub-%03d_%s_kwave_source%s.mat', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    if confirm_overwriting(parameters.kwave_source_filename, parameters)
        [source, source_labels, ~] = setup_source(parameters, kgrid, trans_pos_final, focus_pos_final);
        % save(parameters.kwave_source_filename, 'source', 'source_labels','-v7.3');
    else
        load(parameters.kwave_source_filename);
    end    

    % Creates a sensor that records the the maximum and final pressure
    % values in every point of the grid
    sensor = struct();
    sensor.mask = ones(parameters.grid_dims);
    sensor.record = {'p_max_all','p_final'};

    % Record the last 3 cycles in steady state (when sonic waves have
    % traversed the entire medium)
    num_periods = 3;
    time_points_to_record = round(num_periods * wave_period / kgrid.dt);
    sensor.record_start_index = simulation_time_points - time_points_to_record + 1;

end