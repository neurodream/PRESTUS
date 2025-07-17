function [] = NBM_check_loc(subject_id, focus_side, sham, z_shift, goal_intensity, extra_ID_suffix, parameters_fname, imprecision_modeling)

% function streamlines parameter adjustment and starts the simulation;
% function optimized for NBM study, or for simulating crossbeam TUS for 
% any targets that exist bilaterally and where a symmetric placement of 
% 2 transducers is desired
%
% inputs:
% - subject_id:           integer number, needs to match the existing scans
% - focus_side:           'L', 'R', 'BL'
% - sham:                 bool, whether only near-field (sham == true) or whole field (sham == false) is simulated
% - z_shift:              number in voxels that transducer should be translated upwards (positive values) or downwards (negative values) to correct for skull refractions observed from initial NBM simulations (different values per individual)
% - goal_intensity:       in ?? (TODO)
% - extra_ID_suffix:      string to discriminate/more precisely label this particular sim run
% - parameters_fname:     name of yaml file(s) without path (needs to be located in configs directory) or file extension
% - imprecision_modeling: location where imprecision is modeled ('transducer', 'focus', 'both', 'none')


if size(z_shift, 2) > 1
    xyz_shift = z_shift;
    x_shift = xyz_shift(1);
    y_shift = xyz_shift(2);
    z_shift = xyz_shift(3);
else
    x_shift = 0;
    y_shift = 0;
end

transducer_labels   = {'L',         'R'};
if strcmp(focus_side, 'L')
    contralateral       = [false,        true];
    focus_sides_add_to_ID = {'l', 'l'};
elseif strcmp(focus_side, 'R')
    contralateral       = [true,        false];
    focus_sides_add_to_ID = {'r', 'r'};
else
    contralateral       = [false,        false];
    focus_sides_add_to_ID = {'l', 'r'};
end % TODO add middle focus

x_shift_part_ID = '';
y_shift_part_ID = '';
z_shift_part_ID = '';
% TODO better to add the number
if x_shift > 0
    x_shift_part_ID = '+x';
elseif x_shift < 0
    x_shift_part_ID = '-x';
end
if y_shift > 0
    y_shift_part_ID = '+y';
elseif y_shift < 0
    y_shift_part_ID = '-y';
end
if z_shift > 0
    z_shift_part_ID = '+z';
elseif z_shift < 0
    z_shift_part_ID = '-z';
end

if sham
    sham_part_ID = '-';
else
    sham_part_ID = '--';
end

ID_part = ['L' z_shift_part_ID sham_part_ID focus_sides_add_to_ID{1} '_' 'R' x_shift_part_ID y_shift_part_ID z_shift_part_ID sham_part_ID focus_sides_add_to_ID{2}];
dirs                = {'l',         'r'}; % i.e. target side (note: just for readout of the coordinates that determine the beam vector! the _distance_ will determine ipsi vs. contralateral
angles              = [-1 0 0;    1 0 0];
transd_pos_shift    = [x_shift y_shift z_shift;       x_shift y_shift z_shift];
focus_pos_shift     = [x_shift y_shift z_shift;       x_shift y_shift z_shift];

desired_function = @create_boxcar; % TODO back to established profiles
use_all_phases = true;

% load config(s)
if iscellstr(parameters_fname)
    parameters_fnames = {};
    for i = 1:numel(parameters_fname)
        fname = parameters_fname{i};
        parameters_fnames{end+1} = [fname '.yaml'];
        
    end
    parameters = load_parameters(parameters_fnames{:});
else
    parameters = load_parameters([parameters_fname '.yaml']);
end

ID = [ID_part '_' extra_ID_suffix '_imprecision' imprecision_modeling];


% % TODO with segmentation only, instead of the below, do this quick fix of -
% % comment out of segmentation already exists!
% % transducer placement
% parameters.transducers(1).pos_t1_grid = [10 10 10];
% parameters.transducers(1).focus_pos_t1_grid = [50 50 50];  % flip back to normal space
% parameters.transducers(2).pos_t1_grid = [100 100 100];
% parameters.transducers(2).focus_pos_t1_grid = [60 60 60];  % flip back to normal space

% transducer "preprocessing" (TODO find better word) (TODO disable when running segmentation only)
for i = 1:numel(parameters.transducers)

    parameters.transducers(i).name = transducer_labels{i};
    [parameters, distance] = get_transducer_pos(parameters, subject_id, dirs{i}, i, angles(i,:), transd_pos_shift(i,:), focus_pos_shift(i,:), contralateral(i));
    disp(distance);

    % % TODO debug delete
    % if sham
    %     distance = 20;
    % end

    % TODO figure out which optimization works best
    if contralateral(i)
        parameters = calculate_transducer_phases(parameters, i, distance, 15, goal_intensity, sham, desired_function, use_all_phases); % distance + 30 % 28
    elseif ~contralateral(i)
        parameters = calculate_transducer_phases(parameters, i, distance, 15, goal_intensity, sham, desired_function, use_all_phases);
    end



            % % close all;
            % 
            % figure;
            % 
            % axial_position = (1:parameters.default_grid_dims(3))*(parameters.grid_step_mm);
            % axial_position = axial_position';
            % 
            % transducer = parameters.transducers(1);
            % 
            % dist_to_exit_plane = transducer.curv_radius_mm - transducer.dist_to_plane_mm;
            % ax_pos = axial_position + dist_to_exit_plane;
            % 
            % % TODO readout the intensity correctly
            % p_axial_oneil = focusedAnnulusONeil( ...
            %     transducer.curv_radius_mm/1e3, ...
            %     [transducer.Elements_ID_mm; transducer.Elements_OD_mm]/1e3, ...
            %     transducer.source_amp/(parameters.medium.water.density*parameters.medium.water.sound_speed), ...
            %     transducer.source_phase_rad, ...
            %     transducer.source_freq_hz, ...
            %     parameters.medium.water.sound_speed, ...
            %     parameters.medium.water.density, ...
            %     (ax_pos-0.5)*1e-3 ...
            %     );
            % 
            % i_axial_oneil = p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;
            % 
            % plot(axial_position, i_axial_oneil)
            % % ylim([0 100])
            % 
            % 
            % figure;
            % 
            % axial_position = (1:parameters.default_grid_dims(3))*(parameters.grid_step_mm);
            % axial_position = axial_position';
            % 
            % transducer = parameters.transducers(2);
            % 
            %     dist_to_exit_plane = transducer.curv_radius_mm - transducer.dist_to_plane_mm;
            %     ax_pos = axial_position + dist_to_exit_plane;
            % 
            %     opt_limits = [ax_pos(2,1), ax_pos(end,1)];
            % 
            %     % TODO readout the intensity correctly
            %     p_axial_oneil = focusedAnnulusONeil( ...
            %         transducer.curv_radius_mm/1e3, ...
            %         [transducer.Elements_ID_mm; transducer.Elements_OD_mm]/1e3, ...
            %         transducer.source_amp/(parameters.medium.water.density*parameters.medium.water.sound_speed), ...
            %         transducer.source_phase_rad, ...
            %         transducer.source_freq_hz, ...
            %         parameters.medium.water.sound_speed, ...
            %         parameters.medium.water.density, ...
            %         (ax_pos-0.5)*1e-3 ...
            %         );
            % 
            %     i_axial_oneil = p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;
            % 
            % plot(axial_position, i_axial_oneil)
            % % ylim([0 100])







    % store the indended parameters for later debugging:
    parameters.transducers(i).optim_params = [];
    parameters.transducers(i).optim_params.focal_distance_mm = distance;
    parameters.transducers(i).optim_params.angle = angles(i,:);
    parameters.transducers(i).optim_params.transd_pos_shift = transd_pos_shift(i,:);
    parameters.transducers(i).optim_params.focus_pos_shift = focus_pos_shift(i,:);

end

% Set the results filename
parameters.results_filename_affix = ID;

% TODO debug check the visuals
plot_transducer_pos(parameters, subject_id, false, 'Functional', 'mechanicalindex', 'Structural', 'scalp');

% add field of free water axial intensity to structs
parameters = get_simulated_axial_intensity(parameters);

%% Loop through each transducer and apply the imprecision variations
for i = 1:numel(parameters.transducers)
    % Get current position and focus for the transducer
    t = parameters.transducers(i).pos_t1_grid;
    f = parameters.transducers(i).focus_pos_t1_grid;

    vt = round(-3 + 6 * rand(1, 3));  % Random values for pos_t1_grid
    vf = round(-3 + 6 * rand(1, 3));  % Random values for focus_pos_t1_grid

    % Conditionally set vt and vf based on the 'vary' parameter
    switch imprecision_modeling
        case 'transducer'
            vf = [0, 0, 0];            % No variation in focus_pos_t1_grid
        case 'focus'
            vt = [0, 0, 0];            % No variation in pos_t1_grid
        case 'none'
            vf = [0, 0, 0];
            vt = [0, 0, 0];
    end

    % Update pos_t1_grid and focus_pos_t1_grid with the calculated values
    parameters.transducers(i).pos_t1_grid = t + vt;
    parameters.transducers(i).focus_pos_t1_grid = f + vf;
end

% TODO DELETE!!!! just for iteration 11: check acoustic sims of traveling
% wave
% parameters.transducers(1).source_amp = zeros(1,10);




% % Run the pipeline
% single_subject_pipeline_with_slurm(subject_id, parameters, "08:00:00", 60);
% % single_subject_pipeline_with_qsub(subject_id, parameters);
% % single_subject_pipeline(subject_id, parameters); % TODO change back or keep commented!!
% 
% % store the parameters for debugging
% save(fullfile(parameters.data_path, 'sim_outputs', [sprintf('sub-%03d/sub-%03d', subject_id, subject_id) '_parameters' ID]))

end