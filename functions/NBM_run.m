function [] = NBM_run(subject_id, focus_side, sham, z_shift, goal_intensity, extra_ID_suffix, parameters_fname, imprecision_modeling)

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

if z_shift > 0
    z_shift_part_ID = '+z';
else
    z_shift_part_ID = '';
end

if sham
    sham_part_ID = '-';
else
    sham_part_ID = '--';
end

ID_part = ['L' z_shift_part_ID sham_part_ID focus_sides_add_to_ID{1} '_' 'R' z_shift_part_ID sham_part_ID focus_sides_add_to_ID{2}];
dirs                = {'l',         'r'}; % i.e. target side (note: just for readout of the coordinates that determine the beam vector! the _distance_ will determine ipsi vs. contralateral
angles              = [-1 0 0;    1 0 0];
transd_pos_shift    = [0 0 z_shift;       0 0 z_shift];
focus_pos_shift     = [0 0 z_shift;       0 0 z_shift];

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


% transducer "preprocessing" (TODO find better word)
for i = 1:numel(parameters.transducers)

    parameters.transducers(i).name = transducer_labels{i};
    [parameters, distance] = get_transducer_pos(parameters, subject_id, dirs{i}, i, angles(i,:), transd_pos_shift(i,:), focus_pos_shift(i,:), contralateral(i));
    
    % TODO debug delete
    if sham
        distance = 20;
    end

    % TODO figure out which optimization works best
    if contralateral(i)
        parameters = calculate_transducer_phases(parameters, i, distance, 15, goal_intensity, sham); % distance + 30 % 28
    elseif ~contralateral(i)
        parameters = calculate_transducer_phases(parameters, i, distance, 15, goal_intensity, sham);
    end

    % store the indended parameters for later debugging:
    parameters.transducers(i).optim_params = [];
    parameters.transducers(i).optim_params.focal_distance_mm = distance;
    parameters.transducers(i).optim_params.angle = angles(i,:);
    parameters.transducers(i).optim_params.transd_pos_shift = transd_pos_shift(i,:);
    parameters.transducers(i).optim_params.focus_pos_shift = focus_pos_shift(i,:);

end

% % TODO debug check the visuals
% plot_transducer_pos(parameters, subject_id, true, false, false, false)

% add field of free water axial intensity to structs
parameters = get_simulated_axial_intensity(parameters);

% Set the results filename
parameters.results_filename_affix = ID;

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





% Run the pipeline
single_subject_pipeline_with_slurm(subject_id, parameters, "08:00:00");
% single_subject_pipeline_with_qsub(subject_id, parameters);
% single_subject_pipeline(subject_id, parameters); % TODO change back or keep commented!!

% store the parameters for debugging
save(fullfile(parameters.data_path, 'sim_outputs', [sprintf('sub-%03d/sub-%03d', subject_id, subject_id) '_parameters' ID]))

end