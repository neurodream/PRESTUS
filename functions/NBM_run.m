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

    %% set transducer position
    parameters.transducers(i).name = transducer_labels{i};
    [parameters, distance] = get_transducer_pos(parameters, subject_id, dirs{i}, i, angles(i,:), transd_pos_shift(i,:), focus_pos_shift(i,:), contralateral(i));
    disp(distance);

    % % TODO debug delete
    % if sham
    %     distance = 20;
    % end

    %% set transducer phases

    % comment out when loading config results from PRESTUS acoustic profiling

    % % TODO figure out which optimization works best
    % if contralateral(i)
    %     parameters = calculate_transducer_phases(parameters, i, distance, 15, goal_intensity, sham, desired_function, use_all_phases); % distance + 30 % 28
    % elseif ~contralateral(i)
    %     parameters = calculate_transducer_phases(parameters, i, distance, 15, goal_intensity, sham, desired_function, use_all_phases);
    % end

    %% debug plots



            % % close all;
            % 
            % figure;
            % 
            % axial_position = (1:parameters.default_grid_dims(3))*(parameters.grid_step_mm);
            % axial_position = axial_position';
            % 
            % transducer = parameters.transducers(i);
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





    %% store optimal params

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

% TODO DELETE!!!! hardcoding entry points
parameters.transducers(1).pos_t1_grid = round([19.3652705981541 193.852823508941 150.859066169943]);
parameters.transducers(2).pos_t1_grid = round([204.755585045787 179.393318458352 154.836270127952]);

% TODO DELETE!!!! just for iteration 11: check acoustic sims of traveling
% wave
% parameters.transducers(1).source_amp = zeros(1,10);

% parameters.transducers(1).source_amp = [0 0 0 0 0 0 0 0 0 0]; % TODO debug delete!


% TODO debug remove hardcoding!! below values are for:
% % it 14
%   % set_focus_wrt_exit_plane_mm: 92.5
%   % set_intensity_w_per_cm2: 30.0
% parameters.transducers(1).source_phase_deg = [0.0, 5.827660709755711, 1.9361225024631727, 22.27213694131474, 11.151622815763867, 35.857675317466324, 14.803867876994081, 31.801708627788784, 23.750892901091536, 0.06416493177429752];
% parameters.transducers(1).source_amp = [160014.0, 160014.0, 160014.0, 160014.0, 160014.0, 160014.0, 160014.0, 160014.0, 160014.0, 160014.0];
% parameters.transducers(1).source_phase_rad = [0.0, 0.10171186707434678, 0.033791712389934385, 0.3887221210810013, 0.1946325350756005, 0.6258344964008904, 0.2583762364837696, 0.5550445233148021, 0.4145312814125952, 0.001119889323779019];
% parameters.transducers(1).optim_params.focal_distance_mm = 64.7;
% % set_focus_wrt_exit_plane_mm: 64.7
%   % set_intensity_w_per_cm2: 30.0
% parameters.transducers(2).source_phase_deg = [0.0, 17.536292796409906, 1.5600840619362728, 30.926680645084335, 1.2575583455341289E-4, 168.2734166541756, 71.34989806928922, 135.4452709665064, 99.57003960738558, 105.14036588217982];
% parameters.transducers(2).source_amp = [154100.0, 154100.0, 154100.0, 154100.0, 154100.0, 154100.0, 154100.0, 154100.0, 154100.0, 154100.0];
% parameters.transducers(2).source_phase_rad = [0.0, 0.3060660478911165, 0.0272286034886751, 0.5397724039695255, 2.194853366550308E-6, 2.9369251641956247, 1.245290642271442, 2.363965934621419, 1.7378250274900402, 1.8350455613955503];
% parameters.transducers(2).optim_params.focal_distance_mm = 92.5;

% % it 15
%   % set_focus_wrt_exit_plane_mm: 92.5
%   % set_intensity_w_per_cm2: 30.0 (at target)
% parameters.transducers(1).source_phase_deg = [0.0, 5.827660709755711, 1.9361225024631727, 22.27213694131474, 11.151622815763867, 35.857675317466324, 14.803867876994081, 31.801708627788784, 23.750892901091536, 0.06416493177429752];
% parameters.transducers(1).source_amp = [160014.0, 160014.0, 160014.0, 160014.0, 160014.0, 160014.0, 160014.0, 160014.0, 160014.0, 160014.0];
% parameters.transducers(1).source_phase_rad = [0.0, 0.10171186707434678, 0.033791712389934385, 0.3887221210810013, 0.1946325350756005, 0.6258344964008904, 0.2583762364837696, 0.5550445233148021, 0.4145312814125952, 0.001119889323779019];
% parameters.transducers(1).optim_params.focal_distance_mm = 80;
% % set_focus_wrt_exit_plane_mm: 80
%   % set_intensity_w_per_cm2: 30.0 (at target)
% parameters.transducers(2).source_phase_deg = [0.0, 22.986714640286568, 1.3522658189916815E-5, 34.271896086853005, 2.649515970852959E-6, 70.87202307916031, 7.665689139421615E-6, 88.50704006939364, 40.42272452929138, 47.58669976972975];
% parameters.transducers(2).source_amp = [157741.0, 157741.0, 157741.0, 157741.0, 157741.0, 157741.0, 157741.0, 157741.0, 157741.0, 157741.0];
% parameters.transducers(2).source_phase_rad = [0.0, 0.40119385468938457, 2.3601490903582507E-7, 0.5981574276169455, 4.624277727555825E-8, 1.2369501502807574, 1.3379151491727783E-7, 1.5447392604054695, 0.7055096356628097, 0.8305445911398115];
% parameters.transducers(2).optim_params.focal_distance_mm = 92.5;
% 
% % it16
% % set_focus_wrt_exit_plane_mm: 73
%   % set_intensity_w_per_cm2: 30.0 (at target)
% parameters.transducers(2).source_phase_deg = [0.0, 14.191294489845086, 25.5888333159838, 51.98042830197692, 49.149681722670095, 83.31532820169576, 79.04210285392764, 113.61493450384977, 121.9442997341297, 98.10876821663902];
% parameters.transducers(2).source_amp = [130289.0, 130289.0, 130289.0, 130289.0, 130289.0, 130289.0, 130289.0, 130289.0, 130289.0, 130289.0];
% parameters.transducers(2).source_phase_rad = [0.0, 0.24768481396792574, 0.4466093931079358, 0.9072296204663426, 0.8578237723678717, 1.4541267944992773, 1.379544942501044, 1.9829546865298884, 2.128329534399409, 1.7123210304563163];
% parameters.transducers(2).optim_params.focal_distance_mm = 73;
% 
% source_amp = parameters.transducers(2).source_amp;
% 
% % it17
% % parameters.transducers(2).source_amp = source_amp*1.2;
% 
% % it18
% parameters.transducers(2).source_amp = source_amp*0.8;
% 
% % it19
%   % set_intensity_w_per_cm2: 40.0 (at target)
% parameters.transducers(1).source_phase_deg = [0.0, 343.64904294320337, 0.0028029198905977666, 14.890625268016658, 0.016463736713217453, 28.956908347169765, 11.894798417462951, 28.51115287777094, 359.99811120899284, 1.5344752115046985];
% parameters.transducers(1).source_amp = [185444.0, 185444.0, 185444.0, 185444.0, 185444.0, 185444.0, 185444.0, 185444.0, 185444.0, 185444.0];
% parameters.transducers(1).source_phase_rad = [0.0, 5.997807270686284, 4.8920180760570275E-5, 0.2598904386075538, 2.87346412827114E-4, 0.5053933918563416, 0.20760339624573948, 0.49761349125655924, 6.28315234155763, 0.02678164473099336];
% parameters.transducers(1).optim_params.focal_distance_mm = 92.5;
% 
% parameters.transducers(2).source_phase_deg = [0.0, 18.014986157747888, 0.0030833932256510137, 36.02496863292742, 1.6160459090668836E-5, 159.3958710517933, 64.65566348508041, 137.39509750300215, 99.91101457873015, 106.4379712181883];
% parameters.transducers(2).source_amp = [175765.0, 175765.0, 175765.0, 175765.0, 175765.0, 175765.0, 175765.0, 175765.0, 175765.0, 175765.0];
% parameters.transducers(2).source_phase_rad = [0.0, 0.3144208231539032, 5.381536392129866E-5, 0.628754315572264, 2.8205321976602E-7, 2.7819827639381103, 1.1284542078761248, 2.3979968275260273, 1.7437761634068965, 1.8576930469003456];
% parameters.transducers(2).optim_params.focal_distance_mm = 65;

% TODO why focus_pos_t1_grid between transducers not exactly the same? (at
% least the one we just checked differed in z by 2 voxels)
% more specifically, the left transducer focus_pos_t1_grid is lower by 2
% voxels than set in the LUT

% % Run the pipeline
single_subject_pipeline_with_slurm(subject_id, parameters, "08:00:00", 60);
% single_subject_pipeline_with_qsub(subject_id, parameters);
% % single_subject_pipeline(subject_id, parameters); % TODO change back or keep commented!!

% store the parameters for debugging
save(fullfile(parameters.data_path, 'sim_outputs', [sprintf('sub-%03d/sub-%03d', subject_id, subject_id) '_parameters' ID]))

end