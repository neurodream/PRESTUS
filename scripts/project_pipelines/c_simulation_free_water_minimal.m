clear all; close all; clc;

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..','..'))
rootpath = pwd;

% TODO why is pn != parameters? where are params loaded in from yaml?
% pn.tuSIM = fullfile(rootpath, 'tools', 'PRESTUS');                addpath(pn.tuSIM); % why PRESTUS within PRESTUS??
pn.tuSIM = rootpath;
pn.scripts = fullfile(pn.tuSIM, 'scripts');                       addpath(pn.scripts);
pn.tuSIM_fun = fullfile(pn.tuSIM, 'functions');                   addpath(pn.tuSIM_fun);
pn.tuSIM_tools = fullfile(pn.tuSIM, 'toolboxes');                 addpath(genpath(pn.tuSIM_tools));
pn.kwave = fullfile(rootpath, 'toolboxes', 'k-wave-toolbox-version-1.4', 'k-Wave');
                                                                  addpath(pn.kwave);
pn.minimize = fullfile(rootpath, 'toolboxes', 'FEX-minimize');    addpath(pn.minimize);
pn.configs = fullfile(rootpath, 'configs');                       % why no add path??
pn.profiles = fullfile(rootpath, 'data', 'transducer_profiles');  % why no add path??
pn.data_path = fullfile(rootpath, '..', '..', 'scans');           % why no add path?? % seems like it is "only" output path? or also loading tpars etc.?
% previously: "data_higher_angle"
% it seems like data are saved within the "data_sims" folder
pn.data_seg = fullfile(rootpath, '..', '..', 'scans', 'segmentation_results');         % why no add path??
pn.nifti = (fullfile(rootpath, 'tools', 'nifti_toolbox'));        addpath(pn.nifti); % does not exist
pn.data_sims = fullfile(rootpath, 'data', 'transducer_pos'); % TODO: careful: transducer name not specified in here yet!

%% parameters
transducer_name = 'IMASONIC_IGT_300-15473_100.0mm';
% expected_focal_distance_mm = 100.0; % hardcoding here because don't want to set param in config file, which may lead to incorrect focussing down the road (?)
% TODO check if correct

distinguisher = '_tempwin1_active_noheating_free_water2'; % label to for output folder to distinguish from other sims with otherwise same parameters

targets = [98 142 149];
positions = [9 137 146];

subject_id = 3;

% TODO load them from config yaml!
dutycycle = .45; % actually 50%, but accounting for ramping .01 s each side
trials = 3;

% amp = 200000; % yet another ~average between 94229 (37.8°) and 115608 (95°)

duration = trials.*.2; % in s (aka .2 means 200 ms PRF) % TODO!!! necessary?

for i = 1:size(targets, 1)

    transducer_pos = positions(i, :);
    target_pos = targets(i, :);

    output_folder = [ ...
            'sbj_', num2str(subject_id), ...
            '_duty_', num2str(dutycycle*100), ...
            '_ntrl_', num2str(trials), ...
            '_targ_', [num2str(target_pos(1)) '-' num2str(target_pos(2)) '-' num2str(target_pos(3))], ...
            '_tpos_', [num2str(transducer_pos(1)) '-' num2str(transducer_pos(2)) '-' num2str(transducer_pos(3))], ...,
            distinguisher
            ];
    output_folder_name = fullfile(pwd, 'data', 'sims', output_folder);
        
    %% load parameter file and adjust paths if necessary
    
    % cd(fullfile(rootpath, 'data_higher_angle')) % needs to contain a folder called "configs" in this setup
    cd(fullfile(rootpath)) % needs to contain a folder called "configs" in this setup
    
    parameters = load_parameters(['config_',transducer_name,'.yaml']); % load the configuration file
    
    if ismac
        % reset paths to local mac install; NOTE: horrible not to have this be dynamic
        parameters.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
    else
        parameters.simnibs_bin_path = fullfile('/home', 'neuromod', 'marwim', 'miniconda3', 'envs', 'simnibs_env', 'bin');
    end
    
    parameters.ld_library_path ="/opt/gcc/7.2.0/lib64";
    parameters.data_path = pn.data_path;%fullfile(rootpath, 'data', 'simnibs', ['m2m_sub-00',num2str(subject_id)]);
    parameters.seg_path = pn.data_seg;%fullfile(rootpath, 'data', 'simnibs', ['m2m_sub-00',num2str(subject_id)], 'segmentation');
    parameters.temp_output_dir = output_folder_name; %pn.data_sims;
    parameters.paths_to_add = {pn.kwave, pn.minimize};
    
    if parameters.subject_subfolder
        outputs_folder = fullfile(parameters.temp_output_dir,sprintf('sub-%03d', subject_id));
    else
        outputs_folder = parameters.temp_output_dir;
    end
    if ~exist(outputs_folder, 'dir'); mkdir(outputs_folder); end
    
    
    %% Simulate again in water to check the optimization results
    % Now we will redo the simulations with new parameters. We copy the configuration, 
    % update the settings with the optimized parameters, and rerun the simulations. 

    opt_phases_and_vel = load(fullfile(pn.data_path, 'phases'));
    opt_phases_and_vel = opt_phases_and_vel.opt_phases_and_velocity;
    opt_velocity = opt_phases_and_vel(10); % TODO remove hardcoding
    opt_phases = opt_phases_and_vel(1:9);
    
    opt_source_amp_single = round(opt_velocity*parameters.medium.water.density*parameters.medium.water.sound_speed);
    opt_source_amp = repmat(opt_source_amp_single, 1, 10);

    opt_parameters = load_parameters(['config_',transducer_name,'.yaml']); % load the configuration file
    %     opt_parameters.focus_pos_t1_grid = [locs.targ_x(i), locs.targ_y(i), locs.targ_z(i)];
    %     opt_parameters.transducer.pos_t1_grid = floor([locs.trans_x(i), locs.trans_y(i), locs.trans_z(i)]);
    opt_parameters.run_heating_sims = 0;
    opt_parameters.data_path = parameters.data_path;
    opt_parameters.seg_path = parameters.seg_path;
    opt_parameters.temp_output_dir = output_folder_name;%pn.data_sims;
    opt_parameters.transducer.source_amp = opt_source_amp;
    opt_parameters.transducer.source_phase_rad = [0 opt_phases];
    opt_parameters.transducer.source_phase_deg = [0 opt_phases]/pi*180;
    % opt_parameters.expected_focal_distance_mm = parameters.expected_focal_distance_mm;
    opt_parameters.results_filename_affix = '_optimized';
    opt_parameters.simulation_medium = 'water';
    opt_parameters.overwrite_files = 'always';
    opt_parameters.interactive = 0;

    single_subject_pipeline(subject_id, opt_parameters);
    
    % %%
    % 
    % opt_phases_and_vel = load(fullfile(pn.data_path, 'phases_near_field'));
    % opt_phases_and_vel = opt_phases_and_vel.opt_phases_and_velocity;
    % opt_velocity = opt_phases_and_vel(10); % TODO remove hardcoding
    % opt_phases = opt_phases_and_vel(1:9);
    % 
    % opt_source_amp_single = round(opt_velocity*parameters.medium.water.density*parameters.medium.water.sound_speed);
    % opt_source_amp = repmat(opt_source_amp_single, 1, 10);
    % 
    % segmented_img_head = niftiinfo(fullfile(pn.data_seg, ['m2m_sub-00',num2str(subject_id)], 'final_tissues.nii.gz'));
    % pixel_size = segmented_img_head.PixelDimensions(1); % only for the case of ernie; otherwise: segmented_img_head.PixelDimensions(1)
    % 
    % %         parameters.transducer.pos_t1_grid = floor(table2array(best_trans_pos(1,["trans_x","trans_y","trans_z"])));
    % %         parameters.focus_pos_t1_grid = floor(table2array(best_trans_pos(1,["targ_x","targ_y","targ_z"])));
    % %         parameters.expected_focal_distance_mm = round(best_trans_pos.dist_to_target*pixel_size,1);
    % 
    % opt_parameters = load_parameters(['config_',transducer_name,'.yaml']);
    % % consdier different duty cycles and total durations
    % 
    % opt_parameters.simnibs_bin_path = parameters.simnibs_bin_path;
    % opt_parameters.paths_to_add = parameters.paths_to_add;
    % 
    % % % !!! here load data for insula
    % % target = insula_target_data{subject_id}.(target_name); % TODO not included in code? commented out, see few lines below
    % 
    % opt_parameters.focus_pos_t1_grid = target_pos; %[best_trans_pos.targ_x(1), best_trans_pos.targ_y(1), best_trans_pos.targ_z(1)];
    % opt_parameters.transducer.pos_t1_grid = transducer_pos; %floor([best_trans_pos.trans_x(1), best_trans_pos.trans_y(1), best_trans_pos.trans_z(1)]);
    % % opt_parameters.expected_focal_distance_mm = round(best_trans_pos.dist_to_target*pixel_size,1); % TODO!!! needed?
    % opt_parameters.data_path = pn.data_path;
    % opt_parameters.seg_path = pn.data_seg;
    % opt_parameters.t1_path_template = parameters.t1_path_template;
    % opt_parameters.t2_path_template = parameters.t2_path_template;        
    % opt_parameters.temp_output_dir = output_folder_name;%pn.data_sims;
    % 
    % opt_parameters.transducer.source_amp = opt_source_amp;
    % opt_parameters.transducer.source_phase_rad = [0 opt_phases];
    % opt_parameters.transducer.source_phase_deg = [0 opt_phases]/pi*180;
    % % opt_parameters.expected_focal_distance_mm = parameters.expected_focal_distance_mm;  % recalculated in single_subject_pipeline
    % opt_parameters.simulation_medium = 'layered'; % see default config for the list of mediums possible
    % opt_parameters.run_heating_sims = 0; % this indicates that we want the heating simulations as well
    % opt_parameters.thermal.duty_cycle = dutycycle;
    % opt_parameters.thermal.n_trials = trials;
    % opt_parameters.overwrite_files = 'always';
    % opt_parameters.overwrite_simnibs = 0;
    % opt_parameters.interactive = 0;
    % 
    % single_subject_pipeline_with_qsub(subject_id, opt_parameters, 60*60*12, 64);
    % % single_subject_pipeline(subject_id, opt_parameters);

end