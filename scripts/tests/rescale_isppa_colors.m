%% load paths

clear all;
% close all;
clc;

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
pn.seg_path = fullfile(rootpath, '..', '..', 'scans', 'segmentation_results');         % why no add path??
pn.nifti = (fullfile(rootpath, 'tools', 'nifti_toolbox'));        addpath(pn.nifti); % does not exist
pn.sim_path = fullfile(rootpath, 'data', 'transducer_pos'); % TODO: careful: transducer name not specified in here yet!

cd(fullfile(rootpath)) % needs to contain a folder called "configs" in this setup


%% load data
folder_name = 'sbj_3_duty_45_ntrl_5_amp_200000_targ_115-140-146_tpos_212-152-169';
load('M:\Documents\repos\PRESTUS_forked\data\sims\', folder_name, '\sub-003_layered_results.mat');

%% derive some more vars

% read out info from folder name
subject_id = 3;

targetStr   = regexp(data_folder, 'targ_(\d+)-(\d+)-(\d+)', 'tokens');
target = str2double(targetStr{1});
posStr      = regexp(data_folder, 'tpos_(\d+)-(\d+)-(\d+)', 'tokens');
pos = str2double(posStr{1});

[kgrid, source, sensor, source_labels] = setup_grid_source_sensor(parameters, max_sound_speed, pos, target);

parameters.data_path = 'M:\Documents\scans';
parameters.simnibs_bin_path = 'C:/Users/nicade/SimNIBS-4.0/bin';
parameters.ld_library_path = 'H:/sleep/nicade/lib64';
parameters.segmentation_software = 'charm.cmd';
parameters.output_dir = 'M:\Documents\repos\PRESTUS_forked\scripts\tests';

[medium_masks, segmented_image_cropped, skull_edge, trans_pos_final, ...
            focus_pos_final, t1_image_orig, t1_header, final_transformation_matrix, ...
            inv_final_transformation_matrix] = preprocess_brain(parameters, subject_id);

%% Process results
disp('Processing the results of acoustic simulations...')

% What is the highest pressure level for every gridpoint
data_max = gather(sensor_data.p_max_all); % gather is used since it could be a GPU array
max_pressure = max(data_max(:));

% Calculates the Isppa for every gridpoint
Isppa_map = data_max.^2./(2*(kwave_medium.sound_speed.*kwave_medium.density)).*1e-4;
% Calculates the max Isppa
max_Isppa = max(Isppa_map(:));

% Calculates the Mechanical Index for every gridpoint
MI_map = (data_max/10^6)/sqrt((parameters.transducer.source_freq_hz/10^6));

% Creates the foundation for a mask before the exit plane to calculate max values outside of it
comp_grid_size = size(sensor_data.p_max_all);
after_exit_plane_mask = ones(comp_grid_size);

bowl_depth_grid = round((parameters.transducer.curv_radius_mm-parameters.transducer.dist_to_plane_mm)/parameters.grid_step_mm);

% Places the exit plane mask in the grid, adjusted to the amount of dimensions
if parameters.n_sim_dims == 3
    if trans_pos_final(3) > comp_grid_size(3)/2
        after_exit_plane_mask(:,:,(trans_pos_final(parameters.n_sim_dims)-bowl_depth_grid):end) = 0;
    else
        after_exit_plane_mask(:,:,1:(trans_pos_final(parameters.n_sim_dims)+bowl_depth_grid)) = 0;
    end
else
    if trans_pos_final(2) > comp_grid_size(2)/2
        after_exit_plane_mask(:,(trans_pos_final(parameters.n_sim_dims)-bowl_depth_grid):end) = 0;
    else
        after_exit_plane_mask(:,1:(trans_pos_final(parameters.n_sim_dims)+bowl_depth_grid)) = 0;
    end
end

% Calculates the X, Y and Z coordinates of the max. intensity
[max_Isppa_after_exit_plane, Ix_eplane, Iy_eplane, Iz_eplane] = masked_max_3d(Isppa_map, after_exit_plane_mask);

% Combines these coordinates into a point of max. intensity in the grid
if parameters.n_sim_dims==3
    max_isppa_eplane_pos = [Ix_eplane, Iy_eplane, Iz_eplane];
else 
    max_isppa_eplane_pos = [Ix_eplane, Iy_eplane];
end
disp('Final transducer, expected focus, and max ISPPA positions')

% Calculates the average Isppa within a circle around the target
[trans_pos_final', focus_pos_final', max_isppa_eplane_pos']
real_focal_distance = norm(max_isppa_eplane_pos-trans_pos_final)*parameters.grid_step_mm;
avg_radius = round(parameters.focus_area_radius/parameters.grid_step_mm); %grid
avg_isppa_around_target = Isppa_map((focus_pos_final(1)-avg_radius):(focus_pos_final(1)+avg_radius),...
    (focus_pos_final(2)-avg_radius):(focus_pos_final(2)+avg_radius),...
    (focus_pos_final(3)-avg_radius):(focus_pos_final(3)+avg_radius));
avg_isppa_around_target = mean(avg_isppa_around_target(:));

% Reports the Isppa within the original stimulation target
isppa_at_target = Isppa_map(focus_pos_final(1),focus_pos_final(2),focus_pos_final(3));

% Creates a logical skull mask and register skull_ids
labels = fieldnames(parameters.layer_labels);
skull_i = find(strcmp(labels,  'skull_cortical'));
trabecular_i = find(strcmp(labels,  'skull_trabecular'));
all_skull_ids = [skull_i, trabecular_i];
skull_mask = ismember(medium_masks,all_skull_ids);

% Overwrites the max Isppa by dividing it up into the max Isppa for
% each layer in case a layered simulation_medium was selected
if contains(parameters.simulation_medium, 'skull') || strcmp(parameters.simulation_medium, 'layered')
    [max_Isppa_brain, Ix_brain, Iy_brain, Iz_brain] = masked_max_3d(Isppa_map, medium_masks>0 & medium_masks<3);
    half_max = Isppa_map >= max_Isppa_brain/2 & medium_masks>0 & medium_masks<3;
    half_max_ISPPA_volume_brain = sum(half_max(:))*(parameters.grid_step_mm^3);

    [max_pressure_brain, Px_brain, Py_brain, Pz_brain] = masked_max_3d(data_max, medium_masks>0 & medium_masks<3);
    [max_MI_brain, Px_brain, Py_brain, Pz_brain] = masked_max_3d(MI_map, medium_masks>0 & medium_masks<3);
    [max_Isppa_skull, Ix_skull, Iy_skull, Iz_skull] = masked_max_3d(Isppa_map, skull_mask);
    [max_pressure_skull, Px_skull, Py_skull, Pz_skull] = masked_max_3d(data_max, skull_mask);
    [max_MI_skull, Px_skull, Py_skull, Pz_skull] = masked_max_3d(MI_map, skull_mask);
    [max_Isppa_skin, Ix_skin, Iy_skin, Iz_skin] = masked_max_3d(Isppa_map, medium_masks==5);
    [max_pressure_skin, Px_skin, Py_skin, Pz_skin] = masked_max_3d(data_max, medium_masks==5);
    [max_MI_skin, Px_skin, Py_skin, Pz_skin] = masked_max_3d(MI_map, medium_masks==5);
    highlighted_pos = [Ix_brain, Iy_brain, Iz_brain];
    real_focal_distance = norm(highlighted_pos-trans_pos_final)*parameters.grid_step_mm;

    writetable(table(subject_id, max_Isppa, max_Isppa_after_exit_plane, real_focal_distance, max_Isppa_skin, max_Isppa_skull, max_Isppa_brain, max_pressure_skin, max_pressure_skull, max_pressure_brain, max_MI_skin, max_MI_skull, max_MI_brain, Ix_brain, Iy_brain, Iz_brain, trans_pos_final, focus_pos_final, isppa_at_target, avg_isppa_around_target, half_max_ISPPA_volume_brain), output_pressure_file);
else % If no layered tissue was selected, the max Isppa is highlighted on the plane and written in a table.
    max_Isppa = max(Isppa_map(:)); %  Does this step need to be included? already done at line 225.
    highlighted_pos = max_isppa_eplane_pos;
    writetable(table(subject_id, max_Isppa, max_Isppa_after_exit_plane, max_pressure, real_focal_distance, trans_pos_final, focus_pos_final, isppa_at_target, avg_isppa_around_target), output_pressure_file);
end

% Plots the Isppa on the segmented image
output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_isppa%s.png', subject_id, parameters.simulation_medium, parameters.results_filename_affix));

if parameters.n_sim_dims==3
    [~,~,~,~,~,~,~,h]=plot_isppa_over_image(Isppa_map, segmented_image_cropped, source_labels, parameters, {'y', focus_pos_final(2)}, trans_pos_final, focus_pos_final, highlighted_pos);
else
    [~,~,h]=plot_isppa_over_image_2d(Isppa_map, segmented_image_cropped, source_labels, parameters,  trans_pos_final, focus_pos_final, highlighted_pos);
end
saveas(h, output_plot, 'png')
close(h);