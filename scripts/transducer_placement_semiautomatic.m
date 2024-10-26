%% load data

close all; clear; clc;

% TODO: flip x and y (careful: then also in the other scripts!)

sbj_ID = 5;
angle_L = [0.2 -1 0];
angle_R = [0.2  1 0];

% adjust to your path
cd /home/sleep/nicade/Documents/repos/PRESTUS_forked/;

addpath('functions')
addpath(genpath('toolboxes'))

% automatic readout of coordinates (TODO: read out from table)
all_files = dir(fullfile(pwd, 'data/transducer_pos', 'tpars*.csv'));
sbj_targets = []; %struct(numel(all_files));
for sID = 1:numel(all_files)/2
    matching_files = all_files(contains({all_files.name}, sprintf('%03d', sID)));
    for direction = {'left', 'right'}
        matching_dir_file = matching_files(contains({matching_files.name}, direction));
        T = readtable(fullfile(pwd, 'data/transducer_pos', matching_dir_file.name));
        sbj_targets(sID).(direction{1}) = [T.targ_y(1) T.targ_x(1) T.targ_z(1)];
    end
    disp(['subject ' num2str(sID) ' done.']);
end
clear all_files sID matching_files direction matching_dir_file T

%%

target_L = sbj_targets(sbj_ID).left;
target_R = sbj_targets(sbj_ID).right;

% make sure the subject ID match of seg_file and target:
seg_file = ['/home/sleep/nicade/Documents/scans/segmentation_results/m2m_sub-' sprintf('%03d', sbj_ID) '/final_tissues.nii.gz'];

layers = niftiread(seg_file);
layers_info = niftiinfo(seg_file);
head = layers > 0;

head = fill_head(head);

skull = layers == 7 | layers == 8;

parameters = load_parameters_minimal('config_IMASONIC_IGT_300-15473_100.0mm.yaml');
transformMatrix = layers_info.Transform.T;
parameters.transform = transformMatrix;
parameters.grid_step_mm = mean([transformMatrix(1,1) transformMatrix(2,2) transformMatrix(3,3)]); % TODO check in Julian's code if valid

clear seg_file transformMatrix

%% calculate effective targets

% TODO remove the hardcoding in the following:

% should be target L or R, but remove/add 5 in two dimensions
% defacto_target_L = [target_L(1), target_L(2) - 5, target_L(3) - 5]; %[144 91 144];
% defacto_target_R = [target_R(1), target_R(2) + 5, target_R(3) - 5]; %[144 135 145];
defacto_target_L = [target_L(1), target_L(2), target_L(3)]; %[144 91 144];
defacto_target_R = [target_R(1), target_R(2), target_R(3)]; %[144 135 145];

% location of the center of the exit plane (?) of the transducer
[~, source_L] = get_transducer_voxels(defacto_target_L, angle_L, head, parameters, '', 'red');
[~, source_R] = get_transducer_voxels(defacto_target_R, angle_R, head, parameters, '', 'red');
% source_L = [144.0000 24.6392 169.2794];

% left transducer
% left NBM
% effective_target_L_ipsi = point_line_projection_3D(target_L, source_L, defacto_target_R);
% right NBM
effective_target_L_contra = point_line_projection_3D(defacto_target_R, source_L, target_L);

% right transducer
% left NBM
effective_target_R_contra = point_line_projection_3D(defacto_target_L, source_R, target_R);
% right NBM
% effective_target_R_ipsi = point_line_projection_3D(target_R, source_R, defacto_target_L);

close all; clc;

%% create figure with head/skull and targets

figure;
p = patch(isosurface(head, 0.5)); % Extract and plot outer layer
set(p, 'FaceAlpha', 0.5, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none'); % Customize appearance
isonormals(head, p); % Add normals for proper lighting

hold on;

p_skull = patch(isosurface(skull, 0.5)); % Extract and plot outer layer
set(p_skull, 'FaceAlpha', 0, 'FaceColor', 'black', 'EdgeColor', 'none'); % Customize appearance
isonormals(skull, p_skull); % Add normals for proper lighting

plot3(target_L(1), target_L(2), target_L(3), 'k.', 'MarkerSize', 20); % red point at the target
plot3(target_R(1), target_R(2), target_R(3), 'k.', 'MarkerSize', 20); % green point at the left NBM

xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
axis equal;
view(3);
camlight;
lighting gouraud;

% % plot beam path from left and right transducer, respectively
% plot3([source_L(1), defacto_target_R(1)], [source_L(2), defacto_target_R(2)], [source_L(3), defacto_target_R(3)], 'k-', 'LineWidth', 1);
% plot3([source_R(1), defacto_target_L(1)], [source_R(2), defacto_target_L(2)], [source_R(3), defacto_target_L(3)], 'k-', 'LineWidth', 1);

%% add transducer

% remove a transducer by calling: delete(findobj('Tag', transducer_name));

get_transducer_voxels(defacto_target_L, angle_L, head, parameters, 'transd1', '#0072BD');
get_transducer_voxels(effective_target_R_contra, angle_R, head, parameters, 'transd2', '#A2142F');

view(62, 36);

%% decrepit

% get_transducer_voxels([144, 96, 149], [0.2 -1 0.22], head, parameters, 'suggestion1', [0.1 0.1 0.1]);
% get_transducer_voxels([144, 96, 149], [0.46 -1 0.38], head, parameters, 'suggestion2', [0 0 0]);
% get_transducer_voxels([144, 96, 149], [0.44 -1 0.26], head, parameters, 'suggestion3', [1 0 0]);

% params = load(fullfile('M:\Documents\scans\sim_outputs\sub-001', 'sub-001_parametersnearparallel_var1_151024_1744'));
% get_transducer_voxels_absolute_pos( ...
%     params.parameters.transducers(1).focus_pos_t1_grid([2 1 3]), ...
%     params.parameters.transducers(1).pos_t1_grid([2 1 3]), ...
%     head, parameters, 'transdL1', '#A2142F');
% get_transducer_voxels_absolute_pos( ...
%     params.parameters.transducers(2).focus_pos_t1_grid([2 1 3]), ...
%     params.parameters.transducers(2).pos_t1_grid([2 1 3]), ...
%     head, parameters, 'transdR1', '#0072BD');
% params = load(fullfile('M:\Documents\scans\sim_outputs\sub-001', 'sub-001_parametersnearparallel_var2_151024_1744'));
% get_transducer_voxels_absolute_pos( ...
%     params.parameters.transducers(1).focus_pos_t1_grid([2 1 3]), ...
%     params.parameters.transducers(1).pos_t1_grid([2 1 3]), ...
%     head, parameters, 'transdL2', '#A2142F');
% get_transducer_voxels_absolute_pos( ...
%     params.parameters.transducers(2).focus_pos_t1_grid([2 1 3]), ...
%     params.parameters.transducers(2).pos_t1_grid([2 1 3]), ...
%     head, parameters, 'transdR2', '#0072BD');
% params = load(fullfile('M:\Documents\scans\sim_outputs\sub-001', 'sub-001_parametersnearparallel_var3_151024_1744'));
% get_transducer_voxels_absolute_pos( ...
%     params.parameters.transducers(1).focus_pos_t1_grid([2 1 3]), ...
%     params.parameters.transducers(1).pos_t1_grid([2 1 3]), ...
%     head, parameters, 'transdL3', '#A2142F');
% get_transducer_voxels_absolute_pos( ...
%     params.parameters.transducers(2).focus_pos_t1_grid([2 1 3]), ...
%     params.parameters.transducers(2).pos_t1_grid([2 1 3]), ...
%     head, parameters, 'transdR3', '#0072BD');
% params = load(fullfile('M:\Documents\scans\sim_outputs\sub-001', 'sub-001_parametersnearparallel_var4_151024_1744'));
% get_transducer_voxels_absolute_pos( ...
%     params.parameters.transducers(1).focus_pos_t1_grid([2 1 3]), ...
%     params.parameters.transducers(1).pos_t1_grid([2 1 3]), ...
%     head, parameters, 'transdL4', '#A2142F');
% get_transducer_voxels_absolute_pos( ...
%     params.parameters.transducers(2).focus_pos_t1_grid([2 1 3]), ...
%     params.parameters.transducers(2).pos_t1_grid([2 1 3]), ...
%     head, parameters, 'transdR4', '#0072BD');
% params = load(fullfile('M:\Documents\scans\sim_outputs\sub-001', 'sub-001_parametersnearparallel_var5_151024_1744'));
% get_transducer_voxels_absolute_pos( ...
%     params.parameters.transducers(1).focus_pos_t1_grid([2 1 3]), ...
%     params.parameters.transducers(1).pos_t1_grid([2 1 3]), ...
%     head, parameters, 'transdL5', '#A2142F');
% get_transducer_voxels_absolute_pos( ...
%     params.parameters.transducers(2).focus_pos_t1_grid([2 1 3]), ...
%     params.parameters.transducers(2).pos_t1_grid([2 1 3]), ...
%     head, parameters, 'transdR5', '#0072BD');
% params = load(fullfile('M:\Documents\scans\sim_outputs\sub-001', 'sub-001_parametersnearparallel_var6_151024_1744'));
% get_transducer_voxels_absolute_pos( ...
%     params.parameters.transducers(1).focus_pos_t1_grid([2 1 3]), ...
%     params.parameters.transducers(1).pos_t1_grid([2 1 3]), ...
%     head, parameters, 'transdL6', '#A2142F');
% get_transducer_voxels_absolute_pos( ...
%     params.parameters.transducers(2).focus_pos_t1_grid([2 1 3]), ...
%     params.parameters.transducers(2).pos_t1_grid([2 1 3]), ...
%     head, parameters, 'transdR6', '#0072BD');

% give the transducer a name

% %parameters001
% transd1_pos = [151,213,147];
% transd2_pos = [137,18,146];
% focus1_pos = [144,130,150];
% focus2_pos = [144,96,149];

% %parameters002
% transd1_pos = [144,200,168];
% transd2_pos = [144,25,169];
% focus1_pos = [144,91,144];
% focus2_pos = [144,135,145];

% transducer_voxels = get_transducer_voxels(effective_target_L_contra, [0 -1 0.22], head, parameters, transducer_name);