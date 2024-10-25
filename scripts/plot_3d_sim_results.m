clear all; % clc; % close all

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
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

%% positions

pos_NBM = [98 142 149]; % subject 3 right NBM % TODO make flexible readout based on subjects

%% load file data

global data;
% data_folder = 'sbj_3_duty_45_ntrl_10_amp_500000_targ_128-142-150_tpos_212-152-170';

% fileInfo = dir(fullfile('M:\Documents\repos\PRESTUS_forked\data\sims', data_folder, '*final_isppa_orig_coord.nii.gz'));
% data = niftiread(fullfile(fileInfo.folder, fileInfo.name));

%% load rest
head = niftiread('M:\Documents\scans\segmentation_results\m2m_sub-003\final_tissues.nii.gz');
% skull = head; skull(skull ~= 7) = 0;

% % Flatten the 3D matrix
% flattened_data = data(:);
% 
% % Sort the flattened values
% sorted_data = sort(flattened_data);
% 
% %%
% max_val = max(data, [], 'all');
% cutoff = sorted_data(floor(0.999*numel(sorted_data)));
% cutoff_y = 0.01*max_val;

%% Plot the sorted values
% figure;
% plot(sorted_data);
% title('Ordered Distribution of Flattened Values');
% xlabel('Index');
% ylabel('Value');

%% inspect voxels

% % set all values below cutoff to 0
% data_cut_off = data;
% data_cut_off(data_cut_off < cutoff) = 0;
% 
% volshow(data_cut_off);

%% setup figure

figure;
hold on;
lighting gouraud
axis equal
camlight('infinite');

%% show head
% 
% Ds = smooth3(double(head>0));
% skin_isosurface = isosurface(Ds,0.5);
% hiso = patch(skin_isosurface,...
%    'FaceColor',[1,.75,.65],...
%    'EdgeColor','none',...
%    'facealpha',0.4);

%% show skull/brain

% tricks to fill the cavities
within = ismember(head, [1 2 3]); % brain
within = ismember(head, [4 7 8]); % skull
% within = ismember(head, [1 2 3 4 7 8 9]); % skull
within = ~within;
within = imfill(within, 'holes');
within = ~within;
% make patch
Ds_skull = smooth3(double(within));%smooth3(double(head==7));
skin_isosurface = isosurface(Ds_skull,0.5);
hiso = patch(skin_isosurface,...
   'FaceColor',[0.5 0.5 0.5],...
   'EdgeColor','none',...
   'facealpha',0.05);

hiso.Tag = 'brain';

%% add original target

% add_3d_cross(pos_NBM([2 1 3]), 20);
target_patch = plot_sphere(5, pos_NBM([2 1 3]));

%% add sim results

% good cutoff percentage if not limited to brain: 0.9995
% good cutoff percentage if limited to brain (i.e. when 'HeadData' param is specified): 0.9999


% % z undershoot example + not enough depth?
% data = add_sim_result_patch( ...
%      'sbj_3_duty_45_ntrl_5_amp_500000_targ_128-142-156_tpos_212-144-172', ...
%      'AddTransducer', true, 'Color', colors.orangeish, 'CutoffPerc', 0.9999, 'PatchName', 'sim2', 'HeadData', head);
% 
% % z overshoot example + not enough depth?
% data = add_sim_result_patch( ...
%     'sbj_3_duty_45_ntrl_5_amp_500000_targ_128-142-160_tpos_209-146-180', ...
%     'AddTransducer', true, 'Color', colors.blueish, 'CutoffPerc', 0.9999, 'PatchName', 'sim1', 'HeadData', head);

% % same position, but with IGT transducers
% data = add_sim_result_patch( ...
%     'sbj_3_duty_45_ntrl_400_amp_105000_targ_128-142-156_tpos_212-144-172', ...
%     'AddTransducer', true, 'Color', colors.orangeish, 'CutoffPerc', 0.9999, 'PatchName', 'sim3', 'HeadData', head);
% 
% data = add_sim_result_patch( ...
%     'sbj_3_duty_45_ntrl_400_amp_105000_targ_128-142-160_tpos_209-146-180', ...
%     'AddTransducer', true, 'Color', colors.blueish, 'CutoffPerc', 0.9999, 'PatchName', 'sim4', 'HeadData', head);

% % higher pressure (should not matter) and deeper targeting
% data = add_sim_result_patch( ...
%     'sbj_3_duty_45_ntrl_400_amp_200000_targ_115-140-140_tpos_212-152-159', ...
%     'AddTransducer', true, 'Color', colors.orangeish, 'CutoffPerc', 0.9999, 'PatchName', 'sim3', 'HeadData', head);
% 
% data = add_sim_result_patch( ...
%     'sbj_3_duty_45_ntrl_400_amp_200000_targ_115-140-152_tpos_212-152-169', ...
%     'AddTransducer', true, 'Color', colors.blueish, 'CutoffPerc', 0.9999, 'PatchName', 'sim4', 'HeadData', head);
% 
% data = add_sim_result_patch( ...
%     'sbj_3_duty_45_ntrl_400_amp_200000_targ_115-140-146_tpos_212-152-169', ...
%     'AddTransducer', true, 'Color', colors.purplish, 'CutoffPerc', 0.9999, 'PatchName', 'sim4', 'HeadData', head);

colors = struct();
colors.blueish = [0 0.4470 0.7410];
colors.orangeish = [0.8500 0.3250 0.0980];
colors.purplish = [0.4940 0.1840 0.5560];
colors.greenish = [0.4660 0.6740 0.1880];

% again slightly deeper targeting, but adjusted the expected_focal_distance in this new run (still off!)

% % deepest focus and higher elevation
% data = add_sim_result_patch( ...
%     'sbj_3_duty_45_ntrl_5_amp_200000_targ_61-134-133_tpos_212-152-170', ...
%     'AddTransducer', true, 'Color', colors.blueish, 'CutoffPerc', 0.9999, 'PatchName', 'sim3', 'HeadData', head);
% 
% data = add_sim_result_patch( ...
%     'sbj_3_duty_45_ntrl_5_amp_200000_targ_109-140-145_tpos_195-145-219', ...
%     'AddTransducer', true, 'Color', colors.orangeish, 'CutoffPerc', 0.9999, 'PatchName', 'sim3', 'HeadData', head);

% % transducer positioned lower
% data = add_sim_result_patch( ...
%     'sbj_3_duty_45_ntrl_5_amp_200000_targ_128-142-150_tpos_212-152-150', ...
%     'AddTransducer', true, 'Color', colors.orangeish, 'CutoffPerc', 0.9999, 'PatchName', 'sim3', 'HeadData', head);

% % transducer positioned even lower
% data = add_sim_result_patch( ...
%     'sbj_3_duty_45_ntrl_5_amp_200000_targ_128-142-150_tpos_212-152-130', ...
%     'AddTransducer', true, 'Color', colors.orangeish, 'CutoffPerc', 0.9999, 'PatchName', 'sim3', 'HeadData', head);

% transducer positioned even lower
% add_sim_result_patch
data = add_sim_skull_heating( ...
    'sbj_3_duty_45_ntrl_400_targ_98-142-149_tpos_9-137-146_tempwin1_control', ...
    'AddTransducer', true, 'Color', colors.orangeish, 'CutoffPerc', 0.9995, 'PatchName', 'sim3', 'HeadData', head);

view(57, 17);

% data = add_sim_result_patch( ...
%     'sbj_3_duty_45_ntrl_400_amp_105000_targ_128-142-156_tpos_212-144-172', ...
%     'AddTransducer', true, 'Color', [0.9290 0.6940 0.1250], 'CutoffPerc', 0.9999, 'PatchName', 'sim3', 'HeadData', head);
% 
% data = add_sim_result_patch( ...
%     'sbj_3_duty_45_ntrl_400_amp_105000_targ_128-142-160_tpos_209-146-180', ...
%     'AddTransducer', true, 'Color', [0.9290 0.9 0.1250], 'CutoffPerc', 0.9999, 'PatchName', 'sim4', 'HeadData', head);