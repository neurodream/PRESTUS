% close all;
clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(fileparts(currentFile));
cd(rootpath); % repos/PRESTUS_forked/

addpath('functions');
addpath(genpath('toolboxes'));

% base config ("hard" params)
parameters = load_parameters('nico_test_double_acoustic_100mm_same_temp0_config.yaml');

%% 
% load in the parameter file

cd(fullfile(parameters.data_path, 'sim_outputs'));

subj_ID = 1;
ID_sham   = 'L+z-r_R+z-r_it1_heatingtimeline_imprecisionnone';
ID_active = 'L+z--r_R+z--r_it1_heatingtimeline_imprecisionnone';

files = dir([sprintf('sub-%03d/sub-%03d', subj_ID, subj_ID) '_parameters' ID_sham '*']);
file = files(1);
load(fullfile(file.folder, file.name));
parameters_sham = get_simulated_axial_intensity(parameters);

files = dir([sprintf('sub-%03d/sub-%03d', subj_ID, subj_ID) '_parameters' ID_active '*']);
file = files(1);
load(fullfile(file.folder, file.name));
parameters_active = get_simulated_axial_intensity(parameters);

cd(rootpath);

fig1 = figure('Position', [100, 100, 350, 100]); hold on;
plot(parameters_active.transducers(1).axial_position, parameters_active.transducers(1).axial_intensity_sim_FW, 'LineWidth', 2.5);
plot(parameters_sham.transducers(1).axial_position, parameters_sham.transducers(1).axial_intensity_sim_FW, 'LineWidth', 2.5);
xlim([0 150]);
ylim([0 115]);
legend('active', 'sham');
set(gcf, 'Color', 'w');
exportgraphics(gcf, 'data/transducerL_FW_axial_sbj1.png', 'BackgroundColor', 'white');

fig2 = figure('Position', [100, 100, 350, 100]); hold on;
plot(parameters_active.transducers(2).axial_position, parameters_active.transducers(2).axial_intensity_sim_FW, 'LineWidth', 2.5);
plot(parameters_sham.transducers(2).axial_position, parameters_sham.transducers(2).axial_intensity_sim_FW, 'LineWidth', 2.5);
xlim([0 150]);
ylim([0 115]);
legend('active', 'sham');
set(gcf, 'Color', 'w');
exportgraphics(gcf, 'data/transducerR_FW_axial_sbj1.png', 'BackgroundColor', 'white');

% legend(parameters.transducers(1).name, parameters.transducers(2).name);

% cd ..
% cd ..
% cd scans/sim_outputs/sub-003/debug/
% load('sub-003_parameters_L+y+z--l+z_R+y+z--l+z__251024_2055.mat')
% parameters2 = parameters;