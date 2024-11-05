close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(fileparts(currentFile));
cd(rootpath); % repos/PRESTUS_forked/

addpath('functions');
addpath(genpath('toolboxes'));

%% 
% load in the parameter file

cd(fullfile(parameters.data_path, 'sim_outputs'));

subj_ID = 3;
ID = '_L+y-z--l-z_R+y-z--l-z_test_intensity';
% ID = '_L--r_R--r_var1';
% ID = 'nearparallel_var4_';

files = dir([sprintf('sub-%03d/sub-%03d', subj_ID, subj_ID) '_parameters' ID '*']);
file = files(1);
load(fullfile(file.folder, file.name));

parameters = get_simulated_axial_intensity(parameters);

figure; hold on;
plot(parameters.transducers(1).axial_position, parameters.transducers(1).axial_intensity_sim_FW, 'LineWidth', 2.5);
plot(parameters.transducers(2).axial_position, parameters.transducers(2).axial_intensity_sim_FW, 'LineWidth', 2.5);
legend(parameters.transducers(1).name, parameters.transducers(2).name);

% cd ..
% cd ..
% cd scans/sim_outputs/sub-003/debug/
% load('sub-003_parameters_L+y+z--l+z_R+y+z--l+z__251024_2055.mat')
% parameters2 = parameters;