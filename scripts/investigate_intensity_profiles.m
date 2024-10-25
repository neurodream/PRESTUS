close all; clear; clc;

cd /home/sleep/nicade/Documents/repos/PRESTUS_forked/; % repos/PRESTUS_forked/

addpath('/home/sleep/nicade/Documents/repos/PRESTUS_forked/functions');
addpath(genpath('/home/sleep/nicade/Documents/repos/PRESTUS_forked/toolboxes'));

%% 
% load in the parameter file

cd /home/sleep/nicade/Documents/scans/sim_outputs/;

subj_ID = 1;
ID = '_L--r_R--r_2_';
% ID = '_L--r_R--r_var1';
% ID = 'nearparallel_var4_';

files = dir([sprintf('sub-%03d/sub-%03d', subj_ID, subj_ID) '_parameters' ID '*']);
file = files(1);
load(fullfile(file.folder, file.name));

parameters = get_simulated_axial_intensity(parameters);

figure; hold on;
plot(parameters.transducers(1).axial_position, parameters.transducers(1).axial_intensity_sim_FW);
plot(parameters.transducers(2).axial_position, parameters.transducers(2).axial_intensity_sim_FW);
legend(parameters.transducers(1).name, parameters.transducers(2).name);