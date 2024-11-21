close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(currentFile);
cd(rootpath); % repos/PRESTUS_forked/
cd ..

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

% base config ("hard" params)
parameters = load_parameters('nico_test_double_acoustic_100mm_config.yaml');

sub_id = 3;

% load the MRI and segment
nifti_path = fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', sub_id), 'final_tissues.nii.gz');

layers = niftiread(nifti_path);
info   = niftiinfo(nifti_path);

% simnibs-segmented tissues
skull = layers == 7 | layers == 8;
brain = layers == 1 | layers == 2;
skin  = layers == 5;

% inner brain
SE = strel('sphere', 6);

skull_dilated = imdilate(skull, SE);
brain_eroded  = imerode(brain, SE); % not optimal: erodes around the gyri

brain_inner = brain & ~skull_dilated;

% % targets (careful: hardcoded for subject 1)
% targetL = [96 144 149];
% targetR = [130 144 150];
T = readtable('data/transducer_pos/position_LUT.xlsx');
targetL = [T.y_l(T.sbj_ID == sub_id) T.x_l(T.sbj_ID == sub_id) T.z_l(T.sbj_ID == sub_id)];
targetR = [T.y_r(T.sbj_ID == sub_id) T.x_r(T.sbj_ID == sub_id) T.z_r(T.sbj_ID == sub_id)];

r = round(5/mean(info.PixelDimensions));
% Equation of the sphere: (x - px)^2 + (y - py)^2 + (z - pz)^2 <= r^2
[x, y, z] = ndgrid(1:size(layers,1), 1:size(layers,2), 1:size(layers,3));
ROItarget_L = (x - targetL(1)).^2 + (y - targetL(2)).^2 + (z - targetL(3)).^2 <= r^2;
ROItarget_R = (x - targetR(1)).^2 + (y - targetR(2)).^2 + (z - targetR(3)).^2 <= r^2;

% measure = 'intensity'; % mechanicalindex pressure intensity
% affix = 'it1_FW_imprecisionnone';
targeting = '--'; % -- -

% load the data
filename_pressure = sprintf('sub-%03d/sub-%03d_final_%sL+z%sr_R+z%sr_it1_%s_imprecisionnone.nii.gz', sub_id, sub_id, 'pressure', targeting, targeting, 'heatingtimeline');
pressure = niftiread(fullfile(parameters.data_path, 'sim_outputs', filename_pressure));
pressure = pressure/1000000; % in MegaPascal
filename_pressure_FW = sprintf('sub-%03d/sub-%03d_final_%sL+z%sr_R+z%sr_it1_%s_imprecisionnone.nii.gz', sub_id, sub_id, 'pressure', targeting, targeting, 'FW');
pressure_FW = niftiread(fullfile(parameters.data_path, 'sim_outputs', filename_pressure_FW));
pressure_FW = pressure_FW/1000000; % in MegaPascal

filename_intensity = sprintf('sub-%03d/sub-%03d_final_%sL+z%sr_R+z%sr_it1_%s_imprecisionnone.nii.gz', sub_id, sub_id, 'intensity', targeting, targeting, 'heatingtimeline');
intensity = niftiread(fullfile(parameters.data_path, 'sim_outputs', filename_intensity));
filename_intensity_FW = sprintf('sub-%03d/sub-%03d_final_%sL+z%sr_R+z%sr_it1_%s_imprecisionnone.nii.gz', sub_id, sub_id, 'intensity', targeting, targeting, 'FW');
intensity_FW = niftiread(fullfile(parameters.data_path, 'sim_outputs', filename_intensity_FW));

filename_MI = sprintf('sub-%03d/sub-%03d_final_%sL+z%sr_R+z%sr_it1_%s_imprecisionnone.nii.gz', sub_id, sub_id, 'mechanicalindex', targeting, targeting, 'heatingtimeline');
MI = niftiread(fullfile(parameters.data_path, 'sim_outputs', filename_MI));
filename_MI_FW = sprintf('sub-%03d/sub-%03d_final_%sL+z%sr_R+z%sr_it1_%s_imprecisionnone.nii.gz', sub_id, sub_id, 'mechanicalindex', targeting, targeting, 'FW');
MI_FW = niftiread(fullfile(parameters.data_path, 'sim_outputs', filename_MI_FW));

filename_TR = sprintf('sub-%03d/sub-%03d_layered_%sL+z%sr_R+z%sr_it1_%s_imprecisionnone.nii.gz', sub_id, sub_id, 'heating', targeting, targeting, 'heatingtimeline');
TR = niftiread(fullfile(parameters.data_path, 'sim_outputs', filename_TR));

filename_TD = sprintf('sub-%03d/sub-%03d_layered_%sL+z%sr_R+z%sr_it1_%s_imprecisionnone.nii.gz', sub_id, sub_id, 'maxCEM43', targeting, targeting, 'heatingtimeline');
TD = niftiread(fullfile(parameters.data_path, 'sim_outputs', filename_TD));

% max pressure free water
disp(['max pressure free water: ' num2str(max(pressure_FW(:)))]);
% max pressure brain
disp(['max pressure brain: ' num2str(max(pressure(brain)))]);
% max pressure scalp
% max pressure target
% max pressure target vs. off-target

% avg pressure target


% intensity max

% MI max

% Temperature rise max

% 