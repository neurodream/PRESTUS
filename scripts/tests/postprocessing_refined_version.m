close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(currentFile);
cd(rootpath); % repos/PRESTUS_forked/
cd ..;
cd ..;

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

load('/home/sleep/nicade/Documents/scans/sim_outputs/sub-008/sub-008_parametersL+z-r_R+z-r_it9_just_test_imprecisionnone_100125_1231.mat');
load('/home/sleep/nicade/Documents/scans/sim_outputs/sub-008/sub-008_layered_resultsL+z-r_R+z-r_it9_just_test_imprecisionnone.mat')

output_pressure_file = 'test.xlsx';

data = [];

% What is the highest pressure level for every gridpoint
data.pressure = gather(sensor_data.p_max_all); % gather is used since it could be a GPU array

% Calculates the Isppa for every gridpoint
data.intensity = data.pressure.^2./(2*(kwave_medium.sound_speed.*kwave_medium.density)).*1e-4;

% Calculates the Mechanical Index for every gridpoint
% TODO figure out how to implement different source frequencies if needed
data.mechanicalindex = (data.pressure/10^6)/sqrt((parameters.transducers(1).source_freq_hz/10^6));

[medium_masks, segmented_image_cropped, skull_edge, trans_pos_final, ...
            focus_pos_final, t1_image_orig, t1_header, final_transformation_matrix, ...
            inv_final_transformation_matrix] = preprocess_brain(parameters, 8, 1);

postprocessing_quantification_acoustics(8, parameters, medium_masks, output_pressure_file, data);