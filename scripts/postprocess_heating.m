% assumes a mat file

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

% load('sub-001_layered_heating_resL+z--l_R+z--l_heating_check')

sbj_ID = 1;
% affix = 'L--r_R--r_heating_check_400trials';
% affix = 'L+z--l_R+z--l_heating_check_400trials';
% affix = 'L+z--r_R+z--r_heating_check_400trials';
affix = 'L+z--r_R+z--r_heating_check_400trials_same_temp0_DC50';

data_path = fullfile(parameters.data_path, 'sim_outputs', sprintf('sub-%03d', sbj_ID));
filename_heating = sprintf('sub-%03d_layered_heating_res%s', sbj_ID, affix);
filename_parameters = dir(fullfile(data_path, sprintf('sub-%03d_parameters%s*', sbj_ID, affix)));
filename_parameters = filename_parameters.name;

% % once it did not store the filename, therefore loading it manually here...
% % TODO figure out why not stored
% affix_params = 'L--l_R--l_test_parallel_no_controlling_undershoot';
% filename_parameters = dir(fullfile(data_path, sprintf('sub-%03d_parameters_%s*', sbj_ID, affix_params)));
% filename_parameters = filename_parameters.name;

load(fullfile(data_path, filename_heating));
load(fullfile(data_path, filename_parameters));
parameters.usepseudoCT = 0;
parameters.data_path = '/home/sleep/nicade/Documents/scans/';
parameters.seg_path = '/home/sleep/nicade/Documents/scans/segmentation_results/';
parameters.temp_output_dir = fullfile(parameters.data_path, 'sim_outputs/');
parameters.output_dir = fullfile(parameters.temp_output_dir, sprintf('sub-%03d', sbj_ID));
parameters.debug_dir = fullfile(parameters.output_dir, 'debug');
parameters.subject_id = sbj_ID; % TODO interesting that this is not found

% convert to normal world coordinates & store nifti
% TODO: add this to single_subject_pipeline; store instead of mat

fname_out_heating = sprintf('sub-%03d_%s_heatrise%s', sbj_ID, 'layered', affix);

% if 

[medium_masks, ~, ~, ~, ~, ~, t1_header, final_transformation_matrix, ...
    inv_final_transformation_matrix] = preprocess_brain(parameters, sbj_ID, 1);

parameters.grid_dims = size(medium_masks);
kwave_medium = setup_medium(parameters, medium_masks);

heatrise = gather(maxT - kwave_medium.temp_0); % parameters.thermal.temp_0.brain);

data_heating = tformarray(heatrise, inv_final_transformation_matrix, makeresampler('cubic', 'fill'), [1 2 3], [1 2 3], t1_header.ImageSize, [], 0) ;

t1_header.Datatype = 'single';
niftiwrite(data_heating, fullfile(data_path, fname_out_heating), t1_header, 'Compressed', true);

% missing variables 't1_header', 'inv_final_transformation_matrix'