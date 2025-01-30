close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(currentFile);
cd(rootpath); % repos/PRESTUS_forked/

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

% The parameters of the individual calls are, from left to right:
% 1. sbj_num, 
%  2. focus_side (L or R for unilateral, LR for bilateral), 
%   3. sham or not (true or false), n
%    4. how many z voxels for under/overshoot correction,
%     5. target intensity,
%      6. a string abbreviation which describes more precisely (e.g. the number of iterations, e.g. 'it1'; can also be left blank: '')
%       7. parameters filename (yaml file, but do not add the file extension here)
%        8. and if errors in spatial precision are to be simulated as well ('none' for no imprecision, 'transducer', 'target', or 'both')

% no sham - active
% NBM_run(9, 'R', false,  3, 100, 'intensity_check',       {'nico_test_double_acoustic_same_temp0_config', 'IS_L_100mm_R_75mm', 'no_sim_just_setup'}, 'none');
% sham - active
% NBM_run(8, 'R', false,  3, 100, 'new_sham_test_double100mm',       {'nico_test_double_acoustic_same_temp0_config', 'IS_double_100mm', 'acoustic_sim'}, 'none');
% NBM_run(8, 'R', true,  3, 100, 'new_sham_test_double100mm',       {'nico_test_double_acoustic_same_temp0_config', 'IS_double_100mm', 'acoustic_sim'}, 'none');

NBM_run(4, 'R', false,  2, 100, 'it15_heatingtimeline',       {'nico_test_double_acoustic_same_temp0_config', 'IS_double_100mm', 'heating_sim'}, 'none');
NBM_run(5, 'R', false,  2, 100, 'it15_heatingtimeline',       {'nico_test_double_acoustic_same_temp0_config', 'IS_double_100mm', 'heating_sim'}, 'none');
NBM_run(6, 'R', false,  0, 100, 'it15_heatingtimeline',       {'nico_test_double_acoustic_same_temp0_config', 'IS_double_100mm', 'heating_sim'}, 'none');
NBM_run(4, 'R', true,  2, 100, 'it15_heatingtimeline',       {'nico_test_double_acoustic_same_temp0_config', 'IS_double_100mm', 'heating_sim'}, 'none');
NBM_run(5, 'R', true,  2, 100, 'it15_heatingtimeline',       {'nico_test_double_acoustic_same_temp0_config', 'IS_double_100mm', 'heating_sim'}, 'none');
NBM_run(6, 'R', true,  0, 100, 'it15_heatingtimeline',       {'nico_test_double_acoustic_same_temp0_config', 'IS_double_100mm', 'heating_sim'}, 'none');

% empirical z corrections for the subjects: (TODO: make this programmatically)
% subject 1: 4
% subject 2: 4
% subject 3: 3
% subject 4: 2
% subject 5: 2
% subject 6: 0

% TODO add to single_subject_pipeline, but careful: assumes a stored nifti
% file
% plot_transducer_pos(parameters, sbj_ID, plot_scalp, plot_skull, plot_intensity, save)