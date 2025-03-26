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

% % possible configs:
% nico_test_double_acoustic_same_temp0_config
% acoustic_sim
% heating_sim
% no_sim_just_setup
% IS_double_75mm
% IS_double_100mm
% IS_L_100mm_R_75mm
% debug


% NBM_run( 5, 'R', false, 2, 60, 'checkOffpeaks_it7',       {'nico_test_double_acoustic_same_temp0_config', 'IS_double_100mm', 'acoustic_sim'}, 'none');
% NBM_run( 8, 'R', false, 1, 100, 'fixMI_it4',       {'nico_test_double_acoustic_same_temp0_config', 'IS_double_100mm', 'acoustic_sim'}, 'none');
% NBM_run(10, 'R', false, 4, 100, 'fixMI_it4',       {'nico_test_double_acoustic_same_temp0_config', 'IS_double_100mm', 'acoustic_sim'}, 'none');


NBM_run( 1, 'R', true,  4, 60, 'fixMI_it11',       {'nico_test_double_acoustic_same_temp0_config', 'IS_double_100mm', 'acoustic_sim'}, 'none');
NBM_run( 1, 'R', false, 4, 60, 'fixMI_it11',       {'nico_test_double_acoustic_same_temp0_config', 'IS_double_100mm', 'acoustic_sim'}, 'none');



% empirical z corrections for the subjects: (TODO: make this programmatically)
% subject 1:  4
% subject 2:  4
% subject 3:  3
% subject 4:  2
% subject 5:  2
% subject 6:  0
% subject 7:  4
% subject 8:  1
% subject 9:  - (error because no T2) 
% subject 10: 4

% TODO add to single_subject_pipeline, but careful: assumes a stored nifti
% file
% plot_transducer_pos(parameters, sbj_ID, plot_scalp, plot_skull, plot_intensity, save)