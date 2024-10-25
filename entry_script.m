close all; clear; clc;

cd /home/sleep/nicade/Documents/repos/PRESTUS_Julian/;

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

%%

parameters = load_parameters('nico_test_double_acoustic_100mm_config.yaml');

subject_id = 5;
ID = '_L+y--l_R+y--r_';

parameters.results_filename_affix = ID;

parameters.transducers(1).name = 'L';
parameters = calculate_transducer_phases(parameters, 1, 63.7943, 15, 100);
parameters = get_transducer_pos(parameters, subject_id, 'l', 1, [0.2 -1 0], [0 0 0], [0 0 0]);
parameters.transducers(2).name = 'R';
parameters = calculate_transducer_phases(parameters, 2, 64.691, 15, 100);
parameters = get_transducer_pos(parameters, subject_id, 'r', 2, [0.2 1 0], [0 0 0], [0 0 0]);

update_transducers_and_run(subject_id, parameters, [ID 'var1'], 'transducer');
% update_transducers_and_run(subject_id, parameters, [ID 'var2'], 'transducer');
% update_transducers_and_run(subject_id, parameters, [ID 'var3'], 'focus');
% update_transducers_and_run(subject_id, parameters, [ID 'var4'], 'focus');
% update_transducers_and_run(subject_id, parameters, [ID 'var5'], 'both');
% update_transducers_and_run(subject_id, parameters, [ID 'var6'], 'both');








%% old

% f1 = parameters.transducers(1).focus_pos_t1_grid;
% f2 = parameters.transducers(2).focus_pos_t1_grid;
% t1 = parameters.transducers(1).pos_t1_grid;
% t2 = parameters.transducers(2).pos_t1_grid;
% 
% parameters.results_filename_affix = [ID 'var1'];
% vt1 = [0, -1, -2]; vf1 = [0, 0, 0]; vt2 = [3, -2, -1]; vf2 = [0, 0, 0];
% parameters.transducers(1).pos_t1_grid 		= [t1(1) + vt1(1),t1(2) + vt1(2),t1(3) + vt1(3)];
% parameters.transducers(1).focus_pos_t1_grid 	= [f1(1) + vf1(1),f1(2) + vf1(2),f1(3) + vf1(3)];
% parameters.transducers(2).pos_t1_grid 		= [t2(1) + vt2(1),t2(2) + vt2(2),t2(3) + vt2(3)];
% parameters.transducers(2).focus_pos_t1_grid 	= [f2(1) + vf2(1),f2(2) + vf2(1),f2(3) + vf2(1)];
% single_subject_pipeline_with_slurm(subject_id, parameters);
% 
% parameters.results_filename_affix = [ID 'var2'];
% vt1 = [0, 2, 0]; vf1 = [0, 0, 0]; vt2 = [2, 0, 1]; vf2 = [0, 0, 0];
% parameters.transducers(1).pos_t1_grid 		= [t1(1) + vt1(1),t1(2) + vt1(2),t1(3) + vt1(3)];
% parameters.transducers(1).focus_pos_t1_grid 	= [f1(1) + vf1(1),f1(2) + vf1(2),f1(3) + vf1(3)];
% parameters.transducers(2).pos_t1_grid 		= [t2(1) + vt2(1),t2(2) + vt2(2),t2(3) + vt2(3)];
% parameters.transducers(2).focus_pos_t1_grid 	= [f2(1) + vf2(1),f2(2) + vf2(1),f2(3) + vf2(1)];
% single_subject_pipeline_with_slurm(subject_id, parameters);
% 
% parameters.results_filename_affix = [ID 'var3'];
% vt1 = [0, 0, 0]; vf1 = [-1, 3, 0]; vt2 = [0, 0, 0]; vf2 = [-3, 2, 2];
% parameters.transducers(1).pos_t1_grid 		= [t1(1) + vt1(1),t1(2) + vt1(2),t1(3) + vt1(3)];
% parameters.transducers(1).focus_pos_t1_grid 	= [f1(1) + vf1(1),f1(2) + vf1(2),f1(3) + vf1(3)];
% parameters.transducers(2).pos_t1_grid 		= [t2(1) + vt2(1),t2(2) + vt2(2),t2(3) + vt2(3)];
% parameters.transducers(2).focus_pos_t1_grid 	= [f2(1) + vf2(1),f2(2) + vf2(1),f2(3) + vf2(1)];
% single_subject_pipeline_with_slurm(subject_id, parameters);
% 
% parameters.results_filename_affix = [ID 'var4'];
% vt1 = [0, 0, 0]; vf1 = [1, -3, 1]; vt2 = [0, 0, 0]; vf2 = [-3, -3, 2];
% parameters.transducers(1).pos_t1_grid 		= [t1(1) + vt1(1),t1(2) + vt1(2),t1(3) + vt1(3)];
% parameters.transducers(1).focus_pos_t1_grid 	= [f1(1) + vf1(1),f1(2) + vf1(2),f1(3) + vf1(3)];
% parameters.transducers(2).pos_t1_grid 		= [t2(1) + vt2(1),t2(2) + vt2(2),t2(3) + vt2(3)];
% parameters.transducers(2).focus_pos_t1_grid 	= [f2(1) + vf2(1),f2(2) + vf2(1),f2(3) + vf2(1)];
% single_subject_pipeline_with_slurm(subject_id, parameters);
% 
% parameters.results_filename_affix = [ID 'var5'];
% vt1 = [3, -3, -3]; vf1 = [-1, 0, -1]; vt2 = [-1, 3, 1]; vf2 = [-1, 1, 1];
% parameters.transducers(1).pos_t1_grid 		= [t1(1) + vt1(1),t1(2) + vt1(2),t1(3) + vt1(3)];
% parameters.transducers(1).focus_pos_t1_grid 	= [f1(1) + vf1(1),f1(2) + vf1(2),f1(3) + vf1(3)];
% parameters.transducers(2).pos_t1_grid 		= [t2(1) + vt2(1),t2(2) + vt2(2),t2(3) + vt2(3)];
% parameters.transducers(2).focus_pos_t1_grid 	= [f2(1) + vf2(1),f2(2) + vf2(1),f2(3) + vf2(1)];
% single_subject_pipeline_with_slurm(subject_id, parameters);
% 
% parameters.results_filename_affix = [ID 'var6'];
% vt1 = [-3, -1, -1]; vf1 = [0, 0, 3]; vt2 = [-2, 2, 0]; vf2 = [0, -1, -1];
% parameters.transducers(1).pos_t1_grid 		= [t1(1) + vt1(1),t1(2) + vt1(2),t1(3) + vt1(3)];
% parameters.transducers(1).focus_pos_t1_grid 	= [f1(1) + vf1(1),f1(2) + vf1(2),f1(3) + vf1(3)];
% parameters.transducers(2).pos_t1_grid 		= [t2(1) + vt2(1),t2(2) + vt2(2),t2(3) + vt2(3)];
% parameters.transducers(2).focus_pos_t1_grid 	= [f2(1) + vf2(1),f2(2) + vf2(1),f2(3) + vf2(1)];
% single_subject_pipeline_with_slurm(subject_id, parameters);