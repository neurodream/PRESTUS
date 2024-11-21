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
%   3. sham or not (true or false), 
%    4. how many z voxels for under/overshoot correction,
%     5. target intensity,
%      6. a string abbreviation which describes more precisely (e.g. the number of iterations, e.g. 'it1'; can also be left blank: '')
%       7. and if errors in spatial precision are to be simulated as well ('none' for no imprecision, 'transducer', 'target', or 'both')

% example:
% NBM_run(1, 'R', false, 4, 100, 'it1_FW', 'none');
% NBM_run(1, 'R', true,  4, 100, 'it1_FW', 'none');
% NBM_run(2, 'R', false, 4, 100, 'it1_FW', 'none');
% NBM_run(2, 'R', true,  4, 100, 'it1_FW', 'none');
% NBM_run(3, 'R', false, 3, 100, 'it1_FW', 'none');
% NBM_run(3, 'R', true,  3, 100, 'it1_FW', 'none');
% NBM_run(4, 'R', false, 2, 100, 'it1_FW', 'none');
% NBM_run(4, 'R', true,  2, 100, 'it1_FW', 'none');
% NBM_run(5, 'R', false, 2, 100, 'it1_FW', 'none');
% NBM_run(5, 'R', true,  2, 100, 'it1_FW', 'none');
NBM_run(6, 'R', false, 0, 100, 'it1_heatingtimeline', 'none');
% NBM_run(6, 'R', true,  0, 100, 'it1_FW', 'none');

% empirical z corrections for the subjects:
% subject 1: 4
% subject 2: 4
% subject 3: 3
% subject 4: 2
% subject 5: 2
% subject 6: 0

% TODO add to single_subject_pipeline, but careful: assumes a stored nifti
% file
% plot_transducer_pos(parameters, sbj_ID, plot_scalp, plot_skull, plot_intensity, save)