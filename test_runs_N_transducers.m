close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(currentFile);
cd(rootpath); % repos/PRESTUS_forked/

load test_params2
parameters.sim_path = parameters.temp_output_dir;

addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

% single_subject_pipeline(7, parameters);

single_subject_pipeline_with_slurm(7, parameters)