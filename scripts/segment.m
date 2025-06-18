close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
currentpath = fileparts(currentFile);
cd(currentpath);
cd ..

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

parameters = load_parameters();
parameters.ld_library_path = '/opt/gcc/7.2.0/lib64';

for sub_ind = 102%101:112

    run_segmentation('../../scans', sub_ind, sprintf('sub-%d_T1.nii.gz', sub_ind), [], parameters);

end

disp('')