close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(currentFile);
cd(rootpath); % repos/PRESTUS_forked/
cd ..

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

% load the scans
parameters = load_parameters();

for sbj_i = 1:7

    data = niftiread(fullfile(parameters.seg_path, '..', sprintf('sub-%03d_T1.nii.gz', sbj_i)));
    data_sorted = sort(data(:));
    data_sorted = data_sorted(data_sorted > 0);
    % plot(data_sorted);
    % ylim([0,50]);
    data_binary = data > prctile(data_sorted, 65);
    volshow(data_binary);

    disp('');

end