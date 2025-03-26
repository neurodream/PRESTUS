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
parameters = load_parameters('default_paths_windows.yaml');

for sbj_i = 1:8

    data = niftiread(fullfile(parameters.seg_path, '..', sprintf('sub-%03d_T1.nii.gz', sbj_i)));
    data_sorted = sort(data(:));
    data_sorted = data_sorted(data_sorted > 0);
    % plot(data_sorted);
    % ylim([0,50]);
    data_binary = data > prctile(data_sorted, 65);
    volshow(data_binary);

    disp('');

end

% data = niftiread(fullfile(parameters.seg_path, '..', 'sub-001_ses-mri01_T1w_transformed.nii.gz'));
% data_sorted = sort(data(:));
% data_sorted = data_sorted(data_sorted > 0);
% % plot(data_sorted);
% % ylim([0,50]);
% data_binary = data > prctile(data_sorted, 5);
% volshow(data);