%% (TODO) append to xlsx file; create new if not existing

% filename: tissue_based_postprocessing.xlsx in outpath
% for each nifti (temp, MECH INDEX, ...) its own tab
% row ID columns: sbj_ID, iteration, prefix

%%

% close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(currentFile);
cd(rootpath); % repos/PRESTUS_forked/scripts
cd ..

%% parameters
sbj_ID = 9; % careful: still hardcoded below (really? don't see it anymore)
iteration = 4;
prefix = 'pilot_titration';
filepath = sprintf('../../scans/sim_outputs/sub-%03d', sbj_ID);
outpath = '';
% filepath = 'p/2425076.01/piloting/titration_data';

%%
T = readtable('data/transducer_pos/position_LUT.xlsx');
target_coord = [T.x_r(T.sbj_ID == sbj_ID) T.y_r(T.sbj_ID == sbj_ID) T.z_r(T.sbj_ID == sbj_ID)];

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(currentFile);
cd(rootpath); % repos/PRESTUS_forked/
cd ..

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

% load parameters
load(sprintf(fullfile(filepath, 'sub-%03d_parametersL--r_R--r_%s_it%d_imprecisionnone.mat'), sbj_ID, prefix, iteration), 'parameters');

%% limit to brain
% make sure the subject ID match of seg_file and target:
segmentation_folder = fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', sbj_ID));
filename_segmented = fullfile(segmentation_folder, 'final_tissues.nii.gz');

layers = niftiread(filename_segmented);
layers_info = niftiinfo(filename_segmented);
% [layers, layers_info] = swapNiftiXY(layers, layers_info); % x and y seem swapped in nifti, need to be swapped back

head = layers > 0;

head = fill_head(head);

skull = layers == 7 | layers == 8;

% get the brain (assuming "head" exists globally)
within_brain = ismember(layers, [1 2 3]); % brain

% shrink the shape: conservative brain estimate, so to prevent
% "near field" in estimation of max pos
se = strel('sphere', 1); % A spherical structuring element with a radius of 1
within_brain = imerode(within_brain, se);

% within = ismember(head, [1 2 3 4 7 8 9]); % skull
within_brain = ~within_brain;
within_brain = imfill(within_brain, 'holes');
within_brain = ~within_brain;

%% input params (TODO)

disp('INPUT PARAMS')
disp('left transducer optimal distance');
disp(parameters.transducers(1).optim_params.focal_distance_mm);

%% acoustic postprocessing

disp('PRESSURE (MPa)');
% data = niftiread(sprintf('../../scans/sim_outputs/sub-%03d/sub-%03d_layered_final_pressureL--r_R--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, sbj_ID, prefix, iteration));
data = niftiread(sprintf(fullfile(filepath, 'sub-%03d_layered_final_pressureL--r_R--r_%s_it%d_imprecisionnone.nii.gz'), sbj_ID, prefix, iteration));
data = data / 1000000;
disp('global:');
disp(max(data(:)));
data_brain = data;
data_brain(~within_brain) = 0;
disp('within brain:');
disp(max(data_brain(:)));
disp('at target:');
disp(data(target_coord(1),target_coord(2),target_coord(3)));
disp('around target (max):');
target_cube = data( ...
    target_coord(1)-3:target_coord(1)+3, ...
    target_coord(2)-3:target_coord(2)+3, ...
    target_coord(3)-3:target_coord(3)+3 ...
    );
disp(max(target_cube(:)));
disp('around target (median):');
disp(median(target_cube(:)));
disp('');

disp('INTENSITY (W/cm²)');
% data = niftiread(sprintf('../../scans/sim_outputs/sub-%03d/sub-%03d_layered_final_intensityL--r_R--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, sbj_ID, prefix, iteration));
data = niftiread(sprintf(fullfile(filepath, 'sub-%03d_layered_final_intensityL--r_R--r_%s_it%d_imprecisionnone.nii.gz'), sbj_ID, prefix, iteration));
disp('global:');
disp(max(data(:)));
data_brain = data;
data_brain(~within_brain) = 0;
disp('within brain:');
disp(max(data_brain(:)));
disp('at target:');
disp(data(target_coord(1),target_coord(2),target_coord(3)));
disp('around target (max):');
target_cube = data( ...
    target_coord(1)-3:target_coord(1)+3, ...
    target_coord(2)-3:target_coord(2)+3, ...
    target_coord(3)-3:target_coord(3)+3 ...
    );
disp(max(target_cube(:)));
disp('around target (median):');
disp(median(target_cube(:)));
disp('');

disp('MECHANICAL INDEX');
% data = niftiread(sprintf('../../scans/sim_outputs/sub-%03d/sub-%03d_layered_final_mechanicalindexL--r_R--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, sbj_ID, prefix, iteration));
data = niftiread(sprintf(fullfile(filepath, 'sub-%03d_layered_final_mechanicalindexL--r_R--r_%s_it%d_imprecisionnone.nii.gz'), sbj_ID, prefix, iteration));
disp('global:');
disp(max(data(:)));
data_brain = data;
data_brain(~within_brain) = 0;
disp('within brain:');
disp(max(data_brain(:)));
disp('at target:');
disp(data(target_coord(1),target_coord(2),target_coord(3)));
disp('around target (max):');
target_cube = data( ...
    target_coord(1)-3:target_coord(1)+3, ...
    target_coord(2)-3:target_coord(2)+3, ...
    target_coord(3)-3:target_coord(3)+3 ...
    );
disp(max(target_cube(:)));
disp('around target (median):');
disp(median(target_cube(:)));
disp('');

%% heating

disp('CEM43');

data = niftiread(sprintf('../../scans/sim_outputs/sub-%03d/sub-%03d_final_CEM43L--r_R--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, sbj_ID, prefix, iteration));
data = niftiread(sprintf(fullfile(filepath, 'sub-%03d_layered_final_CEM43L--r_R--r_%s_it%d_imprecisionnone.nii.gz'), sbj_ID, prefix, iteration));
disp('global:');
disp(max(data(:)));
data_brain = data;
data_brain(~within_brain) = 0;
disp('within brain:');
disp(max(data_brain(:)));
% disp('at target:');
% disp(data(target_coord(1),target_coord(2),target_coord(3)));
% disp('around target (max):');
% target_cube = data( ...
%     target_coord(1)-3:target_coord(1)+3, ...
%     target_coord(2)-3:target_coord(2)+3, ...
%     target_coord(3)-3:target_coord(3)+3 ...
%     );
% disp(max(target_cube(:)));
% disp('around target (median):');
% disp(median(target_cube(:)));
% disp('');
data_skull = data(skull);
disp('within skull');
disp(max(data_skull(:)));

disp('temp');

% data = niftiread(sprintf('../../scans/sim_outputs/sub-%03d/sub-%03d_final_tempL--r_R--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, sbj_ID, prefix, iteration));
data = niftiread(sprintf(fullfile(filepath, 'sub-%03d_layered_final_tempL--r_R--r_%s_it%d_imprecisionnone.nii.gz'), sbj_ID, prefix, iteration));
disp('global:');
disp(max(data(:)));
data_brain = data;
data_brain(~within_brain) = 0;
disp('within brain:');
disp(max(data_brain(:)));
% disp('at target:');
% disp(data(target_coord(1),target_coord(2),target_coord(3)));
% disp('around target (max):');
% target_cube = data( ...
%     target_coord(1)-3:target_coord(1)+3, ...
%     target_coord(2)-3:target_coord(2)+3, ...
%     target_coord(3)-3:target_coord(3)+3 ...
%     );
% disp(max(target_cube(:)));
% disp('around target (median):');
% disp(median(target_cube(:)));
% disp('');
data_skull = data(skull);
disp('within skull');
disp(max(data_skull(:)));