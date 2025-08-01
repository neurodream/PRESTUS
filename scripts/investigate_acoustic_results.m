% tissue_based_postprocessing_to_xlsx.m
% One row per run, one sheet per metric
% Creates tissue_based_postprocessing.xlsx if missing, otherwise appends

%% housekeeping
close all; clear; clc

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath    = fileparts(currentFile);
cd(rootpath)                 % .../repos/PRESTUS_forked/scripts
cd ..                        % .../repos/PRESTUS_forked

%% user params
sbj_ID    = 9;
iteration = 5;
prefix    = 'pilot_titration';
filepath  = sprintf('../../scans/sim_outputs/sub-%03d', sbj_ID);
% filepath  = '/project/2425076.01/piloting/titration_data/';
outpath   = '/project/2425076.01/piloting';
% outpath   = 'P:/2425076.01/piloting';
xlsx_name = 'tissue_based_postprocessing.xlsx';

segmentation_folder = sprintf('/home/sleep/nicade/Documents/scans/segmentation_results/m2m_sub-%03d', sbj_ID);
% segmentation_folder = sprintf('M:/Documents/scans/segmentation_results/m2m_sub-%03d', sbj_ID);

if sbj_ID == 7, add_z = '+z'; else, add_z = ''; end

if isempty(outpath), xlsx_file = fullfile(pwd, xlsx_name); else, xlsx_file = fullfile(outpath, xlsx_name); end

%% lookup target coordinates
T = readtable('data/transducer_pos/position_LUT.xlsx');
target_coord = [T.x_r(T.sbj_ID == sbj_ID) T.y_r(T.sbj_ID == sbj_ID) T.z_r(T.sbj_ID == sbj_ID)];
target_coord = round(target_coord);  % just in case

%% paths
addpath('functions')
addpath(genpath('toolboxes'))
addpath('/home/common/matlab/fieldtrip/qsub')  % uncomment for Donders HPC

%% load sim parameters
fname      = sprintf('sub-%03d_parametersL%s--r_R%s--r_%s_it%d_imprecisionnone.mat', sbj_ID, add_z, add_z, prefix, iteration);
param_file = fullfile(filepath, fname);
% param_file = sprintf(fullfile(filepath, 'sub-%03d_parametersL+z--r_R+z--r_%s_it%d_imprecisionnone.mat'), sbj_ID, prefix, iteration);
load(param_file, 'parameters')

%% masks (brain, skull)
% segmentation_folder = fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', sbj_ID));
filename_segmented  = fullfile(segmentation_folder, 'final_tissues.nii.gz');
layers      = niftiread(filename_segmented);
layers_info = niftiinfo(filename_segmented);

head  = layers > 0;
head  = fill_head(head);          % user function assumed on path

skull = layers == 7 | layers == 8;

within_brain = ismember(layers, [1 2 3]);      % brain
se = strel('sphere', 1);
within_brain = imerode(within_brain, se);
within_brain = ~within_brain;
within_brain = imfill(within_brain, 'holes');
within_brain = ~within_brain;

%% helper for target neighborhood
cube_rad = 3;
cx = target_coord(1); cy = target_coord(2); cz = target_coord(3);
target_cube_idx = {cx-cube_rad:cx+cube_rad, cy-cube_rad:cy+cube_rad, cz-cube_rad:cz+cube_rad};

%% PRESSURE (MPa)
fname      = sprintf('sub-%03d_layered_final_pressureL%s--r_R%s--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, add_z, add_z, prefix, iteration);
data = niftiread(fullfile(filepath, fname));
data = data / 1e6;

rec = makeRec(sbj_ID, iteration, prefix, data, within_brain, target_coord, target_cube_idx);
appendRow(xlsx_file, 'PRESSURE_MPa', rec)

%% INTENSITY (W_per_cm2)
% data = niftiread(fullfile(filepath, sprintf('sub-%03d_layered_final_intensityL+z--r_R+z--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, prefix, iteration)));
fname      = sprintf('sub-%03d_layered_final_intensityL%s--r_R%s--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, add_z, add_z, prefix, iteration);
data = niftiread(fullfile(filepath, fname));
rec = makeRec(sbj_ID, iteration, prefix, data, within_brain, target_coord, target_cube_idx);
appendRow(xlsx_file, 'INTENSITY_W_cm2', rec)

%% MECHANICAL_INDEX
% data = niftiread(fullfile(filepath, sprintf('sub-%03d_layered_final_mechanicalindexL+z--r_R+z--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, prefix, iteration)));
fname      = sprintf('sub-%03d_layered_final_mechanicalindexL%s--r_R%s--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, add_z, add_z, prefix, iteration);
data = niftiread(fullfile(filepath, fname));
rec = makeRec(sbj_ID, iteration, prefix, data, within_brain, target_coord, target_cube_idx);
appendRow(xlsx_file, 'MECHANICAL_INDEX', rec)

%% CEM43
% data = niftiread(fullfile(filepath, sprintf('sub-%03d_final_CEM43L+z--r_R+z--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, prefix, iteration)));
fname      = sprintf('sub-%03d_final_CEM43L%s--r_R%s--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, add_z, add_z, prefix, iteration);
data = niftiread(fullfile(filepath, fname));
rec = makeRec(sbj_ID, iteration, prefix, data, within_brain, target_coord, target_cube_idx);
rec.skull_max = max(data(skull));
appendRow(xlsx_file, 'CEM43', rec)

%% TEMP (degC)
% data = niftiread(fullfile(filepath, sprintf('sub-%03d_final_tempL+z--r_R+z--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, prefix, iteration)));
fname      = sprintf('sub-%03d_final_tempL%s--r_R%s--r_%s_it%d_imprecisionnone.nii.gz', sbj_ID, add_z, add_z, prefix, iteration);
data = niftiread(fullfile(filepath, fname));
rec = makeRec(sbj_ID, iteration, prefix, data, within_brain, target_coord, target_cube_idx);
rec.skull_max = max(data(skull));
appendRow(xlsx_file, 'TEMP_degC', rec)

disp('done')

%% -------- local functions --------
function rec = makeRec(sbj_ID, iteration, prefix, data, brain_mask, target_coord, target_cube_idx)
    data_brain = data;
    data_brain(~brain_mask) = 0;

    at_target   = data(target_coord(1), target_coord(2), target_coord(3));
    cube_vals   = data(target_cube_idx{1}, target_cube_idx{2}, target_cube_idx{3});
    around_max  = max(cube_vals(:));
    around_med  = median(cube_vals(:));

    rec = table( ...
        sbj_ID, iteration, string(prefix), ...
        max(data(:)), ...
        max(data_brain(:)), ...
        at_target, ...
        around_max, ...
        around_med, ...
        'VariableNames', {'sbj_ID','iteration','prefix', ...
                          'global_max','brain_max','at_target', ...
                          'around_max','around_med'});
end

function appendRow(fn, sheet, Tnew)
    % create file/sheet -> write header
    if ~isfile(fn)
        writetable(Tnew, fn, 'Sheet', sheet, 'WriteVariableNames', true)
        return
    end

    % check if sheet exists
    try
        sh = sheetnames(fn);
    catch
        [~, sh] = xlsfinfo(fn);
    end
    if ~any(strcmpi(sh, sheet))
        writetable(Tnew, fn, 'Sheet', sheet, 'WriteVariableNames', true)
        return
    end

    % sheet exists -> append without header
    try
        writetable(Tnew, fn, 'Sheet', sheet, ...
            'WriteMode', 'append', 'WriteVariableNames', false)
    catch   % older MATLAB without WriteMode
        Told = readtable(fn, 'Sheet', sheet);
        startRow = height(Told) + 2;          % +1 for header, +1 to start next row
        writetable(Tnew, fn, 'Sheet', sheet, ...
            'Range', sprintf('A%d', startRow), 'WriteVariableNames', false)
    end
end
