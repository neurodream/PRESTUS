close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(currentFile);
cd(rootpath); % repos/PRESTUS_forked/
cd ..

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC




base_dir = 'C:\Users\nicade\Documents\projects\BFTUS\sim_results_buffer'; % Set your base directory
sub_dirs = dir(fullfile(base_dir, 'sub-*')); % List subdirectories matching 'sub-xxx'

for i = 1:length(sub_dirs)
    if ~sub_dirs(i).isdir
        continue;
    end
    sub_name = sub_dirs(i).name;
    if ~regexp(sub_name, '^sub-\d{3}$') % Ensure it matches sub-xxx format
        continue;
    end
    
    sub_path = fullfile(base_dir, sub_name);

    subject_id = sscanf(sub_name, 'sub-%d'); % Extracts the integer part

    if subject_id > 6
        break
    end

    files = dir(fullfile(sub_path, '*mechanicalindex*.nii.gz')); % Get matching files

    for j = 1:length(files)

        filename = files(j).name;
        file_path = fullfile(sub_path, filename);

        data = niftiread(file_path);
        
        max_MI = max(data(:));

        if max_MI > 4

            disp([num2str(max_MI) '; ' filename]);

        end

    end
    


    % files = dir(fullfile(sub_path, '*parameters*.mat')); % Get matching files
    % 
    % % Filter out files containing "--"
    % valid_files = {files(~contains({files.name}, '--')).name};
    % 
    % % Loop through and load valid files
    % fprintf('Loading files in %s:\n', sub_name);
    % 
    % subject_id = sscanf(sub_name, 'sub-%d'); % Extracts the integer part
    % 
    % for j = 1:length(valid_files)
    %     file_path = fullfile(sub_path, valid_files{j});
    %     % fprintf('  Loading: %s\n', valid_files{j});
    %     load(file_path);
    % 
    %     % try
    %         % Process 'data' as needed
    %         filename_cropped_smoothed_skull_data = fullfile([base_dir filesep sub_name filesep 'debug'], ...
    %         sprintf('sub-%03d_%s_after_cropping_and_smoothing_%s%s.mat', ...
    %         subject_id, parameters.simulation_medium, string(parameters.transducers(1).name), parameters.results_filename_affix));
    % 
    %         load(filename_cropped_smoothed_skull_data); % gives medium_masks
    % 
    %     % catch
    %     % 
    %     %     disp('');
    %     % 
    %     %     continue;
    %     % end
    % 
    %     containers = {
    %         'pressure_max_brain', ...
    %         'pressure_max_scalp', ...
    %         'pressure_max_target', ...
    %         'pressure_max_offtarget', ...
    %         'pressure_avg_target', ...
    %         'pressure_avg_scalp', ...
    %         'pressure_avg_offtarget', ...
    %         'pressure_95_brain', ...
    %         'pressure_95_scalp', ...
    %         'mechanicalindex_max_whole', ...
    %         'mechanicalindex_max_brain', ...
    %         'mechanicalindex_max_scalp', ...
    %         'FWHM_in_ROI_perc', ...
    %     };
    % 
    % 
    %     ROI_r = round(parameters.focus_area_radius/parameters.grid_step_mm);
    % 
    %     % create output containers
    %     measure_data = []; % measures = {'pressure', 'intensity', 'mechanicalindex', 'heating', 'maxCEM43'};
    % 
    %     % % targets
    %     T = readtable('data/transducer_pos/position_LUT.xlsx');
    %     % targetL = [T.y_l(T.sbj_ID == sub_id) T.x_l(T.sbj_ID == sub_id) T.z_l(T.sbj_ID == sub_id)];
    %     targetR = [T.x_r(T.sbj_ID == sub_id) T.y_r(T.sbj_ID == sub_id) T.z_r(T.sbj_ID == sub_id)];
    % 
    %     % Equation of the sphere: (x - px)^2 + (y - py)^2 + (z - pz)^2 <= r^2
    %     [x, y, z] = ndgrid(1:size(medium_masks,1), 1:size(medium_masks,2), 1:size(medium_masks,3));
    %     % ROItarget_L = (x - targetL(1)).^2 + (y - targetL(2)).^2 + (z - targetL(3)).^2 <= r^2;
    %     ROItarget_R = (x - targetR(1)).^2 + (y - targetR(2)).^2 + (z - targetR(3)).^2 <= ROI_r^2;
    % 
    %     % segmenting tissues
    % 
    %     % Creates a logical skull mask and register skull_ids
    %     labels = fieldnames(parameters.layer_labels);
    %     skull_i = find(strcmp(labels, 'skull_cortical'));
    %     trabecular_i = find(strcmp(labels, 'skull_trabecular'));
    %     all_skull_ids = [skull_i, trabecular_i];
    %     value_masks.skull = ismember(medium_masks,all_skull_ids);
    %     brain_i = find(strcmp(labels, 'brain'));
    %     value_masks.brain = ismember(medium_masks,brain_i);
    %     skin_i = find(strcmp(labels, 'skin'));
    %     value_masks.scalp = ismember(medium_masks,skin_i);
    %     value_masks.whole = true(size(medium_masks));
    % 
    %     % calculations for inner brain
    %     SE = strel('sphere', 6);
    %     value_masks.skulldilated = imdilate(value_masks.skull, SE);
    %     % braineroded  = imerode(brain, SE); % not optimal: erodes around the gyri
    %     value_masks.braininner = value_masks.brain & ~value_masks.skulldilated;
    %     % manually defined targets
    %     value_masks.target = ROItarget_R;
    %     value_masks.offtarget = value_masks.brain & ~value_masks.target;
    % 
    % end





end