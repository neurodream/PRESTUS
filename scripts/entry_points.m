% script by Kenneth

clc; clear; close all;

%% Remove simnibs from the path to resolve repelem.m conflicts
simnibs_path = '/home/affneu/kenvdzee/.conda/envs/simnibs_env';
matlab_paths = strsplit(path, pathsep);
simnibs_on_path = any(strcmp(matlab_paths, simnibs_path));
if simnibs_on_path
	rmpath(genpath(simnibs_path))
end

%% Set subject id
subject_id = 28;

%% Load csv with all PRESTUS coordinates
all_simulation_coordinates = readtable('/project/3025011.02/TUS_simulations/planning/planning_coordinate_list.csv');
% Filter out everything but the subject's rows
subject_simulation_coordinates = all_simulation_coordinates(all_simulation_coordinates.subject_id == subject_id, :);

%% Retreive the affine matrix from the anatomical file
subject_nifti_file = sprintf('sub-%03d*mprage_T1w.nii.gz', subject_id);
subject_nifti_path = sprintf('/project/3025011.02/bids/sub-%03d/ses-mri01/anat/', subject_id);
subject_nifti_file_and_path = fullfile(subject_nifti_path, subject_nifti_file);
subject_nifti_file_and_path = fullfile(subject_nifti_path, dir(subject_nifti_file_and_path).name);
subject_nifti_info = niftiinfo(subject_nifti_file_and_path);
subject_affine_matrix = subject_nifti_info.Transform.T';
% Validate it
if ~isequal(size(subject_affine_matrix), [4, 4])
    error('Affine matrix must be a 4x4 matrix.');
end

%% Extract voxel coordinates
% Separate coordinates for the entry points from the focus points
voxel_coords_entry = subject_simulation_coordinates{:, 3:5};
voxel_coords_focus = subject_simulation_coordinates{:, 6:8};

%% Apply affine transformation
% Add extra column so that it can be multiplied by the 4x4 affine matrix
number_of_targets = height(subject_simulation_coordinates);
voxel_coords_entry_affine_size = [voxel_coords_entry, ones(number_of_targets, 1)];
voxel_coords_focus_affine_size = [voxel_coords_focus, ones(number_of_targets, 1)];
% Multiply affine matrix with the original matrix
ras_coords_entry = (subject_affine_matrix * voxel_coords_entry_affine_size')';
ras_coords_focus = (subject_affine_matrix * voxel_coords_focus_affine_size')';

%% Extract RAS coordinates
ras_coords_entry = ras_coords_entry(:, 1:3);
ras_coords_focus = ras_coords_focus(:, 1:3);

%% Wrap the new coordinates with the original table
subject_RAS_coordinates = subject_simulation_coordinates;
subject_RAS_coordinates{:, 3:5} = ras_coords_entry;
subject_RAS_coordinates{:, 6:8} = ras_coords_focus;
% Rename the columns
subject_RAS_coordinates.Properties.VariableNames(3:5) = {'entry_x', 'entry_y', 'entry_z'};
subject_RAS_coordinates.Properties.VariableNames(6:8) = {'focus_x', 'focus_y', 'focus_z'};
% Round to 3 decimals
subject_RAS_coordinates{:, 3:8} = round(subject_RAS_coordinates{:, 3:8}, 3);
% Display the updated table
disp(subject_RAS_coordinates);

%% Remove the entrypoint columns
subject_RAS_coordinates(:, 3:5) = [];

%% Write table to subject folder
mkdir(sprintf('/project/3025011.02/localite/sub-%03d', subject_id));
writetable(subject_RAS_coordinates,sprintf('/project/3025011.02/localite/sub-%03d/sub-%03d_localite_coordinates.csv', subject_id, subject_id),'Delimiter',';')  