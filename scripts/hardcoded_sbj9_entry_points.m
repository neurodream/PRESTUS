% sub9
entry_left = [25.8200  163.0000  147.0000];
entry_right = [196.5000  163.0000  149.0000];
% sub7
entry_left = [31.4797 182.0000 145.7532];
entry_right = [196.3276 182.0000 144.8451];

subject_nifti_info = niftiinfo('/home/sleep/nicade/Documents/scans/sub-009_T1.nii.gz');
subject_nifti_info = niftiinfo('/home/sleep/nicade/Documents/scans/sub-007_T1.nii.gz');
subject_affine_matrix = subject_nifti_info.Transform.T';

% Validate it
if ~isequal(size(subject_affine_matrix), [4, 4])
    error('Affine matrix must be a 4x4 matrix.');
end

voxel_coords_entry_L = entry_left;
voxel_coords_entry_R = entry_right;
% voxel_coords_focus = [131 163 149]; % sub9
voxel_coords_focus = [130 182 137]; %sub7


number_of_targets = 1;
voxel_coords_entry_affine_size_L = [voxel_coords_entry_L, ones(number_of_targets, 1)];
voxel_coords_entry_affine_size_R = [voxel_coords_entry_R, ones(number_of_targets, 1)];
voxel_coords_focus_affine_size = [voxel_coords_focus, ones(number_of_targets, 1)];
% Multiply affine matrix with the original matrix
ras_coords_entry_L = (subject_affine_matrix * voxel_coords_entry_affine_size_L')';
ras_coords_entry_R = (subject_affine_matrix * voxel_coords_entry_affine_size_R')';
ras_coords_focus = (subject_affine_matrix * voxel_coords_focus_affine_size')';

%% Extract RAS coordinates
ras_coords_entry_L = ras_coords_entry_L(:, 1:3);
ras_coords_entry_R = ras_coords_entry_R(:, 1:3);
ras_coords_focus = ras_coords_focus(:, 1:3);

disp('left entry RAS:');
disp(round(ras_coords_entry_L, 3));
disp('right entry RAS:');
disp(round(ras_coords_entry_R, 3));
disp('focus pos RAS:');
disp(round(ras_coords_focus, 3));


%% backtransform RAS to MNI

%% RAS -> MNI
% set one of the two blocks below

%% inputs
% inputs
t1 = '/home/sleep/nicade/Documents/scans/sub-009_T1.nii.gz';
left_entry_ras = [-83.362 31.322 -1.336];                 % 1x3 RAS mm, or Nx3
right_entry_ras = [83.153 22.124 13.122];
target = [131 163 149];
ras = right_entry_ras;

% header affine (voxel -> RAS mm)
A = niftiinfo(t1).Transform.T';

% RAS mm -> voxel (floating, i j k)
ijk = (A \ [ras ones(size(ras,1),1)]')';
ijk = ijk(:,1:3);

% integer voxel indices for array addressing (MATLAB is 1-based)
ijk_idx = round(ijk);                            % use round/ceil/floor as needed

