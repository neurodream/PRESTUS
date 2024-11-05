function [parameters, distance] = get_transducer_pos(parameters, sbj_ID, side, transd_ind, transd_angle, transd_pos_add, target_pos_add, to_other)

% overwrites the respective parameters fields

% TODO (identify cause of and) get rid of x y flip

% temporarily add a transducer field so old PRESTUS functions can be used
parameters.transducer = parameters.transducers(transd_ind);

% Read the data from the Excel file
data = readtable('data/transducer_pos/position_LUT.xlsx');

% Find the row with the specified sbj_ID
row_idx = find(data.sbj_ID == sbj_ID);

% Extract the respective coordinates
% TODO careful: x and y are swapped here
target_L = [data.y_l(row_idx), data.x_l(row_idx), data.z_l(row_idx)];
target_L = target_L + target_pos_add([2 1 3]);
target_R = [data.y_r(row_idx), data.x_r(row_idx), data.z_r(row_idx)];
target_R = target_R + target_pos_add([2 1 3]);
if strcmp(side, 'l')
    target = target_L;
    other = target_R;
elseif strcmp(side, 'r')
    target = target_R;
    other = target_L;
else
    error('Invalid side specified. Use ''l'' for left or ''r'' for right.');
end

%%

% make sure the subject ID match of seg_file and target:
seg_file = fullfile(parameters.seg_path, ['m2m_sub-' sprintf('%03d', sbj_ID)], 'final_tissues.nii.gz');

layers = niftiread(seg_file);
layers_info = niftiinfo(seg_file);
head = layers > 0;

head = fill_head(head);

skull = layers == 7 | layers == 8;

transformMatrix = layers_info.Transform.T;
parameters.transform = transformMatrix;
parameters.grid_step_mm = mean([transformMatrix(1,1) transformMatrix(2,2) transformMatrix(3,3)]); % TODO check in Julian's code if valid

%% calculate effective targets

% location of the center of the exit plane (?) of the transducer
[~, source] = get_transducer_voxels(target, transd_angle, head, parameters, '', 'red');

target_contralateral = point_line_projection_3D(other, source, target);

%% add transducer

% remove a transducer by calling: delete(findobj('Tag', transducer_name));

if to_other
    target = target_contralateral;
end
[~, trans_pos, distance] = get_transducer_voxels(target, transd_angle, head, parameters, 'transd1', '#0072BD');

trans_pos = trans_pos([2 1 3]); % flip back to normal space
trans_pos = trans_pos + transd_pos_add;

parameters.transducers(transd_ind).pos_t1_grid = round(trans_pos);
parameters.transducers(transd_ind).focus_pos_t1_grid = round(target([2 1 3]));  % flip back to normal space

parameters = rmfield(parameters, 'transducer');

% TODO careful: assuming no other figures/important command window
% information is open
clc; close all;

end