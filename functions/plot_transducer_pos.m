function [parameters] = plot_transducer_pos(parameters, sbj_ID)

% big TODO!!

% overwrites the respective parameters fields

% temporarily add a transducer field so old PRESTUS functions can be used
parameters.transducer = parameters.transducers(transd_ind);

% Read the data from the Excel file
data = readtable('data/transducer_pos/position_LUT.xlsx');

% Find the row with the specified sbj_ID
row_idx = find(data.sbj_ID == sbj_ID);

% Extract the respective coordinates
if strcmp(side, 'l')
    target = [data.x_l(row_idx), data.y_l(row_idx), data.z_l(row_idx)];
elseif strcmp(side, 'r')
    target = [data.x_r(row_idx), data.y_r(row_idx), data.z_r(row_idx)];
else
    error('Invalid side specified. Use ''l'' for left or ''r'' for right.');
end

defacto_target = target + target_pos_add;
defacto_target = defacto_target([2 1 3]); % flip to plotting space

%%

% make sure the subject ID match of seg_file and target:
seg_file = ['/home/sleep/nicade/Documents/scans/segmentation_results/m2m_sub-' sprintf('%03d', sbj_ID) '/final_tissues.nii.gz'];

layers = niftiread(seg_file);
layers_info = niftiinfo(seg_file);
head = layers > 0;

head = fill_head(head);

skull = layers == 7 | layers == 8;

transformMatrix = layers_info.Transform.T;
parameters.transform = transformMatrix;
parameters.grid_step_mm = mean([transformMatrix(1,1) transformMatrix(2,2) transformMatrix(3,3)]); % TODO check in Julian's code if valid

%% add transducer

% remove a transducer by calling: delete(findobj('Tag', transducer_name));

[~, trans_pos] = get_transducer_voxels(defacto_target, transd_angle, head, parameters, 'transd1', '#0072BD');

trans_pos = trans_pos([2 1 3]); % flip back to normal space
trans_pos = trans_pos + transd_pos_add;

parameters.transducers(transd_ind).pos_t1_grid = trans_pos;
parameters.transducers(transd_ind).focus_pos_t1_grid = defacto_target([2 1 3]);  % flip back to normal space

parameters = rmfield(parameters, 'transducer');

% TODO careful: assuming no other figures/important command window
% information is open
clc; close all;

end