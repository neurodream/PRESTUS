function plot_transducer_pos(parameters, sbj_ID, plot_scalp, plot_skull, plot_intensity, save)

close all; clc;

% Read the data from the Excel file
T = readtable('data/transducer_pos/position_LUT.xlsx');
target_L = [T.y_l(T.sbj_ID == sbj_ID) T.x_l(T.sbj_ID == sbj_ID) T.z_l(T.sbj_ID == sbj_ID)];
target_R = [T.y_r(T.sbj_ID == sbj_ID) T.x_r(T.sbj_ID == sbj_ID) T.z_r(T.sbj_ID == sbj_ID)];

% make sure the subject ID match of seg_file and target:
segmentation_folder = fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', sbj_ID));
filename_segmented = fullfile(segmentation_folder, 'final_tissues.nii.gz');

layers = niftiread(filename_segmented);
layers_info = niftiinfo(filename_segmented);
head = layers > 0;

head = fill_head(head);

skull = layers == 7 | layers == 8;

transformMatrix = layers_info.Transform.T;
parameters.transform = transformMatrix;
parameters.grid_step_mm = mean([transformMatrix(1,1) transformMatrix(2,2) transformMatrix(3,3)]); % TODO check in Julian's code if valid

%% create figure with head/skull and targets

figure;
if plot_scalp
    head_smooth = smooth3(head, 'box', 5);
    p = patch(isosurface(head_smooth, 0.5)); % Extract and plot outer layer
    set(p, 'FaceAlpha', 0.5, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none'); % Customize appearance
    isonormals(head_smooth, p); % Add normals for proper lighting
end

hold on;

if plot_skull
    skull_smooth = smooth3(skull, 'box', 5);
    p_skull = patch(isosurface(skull_smooth, 0.5)); % Extract and plot outer layer
    set(p_skull, 'FaceAlpha', 0.5, 'FaceColor', 'black', 'EdgeColor', 'none'); % Customize appearance
    isonormals(skull_smooth, p_skull); % Add normals for proper lighting
end

plot3(target_L(1), target_L(2), target_L(3), 'k.', 'MarkerSize', 20); % red point at the target
plot3(target_R(1), target_R(2), target_R(3), 'k.', 'MarkerSize', 20); % green point at the left NBM

xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
axis equal;
view(3);
camlight('left');      % Light from the left side
camlight('right');     % Light from the right side
camlight('headlight'); % does not seem to work...
lighting gouraud;

for transducer = parameters.transducers
    
    % % plot beam path from left and right transducer, respectively
    % plot3([source_L(1), defacto_target_R(1)], [source_L(2), defacto_target_R(2)], [source_L(3), defacto_target_R(3)], 'k-', 'LineWidth', 1);
    % plot3([source_R(1), defacto_target_L(1)], [source_R(2), defacto_target_L(2)], [source_R(3), defacto_target_L(3)], 'k-', 'LineWidth', 1);
    
    %% add transducer
    
    % temporarily add transducer as field to parameters, to make old
    % function work
    parameters.transducer = transducer;

    % remove a transducer by calling: delete(findobj('Tag', transducer_name));
    
    target = transducer.focus_pos_t1_grid;
    pos = transducer.pos_t1_grid;
    
    if strcmp(transducer.name, 'L')
        color = '#0072BD';
    elseif strcmp(transducer.name, 'R')
        color = '#A2142F';
    else
        color = 'k';
    end
    get_transducer_voxels_absolute_pos(target([2 1 3]), pos([2 1 3]), head, parameters, transducer.name, color);

end

if plot_intensity
    add_sim_result_patch(parameters, sbj_ID, 'HeadData', layers);
end

view(62, 36);

if save
    data_folder = fullfile(parameters.data_path, 'sim_outputs', sprintf('sub-%03d', sbj_ID));
    figure_file = fullfile(data_folder, sprintf('sub-%03d_3Dplot%s.fig', sbj_ID, parameters.results_filename_affix));
    savefig(figure_file);
end