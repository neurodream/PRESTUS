function plot_transducer_pos(parameters, sbj_ID, save, varargin)

% close all;
% clc;

p = inputParser;
    
% Required positional arguments
addRequired(p, 'parameters', @isstruct);
addRequired(p, 'sbj_ID', @isnumeric);

% Optional positional argument
addOptional(p, 'save', false, @islogical);

% Name-value pair arguments
addParameter(p, 'Structural', 'scalp', @ischar); % 'scalp', 'skull' or 'none'
addParameter(p, 'Functional', 'pressure', @ischar); % 'pressure', 'intensity', 'mechanicalindex', 'maxtemp', 'thermaldose' or 'none'
addParameter(p, 'CutoffPerc', 0.999, @isnumeric); % (TODO maybe add median option)
addParameter(p, 'LowCutoff', 0, @isnumeric); % absolute number; overwrites CutoffPerc if > 0

% Parse input arguments
parse(p, parameters, sbj_ID, varargin{:});



% Read the data from the Excel file
T = readtable('data/transducer_pos/position_LUT.xlsx');
target_L = [T.x_l(T.sbj_ID == sbj_ID) T.y_l(T.sbj_ID == sbj_ID) T.z_l(T.sbj_ID == sbj_ID)];
target_R = [T.x_r(T.sbj_ID == sbj_ID) T.y_r(T.sbj_ID == sbj_ID) T.z_r(T.sbj_ID == sbj_ID)];
% TODO not sure why swapping x and y dimensions necessary here
target_L = target_L([2 1 3]);
target_R = target_R([2 1 3]);

% make sure the subject ID match of seg_file and target:
segmentation_folder = fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', sbj_ID));
filename_segmented = fullfile(segmentation_folder, 'final_tissues.nii.gz');

layers = niftiread(filename_segmented);
layers_info = niftiinfo(filename_segmented);
[layers, layers_info] = swapNiftiXY(layers, layers_info); % x and y seem swapped in nifti, need to be swapped back

head = layers > 0;

head = fill_head(head);

skull = layers == 7 | layers == 8;

transformMatrix = layers_info.Transform.T;
parameters.transform = transformMatrix;
parameters.grid_step_mm = mean([transformMatrix(1,1) transformMatrix(2,2) transformMatrix(3,3)]); % TODO check in Julian's code if valid
% parameters.grid_step_mm = 0.5;

%% create figure with head/skull and targets

figure;
if strcmp(p.Results.Structural, 'scalp')
    head_smooth = smooth3(head, 'box', 5);
    struct_patch = patch(isosurface(head_smooth, 0.5)); % Extract and plot outer layer
    set(struct_patch, 'FaceAlpha', 0.25, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none'); % Customize appearance
    set(struct_patch, 'AmbientStrength', 0.3, 'DiffuseStrength', 0.5, 'SpecularStrength', 0.2, 'SpecularExponent', 1);
    isonormals(head_smooth, struct_patch); % Add normals for proper lighting
end

hold on;

if strcmp(p.Results.Structural, 'skull')
    skull_smooth = smooth3(skull, 'box', 5);
    struct_patch = patch(isosurface(skull_smooth, 0.5)); % Extract and plot outer layer
    set(struct_patch, 'FaceAlpha', 0.25, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none'); % Customize appearance
    set(struct_patch, 'AmbientStrength', 0.3, 'DiffuseStrength', 0.5, 'SpecularStrength', 0.2, 'SpecularExponent', 1);
    isonormals(skull_smooth, struct_patch); % Add normals for proper lighting
end


% entry point left transducer

V = struct_patch.Vertices;    % Nx3
F = struct_patch.Faces;       % Mx3

P1 = parameters.transducers(1).pos_t1_grid;
P2 = parameters.transducers(1).focus_pos_t1_grid;
d  = P2 - P1;                 % direction

X = [];
for k = 1:size(F,1)
    tri = V(F(k,:),:);
    n   = cross(tri(2,:)-tri(1,:),tri(3,:)-tri(1,:));
    denom = n*d';
    if abs(denom)<eps, continue, end
    t = n*(tri(1,:)-P1)'/denom;
    if t<0 || t>1, continue, end
    P = P1 + t*d;
    % inside-triangle test (barycentric)
    u = tri(2,:)-tri(1,:); v = tri(3,:)-tri(1,:);
    w = P - tri(1,:);
    denom2 = dot(u,u)*dot(v,v)-dot(u,v)^2;
    s = (dot(u,u)*dot(w,v)-dot(u,v)*dot(w,u))/denom2;
    r = (dot(v,v)*dot(w,u)-dot(u,v)*dot(w,v))/denom2;
    if s>=0 && r>=0 && s+r<=1
        X = P; break
    end
end

disp('ENTRY LEFT TRANSDUCER')
disp(X);



P1 = parameters.transducers(2).pos_t1_grid;
P2 = parameters.transducers(2).focus_pos_t1_grid;
d  = P2 - P1;                 % direction

X = [];
for k = 1:size(F,1)
    tri = V(F(k,:),:);
    n   = cross(tri(2,:)-tri(1,:),tri(3,:)-tri(1,:));
    denom = n*d';
    if abs(denom)<eps, continue, end
    t = n*(tri(1,:)-P1)'/denom;
    if t<0 || t>1, continue, end
    P = P1 + t*d;
    % inside-triangle test (barycentric)
    u = tri(2,:)-tri(1,:); v = tri(3,:)-tri(1,:);
    w = P - tri(1,:);
    denom2 = dot(u,u)*dot(v,v)-dot(u,v)^2;
    s = (dot(u,u)*dot(w,v)-dot(u,v)*dot(w,u))/denom2;
    r = (dot(v,v)*dot(w,u)-dot(u,v)*dot(w,v))/denom2;
    if s>=0 && r>=0 && s+r<=1
        X = P; break
    end
end

disp('ENTRY RIGHT TRANSDUCER')
disp(X);



xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
axis equal;
view(3);
% camlight('left');      % Light from the left side
% camlight('right');     % Light from the right side
light('Position', [-1, 0, 0], 'Style', 'infinite', 'Color', [1 1 1]*0.5);
light('Position', [1, 0, 0],  'Style', 'infinite', 'Color', [1 1 1]*0.5); 
light('Position', [0, -1, 0], 'Style', 'infinite', 'Color', [1 1 1]*0.5);
light('Position', [0, 1, 0],  'Style', 'infinite', 'Color', [1 1 1]*0.5); 
light('Position', [0, 0, -1], 'Style', 'infinite', 'Color', [1 1 1]*0.5);
light('Position', [0, 0, 1],  'Style', 'infinite', 'Color', [1 1 1]*0.5); 
% camlight('headlight'); % does not seem to work...
lighting gouraud;
material([0.4 0.6 0.2 20]); % Softer shininess with adjusted properties
% material dull;

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
    get_transducer_voxels_absolute_pos(target, pos, head, parameters, transducer.name, color);

end

% if strcmp(p.Results.Functional, 'intensity') % TODO add others
    add_sim_result_patch(parameters, sbj_ID, p.Results.Functional, 'HeadData', layers, 'CutoffPerc', p.Results.CutoffPerc, 'LowCutoff', p.Results.LowCutoff); % 0.9998
% end

% plot3(target_L(1), target_L(2), target_L(3), 'k.', 'MarkerSize', 20);
% plot3(target_R(1), target_R(2), target_R(3), 'k.', 'MarkerSize', 20);

% instead of just a point: plot the quantified ROI

r = round(5/mean(layers_info.PixelDimensions));
% Equation of the sphere: (x - px)^2 + (y - py)^2 + (z - pz)^2 <= r^2
[x, y, z] = ndgrid(1:size(layers,1), 1:size(layers,2), 1:size(layers,3));
ROItarget_L = (x - target_L(1)).^2 + (y - target_L(2)).^2 + (z - target_L(3)).^2 <= r^2;
ROItarget_R = (x - target_R(1)).^2 + (y - target_R(2)).^2 + (z - target_R(3)).^2 <= r^2;

ROItarget_L_smooth = smooth3(ROItarget_L, 'box', 5);
target_patch = patch(isosurface(ROItarget_L_smooth, 0.5)); % Extract and plot outer layer
set(target_patch, 'FaceAlpha', 1, 'FaceColor', 'green', 'EdgeColor', 'none'); % Customize appearance
isonormals(ROItarget_L_smooth, target_patch); % Add normals for proper lighting

ROItarget_R_smooth = smooth3(ROItarget_R, 'box', 5);
target_patch = patch(isosurface(ROItarget_R_smooth, 0.5)); % Extract and plot outer layer
set(target_patch, 'FaceAlpha', 1, 'FaceColor', 'green', 'EdgeColor', 'none'); % Customize appearance
isonormals(ROItarget_R_smooth, target_patch); % Add normals for proper lighting

view(-200,43);

if save
    data_folder = fullfile(parameters.data_path, 'sim_outputs', sprintf('sub-%03d', sbj_ID));
    figure_file = fullfile(data_folder, sprintf('sub-%03d_3Dplot%s.fig', sbj_ID, parameters.results_filename_affix));
    savefig(figure_file);

    image_file = fullfile(data_folder, sprintf('sub-%03d_3Dplot%s.png', sbj_ID, parameters.results_filename_affix));
    % Remove grid
    axis off; % grid off;
    % Set axes and figure background to none
    set(gca, 'Color', 'none'); % Make axes background transparent
    set(gcf, 'Color', 'none'); % Make figure background transparent
    % Save the figure as a PNG image
    exportgraphics(gcf, image_file, 'BackgroundColor', 'none');

end