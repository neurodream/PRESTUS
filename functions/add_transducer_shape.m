function add_transducer_shape(pos, target, varargin)
    % Input: pos - a vector [x, y, z] for translation
    %        target_vec - a target vector for the final direction
    
    default_params = [];
    default_params.transducer.curv_radius_mm = 62.94; % CTX500
    default_params.transducer.dist_to_plane_mm = 52.38; % CTX500
    default_params.grid_step_mm = 0.5;

    % Create an input parser
    p = inputParser;
    
    % Add parameters with default values
    % addRequired(p, 'data_folder', @ischar);
    addParameter(p, 'RealWorldCoords', false, @islogical);
    addParameter(p, 'Casing', false, @islogical);
    addParameter(p, 'TargetType', 'sphere', @ischar); % can be 'sphere', 'line', 'none', 'both'
    addParameter(p, 'PatchName', 'transducerPatch', @ischar);
    addParameter(p, 'Parameters', default_params, @isstruct);
    addParameter(p, 'Color', [0.5, 0.5, 0.5]);
    
    % Parse input arguments
    parse(p, varargin{:});
    
    % Retrieve values
    % data_folder     = p.Results.data_folder;
    real_world_coords = p.Results.RealWorldCoords;
    show_casing       = p.Results.Casing;
    target_type       = p.Results.TargetType;
    patch_name        = p.Results.PatchName;
    parameters        = p.Results.Parameters;
    color             = p.Results.Color;

    % Ensure 'pos' is a vector of length 3
    assert(numel(pos) == 3, 'Input "pos" must be a vector of length 3.');
    assert(numel(target) == 3, 'Input "target_vec" must be a vector of length 3.');

    target_vec = pos - target;

    % Normalize the target vector
    target_vec = target_vec / norm(target_vec);

    % Initial vector (pointing upwards)
    initial_vec = [0, 0, -1];

    % Compute the rotation axis (cross product)
    rotation_axis = cross(initial_vec, target_vec);
    if norm(rotation_axis) == 0  % If the vectors are collinear, no rotation is needed
        rotation_axis = [1, 0, 0];  % Default axis
    end
    rotation_axis = rotation_axis / norm(rotation_axis);

    % Compute the rotation angle (dot product)
    angle = acos(dot(initial_vec, target_vec));

    % Compute the rotation matrix using axis-angle representation
    R = axis_angle_to_rotation_matrix(rotation_axis, angle);

    %% Casing part 1
    % part 2, i.e. the rotation, comes after the bowl part, because the case-bowl
    % connection has to be created with unrotated patches

    % Create casing surface

    % TODO remove hardcode (not too important since CTX and Imasonics
    % transducer casings are relatively similar in diameter)
    [X, Y, Z] = cylinder(81.788/2, 50);
    Z = Z*19.05 - 19.05;
    
    % multiply by grid step size to arrive at voxel coordinates
    % (TODO: check if correct)
    if ~real_world_coords
        X = X*parameters.grid_step_mm;
        Y = Y*parameters.grid_step_mm;
        Z = Z*parameters.grid_step_mm;
    end

    hold on;

    %% Bowl part
    % Parameters
    radius = parameters.transducer.curv_radius_mm;
    dist_to_plane = parameters.transducer.dist_to_plane_mm;
    cut_height = radius - dist_to_plane;

    theta = linspace(0, pi, 50);  % Polar angle
    phi = linspace(0, 2*pi, 50);  % Azimuthal angle

    [theta, phi] = meshgrid(theta, phi);

    % Spherical coordinates
    x = radius * sin(theta) .* cos(phi);
    y = radius * sin(theta) .* sin(phi);
    z = radius * cos(theta) + radius;

    % Apply cut-off
    % mind the resolution of sphere coordinates - move the last z layer to
    % match the exit plane
    z(z > cut_height + radius/50) = NaN;

    % Move so that exit plane is at z = 0
    z = z - cut_height;

    % move the last z row to match the exit plane height
    % z(:, all(isnan(z))) = [];
    % z(:, 1) = zeros(size(z, 1), 1);
    firstNonNanCol = find(any(~isnan(z), 1), 1, 'first');
    z(:, firstNonNanCol) = zeros(size(z, 1), 1);

    % multiply by grid step size to arrive at voxel coordinates
    % (TODO: check if correct)
    if ~real_world_coords
        x = x*parameters.grid_step_mm;
        y = y*parameters.grid_step_mm;
        z = z*parameters.grid_step_mm;
    end

    % bowl-casing connection
    
    x_outer = X(1,1:end-1);
    y_outer = Y(1,1:end-1);
    z_outer = zeros(1, 50);
    x_inner = x(:,firstNonNanCol)';
    y_inner = y(:,firstNonNanCol)';
    z_inner = zeros(1, 50);
    
    % x_bowl2case = [x_outer, fliplr(x_inner)];
    % y_bowl2case = [y_outer, fliplr(y_inner)];
    % z_bowl2case = zeros(1, 2*50);  % Since it's a flat annulus, z = 0 for all points
    % 
    % [x_bowl2case, y_bowl2case, z_bowl2case] = apply_rotation(x_bowl2case, y_bowl2case, z_bowl2case, R);
    % 
    % bowl2case = patch(x_bowl2case + pos(1), y_bowl2case + pos(2), z_bowl2case + pos(3), 'b', 'FaceAlpha', 0.5);

    % Combine the inner and outer circle into one set of vertices for the patch
    x_combined = [x_outer, x_inner];
    y_combined = [y_outer, y_inner];
    z_combined = [z_outer, z_inner];

    [x_combined, y_combined, z_combined] = apply_rotation(x_combined, y_combined, z_combined, R);
    
    % Create vertices matrix (each row is a vertex)
    vertices = [x_combined(:) + pos(1), y_combined(:) + pos(2), z_combined(:) + pos(3)];
    
    % Create faces (connect corresponding points between inner and outer circle)
    faces = [];
    for i = 1:50-1
        % Quad face for each segment of the annulus
        faces = [faces; i, i+1, i+50+1, i+50];
    end
    
    % Close the loop by connecting the last points to the first
    faces = [faces; 50, 1, 50+1, 2*50];
    
    % Create the patch
    bowl2case = patch('Vertices', vertices, 'Faces', faces, ...
                      'FaceColor', color); % , 'FaceAlpha', 0.5


    % remove NaNs from bowl
    x(:, all(isnan(z))) = [];
    y(:, all(isnan(z))) = [];
    z(:, all(isnan(z))) = [];
    
    % Rotate bowl
    [x, y, z] = apply_rotation(x, y, z, R);

    % Create the spherical patch with cut-off
    bowl = surf( ...
        x + pos(1), ...
        y + pos(2), ...
        z + pos(3) ...
        );
    % bowl = surf(x, y, z); % TODO remove debug

    %% casing part 2

    if show_casing

        % Rotate casing
        [X, Y, Z] = apply_rotation(X, Y, Z, R);

        casing = surf(X + pos(1), Y + pos(2), Z + pos(3));
    end

    %% target part

    if strcmp(target_type, 'sphere')
        target_patch = plot_sphere(5, target);
    elseif strcmp(target_type, 'line')
        line = plot3([target(1) pos(1)], [target(2) pos(2)], [target(3) pos(3)], color, 'LineWidth', 2);
    elseif strcmp(target_type, 'both')
        target_patch = plot_sphere(5, target);
        line = plot3([target(1) pos(1)], [target(2) pos(2)], [target(3) pos(3)], color, 'LineWidth', 2);
    end

    %% Set patch properties
    casing.FaceColor = color;
    % casing.FaceAlpha = 0.6;%0.9;
    casing.EdgeColor = 'none';
    casing.Tag = patch_name;
    
    bowl.FaceColor = color;
    % bowl.FaceAlpha = 0.6;%0.9;
    bowl.EdgeColor = 'none';
    bowl.Tag = patch_name;

    bowl2case.FaceColor = color;
    % bowl2case.FaceAlpha = 0.6;%0.9;
    bowl2case.EdgeColor = 'none';
    bowl2case.Tag = patch_name;

    target_patch.Tag = patch_name;

    line.Tag = patch_name;

    axis equal;

    % hold off;
end

function R = axis_angle_to_rotation_matrix(rotation_axis, angle)
    % Axis-angle to rotation matrix
    c = cos(angle);
    s = sin(angle);
    t = 1 - c;
    x = rotation_axis(1);
    y = rotation_axis(2);
    z = rotation_axis(3);

    R = [t*x*x + c,    t*x*y - s*z,  t*x*z + s*y;
         t*x*y + s*z,  t*y*y + c,    t*y*z - s*x;
         t*x*z - s*y,  t*y*z + s*x,  t*z*z + c];
end

function [X, Y, Z] = apply_rotation(X, Y, Z, R)
    % Apply rotation matrix R to the coordinates X, Y, Z
    pts = [X(:), Y(:), Z(:)] * R';
    X = reshape(pts(:, 1), size(X));
    Y = reshape(pts(:, 2), size(Y));
    Z = reshape(pts(:, 3), size(Z));
end