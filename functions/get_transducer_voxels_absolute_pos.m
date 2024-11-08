function [transducer_voxels] = get_transducer_voxels_absolute_pos(target, trans_pos, matter, parameters, patch_name, transducer_color)

    grid_dims = size(matter);
    transducer_voxels = update_transducer_voxels(trans_pos, target, parameters, grid_dims, transducer_color);
    overlap = matter & transducer_voxels;
    collision = any(overlap(:));

    if collision
        disp('WARNING!! transducer within skin');
    end

    % plot
    line = plot3( ...
        [target(1), trans_pos(1)], [target(2), trans_pos(2)], [target(3), trans_pos(3)], ...
        'Color', transducer_color, 'LineWidth', 4); %
    
    line.Tag = patch_name;
    add_transducer_shape(trans_pos, target, 'Parameters', parameters, 'TargetType', 'none', 'PatchName', patch_name, 'Color', transducer_color);
    
    % graphicsArray = findobj('Tag', patch_name);
    % for i = 1:numel(graphicsArray)
    %     set(graphicsArray(i), 'Color', transducer_color); % Set to red (adjust color as needed)
    % end

    % p_trans = patch(isosurface(transducer_voxels, 0.5)); % Extract and plot outer layer
    % set(p_trans, 'FaceAlpha', 0.5, 'FaceColor', transducer_color, 'EdgeColor', 'none'); % Customize appearance
    % isonormals(transducer_voxels, p_trans); % Add normals for proper lighting
    % p_trans.Tag = patch_name;

    %% compute length of line
    p1 = target;
    p2 = trans_pos;
    transformMatrix = parameters.transform;

    % Convert points to homogeneous coordinates (adding a 1 for the 4th coordinate)
    p1_homogeneous = [p1, 1];
    p2_homogeneous = [p2, 1];
    
    % Apply the affine transformation to both points
    p1_transformed = (transformMatrix * p1_homogeneous')';
    p2_transformed = (transformMatrix * p2_homogeneous')';
    
    % Extract the transformed coordinates (first 3 elements, ignore the 4th element)
    p1_transformed = p1_transformed(1:3);
    p2_transformed = p2_transformed(1:3);
    
    % Compute the Euclidean distance (length of the line) between the transformed points
    line_length = sqrt(sum((p2_transformed - p1_transformed).^2));
    
    % Display the length of the line
    disp(['transducer-target distance: ', num2str(line_length), ' mm']);

    % also show the transducer coords
    disp(['transducer coordinates: ', round(num2str(trans_pos(1))), ' ', round(num2str(trans_pos(2))), ' ', round(num2str(trans_pos(3)))]);

end