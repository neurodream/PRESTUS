function [vertices, faces] = get_directions(grid_dims, target)
    
    % find furthest dim edge from target
    r = 0;
    for dim_i = numel(grid_dims)
        dim = grid_dims(dim_i);
        if target(dim_i) - 1 > r
            r = target(dim_i) - 1;
        end
        if dim - target(dim_i) > r
            r = dim - target(dim_i);
        end
    end
    
    % Approximate number of divisions needed in the polar and azimuthal angles
    n_theta = ceil(pi * r);       % Polar angle divisions (0 to pi)
    n_phi = ceil(2 * pi * r); % Azimuthal angle divisions (0 to 2*pi)
    
    % Define spherical coordinates (theta for polar angle, phi for azimuthal angle)
    theta = linspace(0, pi, n_theta);        % Polar angle (0 to pi)
    phi = linspace(0, 2*pi, n_phi);          % Azimuthal angle (0 to 2*pi)
    
    % Create a grid for spherical coordinates (2D grid of theta and phi for fixed radius)
    [Theta, Phi] = meshgrid(theta, phi);

    % Remove the poles (first and last rows of Theta and Phi)
    Theta_no_poles = Theta(2:end-1, :);  % Exclude first and last rows
    Phi_no_poles = Phi(2:end-1, :);
    
    % % Cap the X, Y, and Z coordinates to the valid range in each dimension
    % X = max(1, min(grid_dims(1), X));  % Cap X between 1 and dim_x
    % Y = max(1, min(grid_dims(2), Y));  % Cap Y between 1 and dim_y
    % Z = max(1, min(grid_dims(3), Z));  % Cap Z between 1 and dim_z
    
    % Convert spherical coordinates (Theta_no_poles, Phi_no_poles) to Cartesian coordinates (X, Y, Z)
    X_no_poles = r * sin(Theta_no_poles) .* cos(Phi_no_poles);  % X-coordinates
    Y_no_poles = r * sin(Theta_no_poles) .* sin(Phi_no_poles);  % Y-coordinates
    Z_no_poles = r * cos(Theta_no_poles);                       % Z-coordinates

    % Flatten the matrices to create a list of vertices (after removing poles)
    vertices = [X_no_poles(:), Y_no_poles(:), Z_no_poles(:)];
    
    % Create faces by connecting adjacent vertices (this will create quadrilateral faces)
    faces = [];
    for i = 1:n_phi-1
        for j = 1:n_theta-2  % Adjusted to account for the removed poles
            % Define the four corners of the quadrilateral face
            p1 = (i-1) * (n_theta-2) + j;         % Vertex at (i,j)
            p2 = (i-1) * (n_theta-2) + (j + 1);   % Vertex at (i,j+1)
            p3 = i * (n_theta-2) + (j + 1);       % Vertex at (i+1,j+1)
            p4 = i * (n_theta-2) + j;             % Vertex at (i+1,j)
            
            % Append the face (four vertices)
            faces = [faces; p1, p2, p3, p4];
        end
    end

end