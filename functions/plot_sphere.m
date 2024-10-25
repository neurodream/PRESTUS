function sphere_patch = plot_sphere(r, p)
    % r - radius of the sphere
    % p - center position [x, y, z]
    
    % Ensure 'p' is a vector of length 3
    assert(numel(p) == 3, 'Input "p" must be a vector of length 3.');
    
    % Generate sphere data
    [X, Y, Z] = sphere(50); % You can change 50 to control the resolution
    
    % Scale by radius and translate to the center position
    X = X * r + p(1);
    Y = Y * r + p(2);
    Z = Z * r + p(3);
    
    % Create the spherical patch
    sphere_patch = surf(X, Y, Z);
    
    % Set properties (optional)
    sphere_patch.FaceColor = [1, 0, 0];
    sphere_patch.FaceAlpha = 0.9;
    sphere_patch.EdgeColor = 'none';

    sphere_patch.Tag = 'target_sphere';
end