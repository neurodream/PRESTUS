function [ellipsoid_mask] = get_FWHM(data, threshold, plot, verbose)

% NOTE! assumes axis-aligned ellipsoid.

% Normalize the data
data = data / max(data(:));

% Find voxels above the threshold
above_threshold = data >= threshold;

% axis aligned:

% Get the indices of the voxels above the threshold
[x, y, z] = ind2sub(size(data), find(above_threshold));

% Compute the bounds of the FWHM ellipsoid along each axis
min_x = min(x); max_x = max(x);
min_y = min(y); max_y = max(y);
min_z = min(z); max_z = max(z);

% Calculate the semi-axes lengths of the ellipsoid
a = (max_x - min_x) / 2; % Semi-major axis
b = (max_y - min_y) / 2; % First minor axis
c = (max_z - min_z) / 2; % Second minor axis

% Center of the ellipsoid
center = [(min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2];

% Create a logical 3D matrix for the ellipsoid mask
[X, Y, Z] = ndgrid(1:size(data, 1), 1:size(data, 2), 1:size(data, 3));
ellipsoid_mask = ((X - center(1)) / a).^2 + ...
                 ((Y - center(2)) / b).^2 + ...
                 ((Z - center(3)) / c).^2 <= 1;

if plot
    % Extract a 2D slice (e.g., middle slice along z-axis)
    slice_index = round(center(3)); % Centered in the z-dimension
    slice_2D = squeeze(data(:, :, slice_index));
    
    % Generate 2D projection of the ellipsoid on the slice
    theta = linspace(0, 2*pi, 100);
    ellipsoid_y = a * cos(theta) + center(1); % Align to X in slice
    ellipsoid_x = b * sin(theta) + center(2); % Align to Y in slice
    
    % Plot the 2D slice and overlay the ellipsoid
    figure;
    imagesc(slice_2D);
    hold on;
    plot(ellipsoid_x, ellipsoid_y, 'r', 'LineWidth', 2); % Note the swapped axes
    xlabel('X Axis');
    ylabel('Y Axis');
    title('2D Snapshot with FWHM Ellipsoid Overlay');
    axis equal;
    colorbar;
    hold off;
end

if verbose
    % Display the ellipsoid parameters
    fprintf('Ellipsoid semi-axes: a = %.2f, b = %.2f, c = %.2f\n', a, b, c);
    fprintf('Ellipsoid center: x = %.2f, y = %.2f, z = %.2f\n', center(1), center(2), center(3));
end

end