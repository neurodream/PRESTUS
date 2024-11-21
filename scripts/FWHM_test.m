sbj_ID = 3;
measure = 'intensity'; % mechanicalindex pressure intensity
affix = 'it1_FW_imprecisionnone';
targeting = '--'; % -- -
data = niftiread(sprintf('/home/sleep/nicade/Documents/scans/sim_outputs/sub-%03d/sub-%03d_final_%sL+z%sr_R+z%sr_%s.nii.gz', sbj_ID, sbj_ID, measure, targeting, targeting, affix));

% Assume 'data' is the 3D matrix containing intensity values

% Normalize the data
data = data / max(data(:));

% FWHM threshold
threshold = 0.5;

% Find voxels above the threshold
above_threshold = data >= threshold;

%% axis aligned:

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

% Display the ellipsoid parameters
fprintf('Ellipsoid semi-axes: a = %.2f, b = %.2f, c = %.2f\n', a, b, c);
fprintf('Ellipsoid center: x = %.2f, y = %.2f, z = %.2f\n', center(1), center(2), center(3));



%% non-axis aligned (free rotation) (does not work yet):

% % Get the indices of the voxels above the threshold
% [x, y, z] = ind2sub(size(data), find(above_threshold));
% points = [x, y, z];
% 
% % Perform PCA on the points
% [coeff, ~, ~] = pca(points);
% 
% % Transform points into PCA space
% transformed_points = points * coeff;
% 
% % Get the extents of the ellipsoid in PCA space
% min_vals = min(transformed_points, [], 1);
% max_vals = max(transformed_points, [], 1);
% 
% % Semi-axes lengths
% a = (max_vals(1) - min_vals(1)) / 2; % Semi-major axis
% b = (max_vals(2) - min_vals(2)) / 2; % First minor axis
% c = (max_vals(3) - min_vals(3)) / 2; % Second minor axis
% 
% % Center of the ellipsoid in PCA space
% center = mean(transformed_points, 1);
% 
% % Generate a grid of points for the ellipsoid in 2D (major & first minor axes)
% theta = linspace(0, 2*pi, 100);
% ellipsoid_x = a * cos(theta);
% ellipsoid_y = b * sin(theta);
% 
% % Combine the ellipsoid points into 2D space
% ellipsoid_points = [ellipsoid_x', ellipsoid_y'];
% 
% % Rotate the ellipsoid back to PCA space
% ellipsoid_transformed = ellipsoid_points * coeff(:, 1:2)'; % Transform using PCA
% ellipsoid_transformed(:, 1) = ellipsoid_transformed(:, 1) + center(1); % Shift by center
% ellipsoid_transformed(:, 2) = ellipsoid_transformed(:, 2) + center(2);
% 
% % Swap x and y for correct alignment
% ellipsoid_transformed = ellipsoid_transformed(:, [2, 1]); % Swap for correct alignment
% 
% % Rotate the entire dataset to align with PCA axes
% [x_grid, y_grid, z_grid] = ndgrid(1:size(data, 1), 1:size(data, 2), 1:size(data, 3));
% grid_points = [x_grid(:), y_grid(:), z_grid(:)];
% 
% % Transform grid points into PCA space
% transformed_grid = grid_points * coeff;
% 
% % Interpolate the data onto the rotated grid
% rotated_data = reshape(interp3(data, ...
%     grid_points(:, 2), grid_points(:, 1), grid_points(:, 3), 'linear', 0), ...
%     size(data));
% 
% % Extract a 2D slice in the PCA space (aligned with major and first minor axes)
% slice_index = round(center(3)); % Middle slice in PCA space
% slice_2D = squeeze(rotated_data(:, :, slice_index));
% 
% % Plot the 2D slice and overlay the ellipsoid
% figure;
% imagesc(slice_2D);
% hold on;
% plot(ellipsoid_transformed(:, 1), ellipsoid_transformed(:, 2), 'r', 'LineWidth', 2); % Overlay
% xlabel('Major Axis (PCA)');
% ylabel('First Minor Axis (PCA)');
% title('2D Snapshot with Non-Axis-Aligned FWHM Ellipsoid');
% axis equal;
% colorbar;
% hold off;
% 
% % Display the ellipsoid parameters and PCA orientation
% fprintf('Ellipsoid semi-axes: a = %.2f, b = %.2f, c = %.2f\n', a, b, c);
% fprintf('Ellipsoid center: x = %.2f, y = %.2f, z = %.2f\n', center(1), center(2), center(3));
% disp('Principal axes (columns of coeff):');
% disp(coeff);


