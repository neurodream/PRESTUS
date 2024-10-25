function [p] = NA_plot_mri_patch(mri, resize_factor, color)

% smooth mri
mri_smoothed = smooth3(mri, 'box', 7);

% Compute the isosurface
[f, v] = isosurface(mri_smoothed, 1);
v = v * resize_factor;

% Display the isosurface
p = patch('Faces', f, 'Vertices', v);
p.FaceColor = color;
p.EdgeColor = 'none';
alpha(p, 0.45);

p.SpecularColorReflectance = 0;
p.SpecularExponent = 50;

% Improve the visualization
daspect([1,1,1]);  % Adjust the data aspect ratio
view(3);           % Set the view to 3D
camlight;          % Add a camera light
lighting gouraud;  % Use Gouraud lighting for smooth shading
axis tight;        % Fit the axis tight around the data
grid on;           % Turn the grid on
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('MRI Volume Isosurface');

% % Optionally, add isonormals for smoother appearance
% isonormals(mri_smoothed, p);

end