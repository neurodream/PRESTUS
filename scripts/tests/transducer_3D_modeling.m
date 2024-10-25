%% goal

% just get a patch showing the curvature of the bowl

%% measurements

% measurements from https://neurofus.com/products/
% Transducer Height:   0.75'' == 19.05mm
% Transducer Diameter: 3.22'' == 81.788mm

% measurements from config file:
% curv_radius_mm: 62.94mm
% dist_to_plane_mm: 52.38mm
% max(Elements_OD_mm): 60.83mm

%% casing part

[X,Y,Z] = cylinder(81.788/2);
Z = Z*19.05 - 19.05;

figure;
casing = surf(X/2,Y/2,Z/2);

hold on;

%% bowl part

% Parameters
radius = 62.94;
dist_to_plane = 52.38;
cut_height = radius - dist_to_plane;
theta = linspace(0, pi, 50);  % Polar angle
phi = linspace(0, 2*pi, 50);  % Azimuthal angle

[theta, phi] = meshgrid(theta, phi);

% Spherical coordinates
x = radius * sin(theta) .* cos(phi);
y = radius * sin(theta) .* sin(phi);
z = radius * cos(theta) + radius;

% Apply cut-off
z(z > cut_height) = NaN;

% move so that exit plane is at z = 0
z = z - cut_height;

% Create the spherical patch with cut-off
bowl = surf(x/2, y/2, z/2);

axis equal;
shading interp;

casing.FaceColor = [0.5, 0.5, 0.5];
casing.FaceAlpha = 0.9;
bowl.FaceColor = [0.5, 0.5, 0.5];
bowl.FaceAlpha = 0.9;