clc; clear; close all;

subject = 2;

head      = niftiread(sprintf('M:/Documents/scans/segmentation_results/m2m_sub-%03d/final_tissues.nii.gz', subject));
head_info = niftiinfo(sprintf('M:/Documents/scans/segmentation_results/m2m_sub-%03d/final_tissues.nii.gz', subject));

within = ismember(head, 1:3); % 1:10

within = ~within;
within = imfill(within, 'holes');
within = ~within;

padsize = [256 256 256] - size(within);
within = padarray(within, padsize, 0, 'post');

figure;
hold on;

Ds_skull = smooth3(double(within));%smooth3(double(head==7));
skin_isosurface = isosurface(Ds_skull,0.5);
hiso = patch(skin_isosurface,...
   'FaceColor',[0.5 0.5 0.5],...
   'EdgeColor','none',...
   'facealpha',0.05);

NBM_info = niftiinfo(sprintf('D:/scans/sub-%03d_T1_seg_L.nii.gz', subject));
scaling_factors = head_info.PixelDimensions(1:3) ./ NBM_info.PixelDimensions;

NBM_L = niftiread(sprintf('D:/scans/sub-%03d_T1_seg_L.nii.gz', subject));
Ds_NBM_L = smooth3(NBM_L);
NBM_L_isosurface = isosurface(Ds_NBM_L,0.5);
% NBM_L_isosurface.vertices(:,2) = NBM_L_isosurface.vertices(:,2) + 32;
NBM_L_isosurface.vertices(:,1) = NBM_L_isosurface.vertices(:,1) * scaling_factors(1);
NBM_L_isosurface.vertices(:,2) = NBM_L_isosurface.vertices(:,2) * scaling_factors(2);
NBM_L_isosurface.vertices(:,3) = NBM_L_isosurface.vertices(:,3) * scaling_factors(3);
hiso_NBM_L = patch(NBM_L_isosurface,...
   'FaceColor',[1.0 0.0 0.0],...
   'EdgeColor','none');

NBM_R = niftiread(sprintf('D:/scans/sub-%03d_T1_seg_R.nii.gz', subject));
Ds_NBM_R = smooth3(NBM_R);
NBM_R_isosurface = isosurface(Ds_NBM_R,0.5);
% NBM_R_isosurface.vertices(:,2) = NBM_R_isosurface.vertices(:,2) + 32;
NBM_R_isosurface.vertices(:,1) = NBM_R_isosurface.vertices(:,1) * scaling_factors(1);
NBM_R_isosurface.vertices(:,2) = NBM_R_isosurface.vertices(:,2) * scaling_factors(2);
NBM_R_isosurface.vertices(:,3) = NBM_R_isosurface.vertices(:,3) * scaling_factors(3);
hiso_NBM_R = patch(NBM_R_isosurface,...
   'FaceColor',[1.0 0.0 0.0],...
   'EdgeColor','none');

% prevent the ears
sphere_radius = 30;
center_L = [25, 120, 100];
center_R = [200, 120, 100];
[x, y, z] = ndgrid(1:size(Ds_NBM_L, 1), 1:size(Ds_NBM_L, 2), 1:size(Ds_NBM_L, 3));
sphere_voxels_L = ((x - center_L(1)).^2 + (y - center_L(2)).^2 + (z - center_L(3)).^2) <= sphere_radius^2;
sphere_voxels_R = ((x - center_R(1)).^2 + (y - center_R(2)).^2 + (z - center_R(3)).^2) <= sphere_radius^2;
spheres_voxels = sphere_voxels_L | sphere_voxels_R;
[xv, yv, zv] = ind2sub(size(spheres_voxels), find(spheres_voxels));
% scatter3(yv, xv, zv, 1, 'filled'); % in this space, x and y need to be swapped

% [xs, ys, zs] = sphere;
% patch(surf2patch(xs*10 + 137, ys*10 + 18, zs*10 + 146), 'FaceColor', 'red', 'EdgeColor', 'none');

% % plot_sphere(5, [98 142 149]);
% plot_sphere(5, [142 98 149]);
% % plot_sphere(5, [128	142	150]);
% plot_sphere(5, [142 128 150]);

axis equal;
view(90, 90);

%% store data for TUS_entry

info = niftiinfo(sprintf('D:/scans/sub-%03d_T1_seg_L.nii.gz', subject));
info.Datatype = 'single';
info.PixelDimensions = head_info.PixelDimensions(1:3); % [1 1 1];
% info.ImageSize = NBM_info.ImageSize;

% scalp data:

fname_scalp = sprintf('M:/Documents/scans/segmentation_results/m2m_sub-%03d/scalp.nii', subject);
% info.PixelDimensions = head_info.PixelDimensions;
% info.Filename = fname_scalp;
% niftiwrite(single(within), fname_scalp, info);
niftiwrite(single(within), fname_scalp);

% target data:

resized_L = imresize3(NBM_L, 'Scale', scaling_factors);
padsize = size(NBM_L) - size(resized_L);% / 2;
% resized_L_padded = padarray(resized_L, floor(padsize), 0, 'pre');
resized_L_padded = padarray(resized_L, ceil(padsize), 0, 'post');
% Ds_NBM_L_padded = smooth3(resized_L_padded);
% NBM_L_isosurface_padded = isosurface(Ds_NBM_L_padded,0.5);
% hiso_NBM_L_padded = patch(NBM_L_isosurface_padded,...
%    'FaceColor',[0.0 0.7 0.0],...
%    'EdgeColor','none');

fname_NBM_L = sprintf('M:/Documents/scans/segmentation_results/m2m_sub-%03d/NBM_L.nii', subject);
% niftiwrite(resized_L_padded, fname_NBM_L, info);
niftiwrite(resized_L_padded, fname_NBM_L);

resized_R = imresize3(NBM_R, 'Scale', scaling_factors);
padsize = size(NBM_R) - size(resized_R);% / 2;
% resized_R_padded = padarray(resized_R, floor(padsize), 0, 'pre');
resized_R_padded = padarray(resized_R, ceil(padsize), 0, 'post');
% Ds_NBM_R_padded = smooth3(resized_R_padded);
% NBM_R_isosurface_padded = isosurface(Ds_NBM_R_padded,0.5);
% hiso_NBM_R_padded = patch(NBM_R_isosurface_padded,...
%    'FaceColor',[0.0 0.7 0.0],...
%    'EdgeColor','none');

fname_NBM_R = sprintf('M:/Documents/scans/segmentation_results/m2m_sub-%03d/NBM_R.nii', subject);
% niftiwrite(resized_R_padded, fname_NBM_R, info);
niftiwrite(resized_R_padded, fname_NBM_R);

% exclude data:

fname_exclude = sprintf('M:/Documents/scans/segmentation_results/m2m_sub-%03d/ears.nii', subject);
% niftiwrite(single(spheres_voxels), fname_exclude, info);
niftiwrite(single(spheres_voxels), fname_exclude);

%% visual confirmation

figure;
hold on;
axis equal;

head = niftiread(fname_scalp);
NBM_L = niftiread(fname_NBM_L);
NBM_R = niftiread(fname_NBM_R);
ears = niftiread(fname_exclude);

Ds = smooth3(head);
skin_isosurface = isosurface(Ds,0.5);
hiso = patch(skin_isosurface,...
   'FaceColor',[0.5 0.5 0.5],...
   'EdgeColor','none',...
   'facealpha',0.05);
Ds = smooth3(NBM_L);
skin_isosurface = isosurface(Ds,0.5);
hiso = patch(skin_isosurface,...
   'FaceColor',[0.0 0.8 0.0],...
   'EdgeColor','none',...
   'facealpha',0.5);
Ds = smooth3(NBM_R);
skin_isosurface = isosurface(Ds,0.5);
hiso = patch(skin_isosurface,...
   'FaceColor',[0.0 0.8 0.0],...
   'EdgeColor','none',...
   'facealpha',0.5);
Ds = smooth3(ears);
skin_isosurface = isosurface(Ds,0.5);
hiso = patch(skin_isosurface,...
   'FaceColor',[1 0 0],...
   'EdgeColor','none',...
   'facealpha',0.5);

[xv, yv, zv] = ind2sub(size(ears), find(ears));
scatter3(yv, xv, zv, 1, 'filled');