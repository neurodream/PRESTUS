heating = niftiread('D:\sim_outputs\sub-003\sub-003_final_heating_orig_coordsingletransducer.nii.gz');
CEM43 = niftiread('D:\sim_outputs\sub-003\sub-003_final_CEM43_orig_coordsingletransducer.nii.gz');
isppa = niftiread('D:\sim_outputs\sub-003\sub-003_final_isppa_orig_coordsingletransducer.nii.gz');
MI = niftiread('D:\sim_outputs\sub-003\sub-003_final_MI_orig_coordsingletransducer.nii.gz');
pressure = niftiread('D:\sim_outputs\sub-003\sub-003_final_pressure_orig_coordsingletransducer.nii.gz');

plot(sort(heating(:)))
yline(37, 'LineStyle', '--')

%% get the heating values deeper in the brain

% first dilate the skull, then subtract it from brain

layers = niftiread('M:\Documents\scans\segmentation_results\m2m_sub-003\final_tissues.nii.gz');

skull = layers == 7 | layers == 8;
brain = layers == 1 | layers == 2;
skin  = layers == 5;

SE = strel('sphere', 6);

skull_dilated = imdilate(skull, SE);
brain_eroded  = imerode(brain, SE); % not optimal: erodes around the gyri

brain_inner = brain & ~skull_dilated;

target_coords = [98 142 149];

%% test

volshow(brain)
volshow(brain_inner)

%% get the temperature values

heating_brain       = max(heating(brain));
heating_skull       = max(heating(skull));
heating_skin        = max(heating(skin));
heating_target      = heating(target_coords(1),target_coords(2),target_coords(3));


isppa_brain_inner = max(isppa(brain_inner))
pressure_brain_inner = max(pressure(brain_inner))
MI_brain_inner = max(MI(brain_inner))
heating_brain_inner = max(heating(brain_inner))
CEM43_brain_inner = max(CEM43(brain_inner))


%%

isppa_target = isppa(target_coords(1),target_coords(2),target_coords(3))
pressure_target = pressure(target_coords(1),target_coords(2),target_coords(3))
MI_target = MI(target_coords(1),target_coords(2),target_coords(3))
heating_target = heating(target_coords(1),target_coords(2),target_coords(3))
CEM43_target = CEM43(target_coords(1),target_coords(2),target_coords(3))


% good news: rather low inner brain heating!







%% not really needed: display the heating in 3D

% heating_thresholded = heating;
% 
% % TODO threshold should be based on 
% heating_thresholded(heating_thresholded <= 37.1) = 0;
% 
% % remove borders
% border_thickness = 20;
% 
% heating_thresholded(:,1:border_thickness,:) = zeros(size(heating_thresholded,1), border_thickness, size(heating_thresholded,3));
% heating_thresholded(1:border_thickness,:,:) = zeros(border_thickness, size(heating_thresholded,2), size(heating_thresholded,3));
% heating_thresholded(:,:,1:border_thickness) = zeros(size(heating_thresholded,1), size(heating_thresholded,2), border_thickness);
% 
% heating_thresholded(:,end-border_thickness+1:end,:) = zeros(size(heating_thresholded,1), border_thickness, size(heating_thresholded,3));
% heating_thresholded(end-border_thickness+1:end,:,:) = zeros(border_thickness, size(heating_thresholded,2), size(heating_thresholded,3));
% heating_thresholded(:,:,end-border_thickness+1:end) = zeros(size(heating_thresholded,1), size(heating_thresholded,2), border_thickness);
% 
% volshow(heating_thresholded)