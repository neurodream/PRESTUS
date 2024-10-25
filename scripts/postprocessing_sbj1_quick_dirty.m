close all; clear; clc;

cd M:/Documents/; % repos/PRESTUS_forked/

% heating     = niftiread('scans/sim_outputs/sub-001/sub-001_final_isppa_orig_coord_crossbeam01_.nii.gz');
% CEM43       = niftiread('scans/sim_outputs/sub-001/sub-001_final_isppa_orig_coord_crossbeam01_.nii.gz');
isppa       = niftiread('scans/sim_outputs/sub-001/sub-001_final_isppa_orig_coord_crossbeam01_.nii.gz');
MI          = niftiread('scans/sim_outputs/sub-001/sub-001_final_MI_orig_coord_crossbeam01_.nii.gz');
pressure    = niftiread('scans/sim_outputs/sub-001/sub-001_final_pressure_orig_coord_crossbeam01_.nii.gz');

%%

layers = niftiread('scans/segmentation_results/m2m_sub-001/final_tissues.nii.gz');

skull = layers == 7 | layers == 8;
brain = layers == 1 | layers == 2;
skin  = layers == 5;

SE = strel('sphere', 6);

skull_dilated = imdilate(skull, SE);
brain_eroded  = imerode(brain, SE); % not optimal: erodes around the gyri

brain_inner = brain & ~skull_dilated;

%%
tL = [96 144 149];
tR = [130 144 150];

%% test

volshow(brain)
volshow(brain_inner)

%% at target 1

isppa_target    = isppa(tL(1),tL(2),tL(3))
pressure_target = pressure(tL(1),tL(2),tL(3))/1000000
% MI_target       = MI(tL(1),tL(2),tL(3))

%% at target 2

isppa_target    = isppa(tR(1),tR(2),tR(3))
pressure_target = pressure(tR(1),tR(2),tR(3))/1000000
% MI_target       = MI(tR(1),tR(2),tR(3))

%% max in skin

isppa_skin = isppa(skin);
isppa_skin_max = max(isppa_skin(:))
pressure_skin = pressure(skin);
pressure_skin_max = max(pressure_skin(:))/1000000

%% max in skull

isppa_skull = isppa(skull);
isppa_skull_max = max(isppa_skull(:))
pressure_skull = pressure(skull);
pressure_skull_max = max(pressure_skull(:))/1000000

%% max in brain

isppa_brain = isppa(brain);
isppa_brain_max = max(isppa_brain(:))
pressure_brain = pressure(brain);
pressure_brain_max = max(pressure_brain(:))/1000000

%% max in inner brain

isppa_brain_inner = isppa(brain_inner);
isppa_brain_inner_max = max(isppa_brain_inner(:))
pressure_brain_inner = pressure(brain_inner);
pressure_brain_inner_max = max(pressure_brain_inner(:))/1000000