%% set the layers
% of interest: skin, skull, brain, brain_inner, target (x2), ROItarget (x2)
% nice to have: inner brain minus target ROI

close all; clear; clc;

cd /home/sleep/nicade/Documents/repos/PRESTUS_forked/; % repos/PRESTUS_forked/

addpath('/home/sleep/nicade/Documents/repos/PRESTUS_forked/functions');

layers = niftiread('/home/sleep/nicade/Documents/scans/segmentation_results/m2m_sub-001/final_tissues.nii.gz');
info   = niftiinfo('/home/sleep/nicade/Documents/scans/segmentation_results/m2m_sub-001/final_tissues.nii.gz');

% simnibs-segmented tissues
skull = layers == 7 | layers == 8;
brain = layers == 1 | layers == 2;
skin  = layers == 5;

% inner brain
SE = strel('sphere', 6);

skull_dilated = imdilate(skull, SE);
brain_eroded  = imerode(brain, SE); % not optimal: erodes around the gyri

brain_inner = brain & ~skull_dilated;

% targets (careful: hardcoded for subject 1)
targetL = [96 144 149];
targetR = [130 144 150];

r = round(5/mean(info.PixelDimensions));
% Equation of the sphere: (x - px)^2 + (y - py)^2 + (z - pz)^2 <= r^2
[x, y, z] = ndgrid(1:size(layers,1), 1:size(layers,2), 1:size(layers,3));
ROItarget_L = (x - targetL(1)).^2 + (y - targetL(2)).^2 + (z - targetL(3)).^2 <= r^2;
ROItarget_R = (x - targetR(1)).^2 + (y - targetR(2)).^2 + (z - targetR(3)).^2 <= r^2;

%% load in the results and readout

filename = '_distance_equal.csv';

% List of IDs to loop over
IDs = {
    '_L--l_R--r_', ...
    '_L--l_R--r_var1_', ...
    '_L--l_R--r_var2_', ...
    '_L--l_R--r_var3_', ...
    '_L--l_R--r_var4_', ...
    '_L--l_R--r_var5_' ...
    % '_L--l_R--r_var6_' ...
    };
IDs = {'_L+y--l_R+y--r_'};

sub_id = 5;

for id_idx = 1:numel(IDs)
    ID = IDs{id_idx};  % Current ID

    [isppa, p, mi] = acoustic_mat_to_nifti(sub_id, ID, sprintf('/home/sleep/nicade/Documents/scans/sim_outputs/sub-%03d', sub_id), true);
    p = p/1e6; % in megapascal
    i = isppa;
    
    % skin, skull, brain, brain_inner, target (x2), ROItarget (x2)
    
    p_skin = max(p(skin));
    p_skull = max(p(skull));
    p_brain = max(p(brain));
    p_brain_inner = max(p(brain_inner));
    p_targetL = p(targetL(1), targetL(2), targetL(3));
    p_targetR = p(targetR(1), targetR(2), targetR(3));
    p_ROItarget_L_max = max(p(ROItarget_L));
    p_ROItarget_R_max = max(p(ROItarget_R));
    p_ROItarget_L_avg = mean(p(ROItarget_L));
    p_ROItarget_R_avg = mean(p(ROItarget_R));
    
    i_skin = max(i(skin));
    i_skull = max(i(skull));
    i_brain = max(i(brain));
    i_brain_inner = max(i(brain_inner));
    i_targetL = i(targetL(1), targetL(2), targetL(3));
    i_targetR = i(targetR(1), targetR(2), targetR(3));
    i_ROItarget_L_max = max(i(ROItarget_L));
    i_ROItarget_R_max = max(i(ROItarget_R));
    i_ROItarget_L_avg = mean(i(ROItarget_L));
    i_ROItarget_R_avg = mean(i(ROItarget_R));

    mi_skin = max(mi(skin));
    mi_skull = max(mi(skull));
    mi_brain = max(mi(brain));
    mi_brain_inner = max(mi(brain_inner));
    mi_targetL = mi(targetL(1), targetL(2), targetL(3));
    mi_targetR = mi(targetR(1), targetR(2), targetR(3));
    mi_ROItarget_L_max = max(mi(ROItarget_L));
    mi_ROItarget_R_max = max(mi(ROItarget_R));
    mi_ROItarget_L_avg = mean(mi(ROItarget_L));
    mi_ROItarget_R_avg = mean(mi(ROItarget_R));

    % Create a table for the current ID
    T = table( ...
        [p_skin; i_skin, mi_skin], ...
        [p_skull; i_skull, mi_skull], ...
        [p_brain; i_brain, mi_brain], ...
        [p_brain_inner; i_brain_inner, mi_brain_inner], ...
        [p_targetL; i_targetL, mi_targetL], ...
        [p_targetR; i_targetR, mi_targetR], ...
        [p_ROItarget_L_max; i_ROItarget_L_max, mi_ROItarget_L_max], ...
        [p_ROItarget_R_max; i_ROItarget_R_max, mi_ROItarget_R_max], ...
        [p_ROItarget_L_avg; i_ROItarget_L_avg, mi_ROItarget_L_avg], ...
        [p_ROItarget_R_avg; i_ROItarget_R_avg, mi_ROItarget_R_avg], ...
        'VariableNames', {'skin', 'skull', 'brain', 'brain_inner', 'targetL', 'targetR', 'ROItarget_L_max', 'ROItarget_R_max', 'ROItarget_L_avg', 'ROItarget_R_avg'}, ...
        'RowNames', {['pressure_' ID], ['intensity_' ID], ['mechanicalindex_' ID]} ...
    );

    % Add the sub_id column
    T.sub_id = repmat(sub_id, height(T), 1);
    
    % Reorder columns to put sub_id first if needed
    T = movevars(T, 'sub_id', 'Before', 'skin');
    
    % Check if the file exists
    if isfile(filename)
        % Read the existing table
        T_existing = readtable(filename, 'ReadRowNames', true);
        
        % Concatenate the new data to the existing table
        T_combined = [T_existing; T];
        
        % Write the combined table back to the file
        writetable(T_combined, filename, 'WriteRowNames', true);
    else
        % If the file does not exist, write the new table
        writetable(T, filename, 'WriteRowNames', true);
    end
end
