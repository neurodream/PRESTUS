
close all; clear; clc;

cd(fileparts(mfilename('fullpath')));
cd ..

% add paths
addpath('functions')
addpath(genpath('toolboxes'))

%% adjust params here:

% base config ("hard" params)
parameters = load_parameters('nico_test_double_acoustic_100mm_config.yaml');

expected_focal_distances_mm = [70 75 80 85 90 95 100 105];
ROI_widths_mm = [5 10 15 20 25 30 35];
i = 1;
figure;
for expected_focal_distance_mm = expected_focal_distances_mm
    for ROI_width_mm = ROI_widths_mm
        [parameters, axial_position, i_axial_oneil] = calculate_transducer_phases(parameters, 1, expected_focal_distance_mm, ROI_width_mm, 100, false);
        subplot(numel(expected_focal_distances_mm), numel(ROI_widths_mm), i);
        plot(axial_position, i_axial_oneil);
        i = i + 1;
    end
end
ylim([0 120]);