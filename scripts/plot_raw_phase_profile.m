% close all;
clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(fileparts(currentFile));
cd(rootpath); % repos/PRESTUS_forked/

addpath('functions');
addpath(genpath('toolboxes'));

% base config ("hard" params)
parameter_fids = {'nico_test_double_acoustic_same_temp0_config', 'IS_double_75mm'};
parameters_fnames = {};
for i = 1:numel(parameter_fids)
    fname = parameter_fids{i};
    parameters_fnames{end+1} = [fname '.yaml'];
    
end
parameters = load_parameters(parameters_fnames{:});

%% 

transducer = parameters.transducers(1);

axial_position = (1:parameters.default_grid_dims(3))*(parameters.grid_step_mm);
axial_position = axial_position';

velocity = 0.02;

p_axial_oneil = focusedAnnulusONeil( ...
    transducer.curv_radius_mm/1e3, ...
    [transducer.Elements_ID_mm; transducer.Elements_OD_mm]/1e3, ...
    repmat(velocity, 1, transducer.n_elements), ...
    [0 159.8 287.37 171.01 142.95 252.9 258.49 234.1 252.3 258.09], ...
    transducer.source_freq_hz, ...
    parameters.medium.water.sound_speed, ...
    parameters.medium.water.density, ...
    (axial_position-0.5)*1e-3 ... % before here was ax_pos
    );

cd(rootpath);

i_axial_oneil = p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;

fig1 = figure('Position', [400, 300, 600, 350]); hold on;

plot(axial_position, i_axial_oneil, 'LineWidth', 2.5, 'Color', [0.6 0.6 0.6]);

load('data/measures/sham_hydrophone_measure.mat');
axial_position_hydrophone = x_data;
p_axis_hydrophone = y_data*1000000;

i_axis_hydrophone = p_axis_hydrophone.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;

plot(axial_position_hydrophone, i_axis_hydrophone, 'LineWidth', 2.5, 'Color', [0.9 0 0]);

load('data/measures/ref_hydrophone_measure.mat');
axial_position_hydrophone_ref = x_data;
p_axis_hydrophone_ref = y_data*1000000;

i_axis_hydrophone_ref = p_axis_hydrophone_ref.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;

plot(axial_position_hydrophone_ref, i_axis_hydrophone_ref, 'LineWidth', 2.5, 'Color', [0.5 0.5 0.9]);

xlabel('axial position (mm)');
ylabel('intensity (W/cm²)');
legend('simulated sham', 'measured sham', 'measured ref');

exportgraphics(gcf,"myplot.png", "BackgroundColor", 'white');