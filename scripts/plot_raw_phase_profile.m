% close all;
clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(fileparts(currentFile));
cd(rootpath); % repos/PRESTUS_forked/

addpath('functions');
addpath(genpath('toolboxes'));

% base config ("hard" params)
parameter_fids = {'nico_test_double_acoustic_same_temp0_config', 'IS_L_100mm_R_75mm'}; % IS_double_100mm
parameters_fnames = {};
for i = 1:numel(parameter_fids)
    fname = parameter_fids{i};
    parameters_fnames{end+1} = [fname '.yaml'];
    
end
parameters = load_parameters(parameters_fnames{:});
parameters_win_paths = load_parameters('default_paths_windows.yaml');
parameters.seg_path = parameters_win_paths.seg_path;
parameters.data_path = parameters_win_paths.data_path;

%% 

transducer = parameters.transducers(1);

axial_position = (1:parameters.default_grid_dims(3))*(parameters.grid_step_mm);
axial_position = axial_position';

folders = dir(fullfile(parameters.data_path, 'sim_outputs', 'sub-*'));
folders = folders(~ismember({folders.name}, {'.', '..'})); % TODO ? necessary?
for i = 1:numel(folders)
    sub_id = sscanf(folders(i).name, 'sub-%d');
    files = dir(fullfile(parameters.data_path, 'sim_outputs', folders(i).name, '*parametersL*-*_R*-*_it2_*'));
    % if numel(files) > 2
    %     files = files(~contains({files.name}, 'parallel'));
    % end
    % for j = 1:numel(files)
    %     file = fullfile(files(j).folder, files(j).name);
        file = fullfile(files(1).folder, files(1).name);
        parameters = load(file);
        parameters = parameters.parameters;
    % end

    fig1 = figure('Position', [400, 300, 600, 350]); hold on;

    for transducer = parameters.transducers
        velocity = 0.15;
        
        phases = transducer.source_phase_rad;
        
        if strcmp(transducer.name, 'L')
            % phases = deg2rad([163.85, 148.61, 133.58, 117.97, 101.91, 85.67, 68.80, 51.59, 34.31, 16.72]); % phases for 96 mm wrt
        elseif strcmp(transducer.name, 'R')
            % phases = deg2rad([232.81, 230.72, 228.65, 226.51, 224.31, 222.10, 219.80, 217.46, 215.11, 212.73]); % phases for 67 mm wrt
            velocity = velocity*0.75;
        end
        
        
        p_axial_oneil = focusedAnnulusONeil( ...
            transducer.curv_radius_mm/1e3, ...
            [transducer.Elements_ID_mm; transducer.Elements_OD_mm]/1e3, ...
            repmat(velocity, 1, transducer.n_elements), ...
            phases, ...
            transducer.source_freq_hz, ...
            parameters.medium.water.sound_speed, ...
            parameters.medium.water.density, ...
            (axial_position-0.5)*1e-3 ... % before here was ax_pos
            );
        
        cd(rootpath);
        
        i_axial_oneil = p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;
        
        plot(axial_position, i_axial_oneil, 'LineWidth', 2.5);
        
        % load('data/measures/sham_hydrophone_measure.mat');
        % axial_position_hydrophone = x_data;
        % p_axis_hydrophone = y_data*1000000;
        % 
        % i_axis_hydrophone = p_axis_hydrophone.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;
        % 
        % plot(axial_position_hydrophone, i_axis_hydrophone, 'LineWidth', 2.5, 'Color', [0.9 0 0]);
        % 
        % load('data/measures/ref_hydrophone_measure.mat');
        % axial_position_hydrophone_ref = x_data;
        % p_axis_hydrophone_ref = y_data*1000000;
        % 
        % i_axis_hydrophone_ref = p_axis_hydrophone_ref.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;
        % 
        % plot(axial_position_hydrophone_ref, i_axis_hydrophone_ref, 'LineWidth', 2.5, 'Color', [0.5 0.5 0.9]);
            

    end

    xlabel('axial position (mm)');
    ylabel('intensity (W/cm²)');

    legend({parameters.transducers.name});
    
    % exportgraphics(gcf,"myplot.png", "BackgroundColor", 'white');

end

