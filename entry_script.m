close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(currentFile);
cd(rootpath); % repos/PRESTUS_forked/

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

%% adjust params here:

subject_id = 5;

% adjust as needed
transducer_labels   = {'L',         'R'};
ID_parts            = {'L+z--l+z',      'R+z--l+z'}; % TODO set label automatically
contralateral       = [false,        true];
dirs                = {'l',         'r'}; % i.e. target side
sham                = [false        false]; % CAREFUL!!!!!! TODO just one; and add to ID!
angles              = [0 -1 0;    0 1 0]; % careful: x and y swapped here!! TODO!
transd_pos_shift    = [0 0 2;       0 0 2];
focus_pos_shift     = [0 0 2;       0 0 2];

% base config ("hard" params)
parameters = load_parameters('nico_test_double_acoustic_100mm_config.yaml');



%%

ID = ['_' strjoin(ID_parts, '_') '']; % TODO make sham a variable!
parameters.results_filename_affix = ID; % TODO this line needed? but also not dangerous

% TODO enable one transducer or more than two transducers
% if numel(transducer_labels)

for i = 1:numel(parameters.transducers)

    parameters.transducers(i).name = transducer_labels{i};
    [parameters, distance] = get_transducer_pos(parameters, subject_id, dirs{i}, i, angles(i,:), transd_pos_shift(i,:), focus_pos_shift(i,:), contralateral(i));
    parameters = calculate_transducer_phases(parameters, i, distance, 15, 100, sham(i));

    % TODO debug delete
    % TODO note!! distance maximized with +30 and -10 for contra- and
    % ipsilateral, respectively
    if contralateral(i)
        parameters = calculate_transducer_phases(parameters, i, distance, 40, 100, sham(i)); % distance + 30
    elseif ~contralateral(i)
        parameters = calculate_transducer_phases(parameters, i, distance, 15, 100, sham(i));
    end

    % parameters = calculate_transducer_phases(parameters, i, focal_distances_mm(i), 15, 100, sham(i));

    parameters = get_simulated_axial_intensity(parameters);
    
    % store the indended parameters for later debugging:
    parameters.transducers(i).optim_params = [];
    % parameters.transducers(i).optim_params.focal_distance_mm = focal_distances_mm(i);
    parameters.transducers(i).optim_params.focal_distance_mm = distance;
    parameters.transducers(i).optim_params.angle = angles(i,:);
    parameters.transducers(i).optim_params.transd_pos_shift = transd_pos_shift(i,:);
    parameters.transducers(i).optim_params.focus_pos_shift = focus_pos_shift(i,:);

end

% adjust the source_amp of contralateral-focused transducer such that
% intensity peaks match
% TODO handle cases of more than 2 transducers
peaks_list1 = findpeaks(parameters.transducers(1).axial_intensity_sim_FW); 
peak1 = peaks_list1(end); % max(parameters.transducers(1).axial_intensity_sim_FW);
peaks_list2 = findpeaks(parameters.transducers(2).axial_intensity_sim_FW); 
peak2 = peaks_list2(end); % max(parameters.transducers(2).axial_intensity_sim_FW);
source_amp1 = parameters.transducers(1).source_amp;
source_amp2 = parameters.transducers(2).source_amp;

if peak1 > peak2
    parameters.transducers(2).source_amp = parameters.transducers(2).source_amp * sqrt(peak1/peak2);
else
    parameters.transducers(1).source_amp = parameters.transducers(1).source_amp * sqrt(peak2/peak1);
end
% TODO handle case of exactly equal peaks

parameters = get_simulated_axial_intensity(parameters);

% figure; hold on;
% plot(parameters.transducers(1).axial_position, parameters.transducers(1).axial_intensity_sim_FW, 'LineWidth', 2.5);
% plot(parameters.transducers(2).axial_position, parameters.transducers(2).axial_intensity_sim_FW, 'LineWidth', 2.5);
% legend(parameters.transducers(1).name, parameters.transducers(2).name);

% plot_transducer_pos(parameters, subject_id, true, false);
% last parameter: where to add position noise ('none' means perfect
% positioning)
update_transducers_and_run(subject_id, parameters, ID, 'none');
% update_transducers_and_run(subject_id, parameters, [ID 'var1'], 'transducer');
% update_transducers_and_run(subject_id, parameters, [ID 'var2'], 'transducer');
% update_transducers_and_run(subject_id, parameters, [ID 'var3'], 'focus');
% update_transducers_and_run(subject_id, parameters, [ID 'var4'], 'focus');
% update_transducers_and_run(subject_id, parameters, [ID 'var5'], 'both');
% update_transducers_and_run(subject_id, parameters, [ID 'var6'], 'both');








%% old

% parameters.transducers(2).name = 'R';
% parameters = calculate_transducer_phases(parameters, 2, 64.691, 15, 100);
% parameters = get_transducer_pos(parameters, subject_id, 'r', 2, [0.2 1 0], [0 0 0], [0 0 0]);

% f1 = parameters.transducers(1).focus_pos_t1_grid;
% f2 = parameters.transducers(2).focus_pos_t1_grid;
% t1 = parameters.transducers(1).pos_t1_grid;
% t2 = parameters.transducers(2).pos_t1_grid;
% 
% parameters.results_filename_affix = [ID 'var1'];
% vt1 = [0, -1, -2]; vf1 = [0, 0, 0]; vt2 = [3, -2, -1]; vf2 = [0, 0, 0];
% parameters.transducers(1).pos_t1_grid 		= [t1(1) + vt1(1),t1(2) + vt1(2),t1(3) + vt1(3)];
% parameters.transducers(1).focus_pos_t1_grid 	= [f1(1) + vf1(1),f1(2) + vf1(2),f1(3) + vf1(3)];
% parameters.transducers(2).pos_t1_grid 		= [t2(1) + vt2(1),t2(2) + vt2(2),t2(3) + vt2(3)];
% parameters.transducers(2).focus_pos_t1_grid 	= [f2(1) + vf2(1),f2(2) + vf2(1),f2(3) + vf2(1)];
% single_subject_pipeline_with_slurm(subject_id, parameters);
% 
% parameters.results_filename_affix = [ID 'var2'];
% vt1 = [0, 2, 0]; vf1 = [0, 0, 0]; vt2 = [2, 0, 1]; vf2 = [0, 0, 0];
% parameters.transducers(1).pos_t1_grid 		= [t1(1) + vt1(1),t1(2) + vt1(2),t1(3) + vt1(3)];
% parameters.transducers(1).focus_pos_t1_grid 	= [f1(1) + vf1(1),f1(2) + vf1(2),f1(3) + vf1(3)];
% parameters.transducers(2).pos_t1_grid 		= [t2(1) + vt2(1),t2(2) + vt2(2),t2(3) + vt2(3)];
% parameters.transducers(2).focus_pos_t1_grid 	= [f2(1) + vf2(1),f2(2) + vf2(1),f2(3) + vf2(1)];
% single_subject_pipeline_with_slurm(subject_id, parameters);
% 
% parameters.results_filename_affix = [ID 'var3'];
% vt1 = [0, 0, 0]; vf1 = [-1, 3, 0]; vt2 = [0, 0, 0]; vf2 = [-3, 2, 2];
% parameters.transducers(1).pos_t1_grid 		= [t1(1) + vt1(1),t1(2) + vt1(2),t1(3) + vt1(3)];
% parameters.transducers(1).focus_pos_t1_grid 	= [f1(1) + vf1(1),f1(2) + vf1(2),f1(3) + vf1(3)];
% parameters.transducers(2).pos_t1_grid 		= [t2(1) + vt2(1),t2(2) + vt2(2),t2(3) + vt2(3)];
% parameters.transducers(2).focus_pos_t1_grid 	= [f2(1) + vf2(1),f2(2) + vf2(1),f2(3) + vf2(1)];
% single_subject_pipeline_with_slurm(subject_id, parameters);
% 
% parameters.results_filename_affix = [ID 'var4'];
% vt1 = [0, 0, 0]; vf1 = [1, -3, 1]; vt2 = [0, 0, 0]; vf2 = [-3, -3, 2];
% parameters.transducers(1).pos_t1_grid 		= [t1(1) + vt1(1),t1(2) + vt1(2),t1(3) + vt1(3)];
% parameters.transducers(1).focus_pos_t1_grid 	= [f1(1) + vf1(1),f1(2) + vf1(2),f1(3) + vf1(3)];
% parameters.transducers(2).pos_t1_grid 		= [t2(1) + vt2(1),t2(2) + vt2(2),t2(3) + vt2(3)];
% parameters.transducers(2).focus_pos_t1_grid 	= [f2(1) + vf2(1),f2(2) + vf2(1),f2(3) + vf2(1)];
% single_subject_pipeline_with_slurm(subject_id, parameters);
% 
% parameters.results_filename_affix = [ID 'var5'];
% vt1 = [3, -3, -3]; vf1 = [-1, 0, -1]; vt2 = [-1, 3, 1]; vf2 = [-1, 1, 1];
% parameters.transducers(1).pos_t1_grid 		= [t1(1) + vt1(1),t1(2) + vt1(2),t1(3) + vt1(3)];
% parameters.transducers(1).focus_pos_t1_grid 	= [f1(1) + vf1(1),f1(2) + vf1(2),f1(3) + vf1(3)];
% parameters.transducers(2).pos_t1_grid 		= [t2(1) + vt2(1),t2(2) + vt2(2),t2(3) + vt2(3)];
% parameters.transducers(2).focus_pos_t1_grid 	= [f2(1) + vf2(1),f2(2) + vf2(1),f2(3) + vf2(1)];
% single_subject_pipeline_with_slurm(subject_id, parameters);
% 
% parameters.results_filename_affix = [ID 'var6'];
% vt1 = [-3, -1, -1]; vf1 = [0, 0, 3]; vt2 = [-2, 2, 0]; vf2 = [0, -1, -1];
% parameters.transducers(1).pos_t1_grid 		= [t1(1) + vt1(1),t1(2) + vt1(2),t1(3) + vt1(3)];
% parameters.transducers(1).focus_pos_t1_grid 	= [f1(1) + vf1(1),f1(2) + vf1(2),f1(3) + vf1(3)];
% parameters.transducers(2).pos_t1_grid 		= [t2(1) + vt2(1),t2(2) + vt2(2),t2(3) + vt2(3)];
% parameters.transducers(2).focus_pos_t1_grid 	= [f2(1) + vf2(1),f2(2) + vf2(1),f2(3) + vf2(1)];
% single_subject_pipeline_with_slurm(subject_id, parameters);