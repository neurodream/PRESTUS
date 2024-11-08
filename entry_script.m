close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(currentFile);
cd(rootpath); % repos/PRESTUS_forked/

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

% The parameters of the individual calls are, from left to right:
% sbj_num, 
% focus_side (L or R for unilateral, LR for bilateral), 
% sham or not (true or false), 
% how many z voxels for under/overshoot correction,
% a string abbreviation which describes more precisely (can also be left blank: '')
% and if errors in spatial precision are to be simulated as well (true or false)
NBM_run(1, 'L', false, 3, 'all_phases', false);
NBM_run(2, 'L', false, 4, 'all_phases', false);
NBM_run(3, 'L', false, 3, 'all_phases', false);
NBM_run(4, 'L', false, 2, 'all_phases', false);
NBM_run(5, 'L', false, 2, 'all_phases', false);
NBM_run(6, 'L', false, 0, 'all_phases', false);

NBM_run(1, 'R', false, 3, 'all_phases', false);
NBM_run(2, 'R', false, 4, 'all_phases', false);
NBM_run(3, 'R', false, 3, 'all_phases', false);
NBM_run(4, 'R', false, 2, 'all_phases', false);
NBM_run(5, 'R', false, 2, 'all_phases', false);
NBM_run(6, 'R', false, 0, 'all_phases', false);

% TODO add to single_subject_pipeline, but careful: assumes a stored nifti
% file
% plot_transducer_pos(parameters, sbj_ID, plot_scalp, plot_skull, plot_intensity, save)

% 
% %% adjust params here:
% 
% subject_id = 6;
% 
% % adjust as needed
% % TODO also give transducer hardware specs
% transducer_labels   = {'L',         'R'};
% ID_parts            = {'L+z--l+z',      'R+z--l+z'}; % TODO set label automatically
% contralateral       = [false,        true];
% dirs                = {'l',         'r'}; % i.e. target side
% sham                = [false        false]; % CAREFUL!!!!!! TODO just one; and add to ID!
% angles              = [0 -1 0;    0 1 0]; % careful: x and y swapped here!! TODO!
% transd_pos_shift    = [0 0 2;       0 0 2];
% focus_pos_shift     = [0 0 2;       0 0 2];
% 
% % transducer_labels   = {'L'};
% % ID_parts            = {'L--r'}; % TODO set label automatically
% % contralateral       = [true];
% % dirs                = {'r'}; % i.e. target side
% % sham                = [false]; % CAREFUL!!!!!! TODO just one; and add to ID!
% % angles              = [0 -1 0]; % careful: x and y swapped here!! TODO!
% % transd_pos_shift    = [0 0 0];
% % focus_pos_shift     = [0 0 0];
% 
% % base config ("hard" params)
% parameters = load_parameters('nico_test_double_acoustic_100mm_config.yaml');
% 
% 
% 
% %%
% 
% ID = ['_' strjoin(ID_parts, '_') '']; % TODO make sham a variable!
% parameters.results_filename_affix = ID; % TODO this line needed? but also not dangerous
% 
% % TODO enable one transducer or more than two transducers
% if numel(transducer_labels) == 1 %#ok<ISCL>
%     parameters.transducers = parameters.transducers(1);
% end
% 
% 
% for i = 1:numel(parameters.transducers)
% 
%     parameters.transducers(i).name = transducer_labels{i};
%     [parameters, distance] = get_transducer_pos(parameters, subject_id, dirs{i}, i, angles(i,:), transd_pos_shift(i,:), focus_pos_shift(i,:), contralateral(i));
%     % parameters = calculate_transducer_phases(parameters, i, distance, 15, 100, sham(i));
% 
%     % TODO figure out which optimization works best
%     if contralateral(i)
%         parameters = calculate_transducer_phases(parameters, i, distance, 40, 100, sham(i)); % distance + 30
%     elseif ~contralateral(i)
%         parameters = calculate_transducer_phases(parameters, i, distance, 15, 100, sham(i));
%     end
% 
%     % parameters = calculate_transducer_phases(parameters, i, focal_distances_mm(i), 15, 100, sham(i));
% 
%     % store the indended parameters for later debugging:
%     parameters.transducers(i).optim_params = [];
%     parameters.transducers(i).optim_params.focal_distance_mm = distance;
%     parameters.transducers(i).optim_params.angle = angles(i,:);
%     parameters.transducers(i).optim_params.transd_pos_shift = transd_pos_shift(i,:);
%     parameters.transducers(i).optim_params.focus_pos_shift = focus_pos_shift(i,:);
% 
% end
% 
% % TODO check: this should not be necessary anymore
% 
% % % adjust the source_amp of contralateral-focused transducer such that
% % % intensity peaks match
% % % TODO handle cases of more than 2 transducers
% % peaks_list1 = findpeaks(parameters.transducers(1).axial_intensity_sim_FW); 
% % peak1 = peaks_list1(end); % max(parameters.transducers(1).axial_intensity_sim_FW);
% % peaks_list2 = findpeaks(parameters.transducers(2).axial_intensity_sim_FW); 
% % peak2 = peaks_list2(end); % max(parameters.transducers(2).axial_intensity_sim_FW);
% % source_amp1 = parameters.transducers(1).source_amp;
% % source_amp2 = parameters.transducers(2).source_amp;
% % 
% % if peak1 > peak2
% %     parameters.transducers(2).source_amp = parameters.transducers(2).source_amp * sqrt(peak1/peak2);
% % else
% %     parameters.transducers(1).source_amp = parameters.transducers(1).source_amp * sqrt(peak2/peak1);
% % end
% % % TODO handle case of exactly equal peaks
% 
% % add field of free water axial intensity to strucuts
% parameters = get_simulated_axial_intensity(parameters);
% 
% % plot_transducer_pos(parameters, subject_id, true, false);
% % last parameter: where to add position noise ('none' means perfect
% % positioning)
% update_transducers_and_run(subject_id, parameters, ID, 'none');
% % update_transducers_and_run(subject_id, parameters, [ID 'var1'], 'transducer');
% % update_transducers_and_run(subject_id, parameters, [ID 'var2'], 'transducer');
% % update_transducers_and_run(subject_id, parameters, [ID 'var3'], 'focus');
% % update_transducers_and_run(subject_id, parameters, [ID 'var4'], 'focus');
% % update_transducers_and_run(subject_id, parameters, [ID 'var5'], 'both');
% % update_transducers_and_run(subject_id, parameters, [ID 'var6'], 'both');