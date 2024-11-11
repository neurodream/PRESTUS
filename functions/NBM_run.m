function [] = NBM_run(subject_id, focus_side, is_sham, z_shift, extra_ID_suffix, with_variations)

% focus_side: L, R, BL

transducer_labels   = {'L',         'R'};
if strcmp(focus_side, 'L')
    contralateral       = [false,        true];
    focus_sides_add_to_ID = {'l', 'l'};
elseif strcmp(focus_side, 'R')
    contralateral       = [true,        false];
    focus_sides_add_to_ID = {'r', 'r'};
else
    contralateral       = [false,        false];
    focus_sides_add_to_ID = {'l', 'r'};
end % TODO add middle focus

if z_shift > 0
    z_shift_part_ID = '+z';
else
    z_shift_part_ID = '';
end

if is_sham
    sham_part_ID = '-';
else
    sham_part_ID = '--';
end

ID_part = ['L' z_shift_part_ID sham_part_ID focus_sides_add_to_ID{1} '_' 'R' z_shift_part_ID sham_part_ID focus_sides_add_to_ID{2}];
%= {'L+z--l+z',      'R+z--l+z'}; % TODO set label automatically
dirs                = {'l',         'r'}; % i.e. target side (note: just for readout of the coordinates that determine the beam vector! the _distance_ will determine ipsi vs. contralateral
sham                = [is_sham        is_sham]; % CAREFUL!!!!!! TODO just one; and add to ID!
angles              = [0 -1 0;    0 1 0]; % careful: x and y swapped here!! TODO!
transd_pos_shift    = [0 0 z_shift;       0 0 z_shift];
focus_pos_shift     = [0 0 z_shift;       0 0 z_shift];

% base config ("hard" params)
parameters = load_parameters('nico_test_double_acoustic_100mm_config.yaml');

ID = [ID_part '_' extra_ID_suffix]; % TODO make sham a variable!
parameters.results_filename_affix = ID; % TODO this line needed? but also not dangerous

for i = 1:numel(parameters.transducers)

    parameters.transducers(i).name = transducer_labels{i};
    [parameters, distance] = get_transducer_pos(parameters, subject_id, dirs{i}, i, angles(i,:), transd_pos_shift(i,:), focus_pos_shift(i,:), contralateral(i));
    % parameters = calculate_transducer_phases(parameters, i, distance, 15, 100, sham(i));

    % TODO figure out which optimization works best
    if contralateral(i)
        parameters = calculate_transducer_phases(parameters, i, distance, 40, 100, sham(i)); % distance + 30
    elseif ~contralateral(i)
        parameters = calculate_transducer_phases(parameters, i, distance, 15, 100, sham(i));
    end

    % parameters = calculate_transducer_phases(parameters, i, focal_distances_mm(i), 15, 100, sham(i));
    
    % store the indended parameters for later debugging:
    parameters.transducers(i).optim_params = [];
    parameters.transducers(i).optim_params.focal_distance_mm = distance;
    parameters.transducers(i).optim_params.angle = angles(i,:);
    parameters.transducers(i).optim_params.transd_pos_shift = transd_pos_shift(i,:);
    parameters.transducers(i).optim_params.focus_pos_shift = focus_pos_shift(i,:);

end

% add field of free water axial intensity to strucuts
parameters = get_simulated_axial_intensity(parameters);

update_transducers_and_run(subject_id, parameters, ID, 'none');
if with_variations
    update_transducers_and_run(subject_id, parameters, [ID 'var1'], 'transducer');
    % update_transducers_and_run(subject_id, parameters, [ID 'var2'], 'transducer');
    % update_transducers_and_run(subject_id, parameters, [ID 'var3'], 'focus');
    % update_transducers_and_run(subject_id, parameters, [ID 'var4'], 'focus');
    % update_transducers_and_run(subject_id, parameters, [ID 'var5'], 'both');
    % update_transducers_and_run(subject_id, parameters, [ID 'var6'], 'both');
end

end