function [source, source_labels, transducer_pars] = setup_source_two_transducers(parameters, kgrid, trans_pos, focus_pos)

% Function that creates a source for kwave based on the parameters structure, kwave grid, the transducer coordinates
% in the computational grid, and the geometric focus coordinates. Also returns
% a 3d array with different numbers indicating different transducer
% elements and the updated transducer parameters for two transducers.

% Define transducer parameters for both transducers
transducer_pars_all = {parameters.transducer, parameters.transducer2};
cw_signals_all = cell(1, numel(parameters.transducers));

% TODO think if the source labels could be optimized;
% maybe two matrices, one storing transducer number and the other storing
% element number

% Convert from [mm] to [grid points] and round diameters to nearest odd integer for both transducers
for transducer_idx = 1:numel(parameters.transducers)
    parameters.transducers(transducer_idx).Elements_OD = 2 * floor(parameters.transducers(transducer_idx).Elements_OD_mm / parameters.grid_step_mm / 2) + 1; % [grid points]
    parameters.transducers(transducer_idx).Elements_ID = 2 * floor(parameters.transducers(transducer_idx).Elements_ID_mm / parameters.grid_step_mm / 2) + 1; % [grid points]
    parameters.transducers(transducer_idx).Elements_ID(parameters.transducers(transducer_idx).Elements_ID_mm == 0) = 0;
    parameters.transducers(transducer_idx).radius_grid = round(parameters.transducers(transducer_idx).curv_radius_mm / parameters.grid_step_mm); % [grid points]

    % Create the CW signals for both transducers
    cw_signals_all{transducer_idx} = createCWSignals(kgrid.t_array, parameters.transducers(transducer_idx).source_freq_hz, parameters.transducers(transducer_idx).source_amp, parameters.transducers(transducer_idx).source_phase_rad);
end

% Initialize source and labels
source = struct();
grid_dims = parameters.grid_dims;
transducer_mask = zeros(grid_dims);
source_labels = [];
source_labels.transducer = zeros(grid_dims);
source_labels.element = zeros(grid_dims);

% Loop over both transducers to create the source and mask
for transducer_idx = 1:numel(parameters.transducers)
    trans_pos_current = parameters.transducers(transducer_idx).trans_pos_final;
    focus_pos_current = parameters.transducers(transducer_idx).focus_pos_final;
    transducer_pars = parameters.transducers(transducer_idx);
    
    % Create element bowls/arcs one by one
    for el_i = 1:transducer_pars.n_elements
        if parameters.n_sim_dims == 3
            bowl = makeBowl(grid_dims, trans_pos_current, transducer_pars.radius_grid, transducer_pars.Elements_OD(el_i), focus_pos_current(transducer_idx,:));
        else
            bowl = makeArc(grid_dims, trans_pos_current, transducer_pars.radius_grid, transducer_pars.Elements_OD(el_i), focus_pos_current(transducer_idx,:));
        end
        if transducer_pars.Elements_ID(el_i) > 0
            if parameters.n_sim_dims == 3
                bowl = bowl - makeBowl(grid_dims, trans_pos_current, transducer_pars.radius_grid, transducer_pars.Elements_ID(el_i), focus_pos_current(transducer_idx,:));
            else
                bowl = bowl - makeArc(grid_dims, trans_pos_current, transducer_pars.radius_grid, transducer_pars.Elements_ID(el_i), focus_pos_current(transducer_idx,:));
            end
        end

        % Define the binary source mask and label the elements
        transducer_mask = transducer_mask + bowl;
        source_labels.transducer = bowl * transducer_idx;
        source_labels.element = bowl * el_i;
        % source_labels + (el_i + (transducer_idx - 1) * transducer_pars.n_elements) * bowl;
    end
end

% Ensure each source point has an assigned time series of the source signal
p_mask_source_p = source_labels.element;
p_mask_source_p(p_mask_source_p(:,:,:) == 0) = [];
p_mask_source_p = reshape(p_mask_source_p, [], 1);

% Initialize source signal
source.p = zeros(length(p_mask_source_p), length(cw_signals_all{1}));

for ii = 1:length(p_mask_source_p)
    [x, y, z] = ind2sub(size(p_mask_source_p), ii);
    transducer_idx = source_labels.transducer(x,y,z);
    element_idx = source_labels.element(x,y,z);
    source.p(ii, :) = cw_signals_all{transducer_idx}(element_idx, :);
    % if p_mask_source_p(ii) <= parameters.transducers(1).n_elements
    %     source.p(ii, :) = cw_signals_all{1}(p_mask_source_p(ii), :);
    % else
    %     mod_index = mod(p_mask_source_p(ii) - 1, parameters.transducers(1).n_elements) + 1;
    %     source.p(ii, :) = cw_signals_all{2}(mod_index, :);
    % end
end

% Assign the binary source mask
source.p_mask = transducer_mask;

% If using kWaveArray, handle it here
% TODO multiple transducer setup not fully integrated here, won't work!
if parameters.use_kWaveArray ~= 0
    disp('Setting up kWaveArray (might take a bit of time)');
    % create empty kWaveArray
    karray = kWaveArray('BLITolerance', 0.1, 'UpsamplingRate', 10, 'BLIType', 'sinc');

    n_elements_total = 0;

    % Loop over transducers to add to kWaveArray
    for transducer_idx = 1:numel(parameters.transducers)

        trans_pos = parameters.transducers(transducer_idx).trans_pos_final;
        transducer_pars = parameters.transducers(transducer_idx);

        n_elements_total = n_elements_total + transducer_pars.n_elements;
        
        % Add bowl shaped element
        karray.addAnnularArray([kgrid.x_vec(trans_pos(1)) kgrid.y_vec(trans_pos(2)) kgrid.z_vec(trans_pos(3))], transducer_pars.curv_radius_mm * 1e-3, [transducer_pars.Elements_ID_mm; transducer_pars.Elements_OD_mm] * 1e-3, [kgrid.x_vec(focus_pos(transducer_idx,1)) kgrid.y_vec(focus_pos(transducer_idx,2)) kgrid.z_vec(focus_pos(transducer_idx,3))]);
    end

    % Calculate grid weights for both transducers
    grid_weights_4d = zeros(n_elements_total, kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1));
    
    for transducer_idx = 1:numel(parameters.transducers)
        for el_i = 1:parameters.transducers(transducer_idx).n_elements
            idx = el_i + (transducer_idx - 1) * parameters.transducers(1).n_elements; % TODO check if correct... probably not!
            fprintf('Computing weights for transducer %i, element %i...', transducer_idx, el_i); 
            grid_weights_4d(idx,:,:,:) = karray.getElementGridWeights(kgrid, idx);                
            fprintf(' done\n');
        end
    end
    
    % Get the binary mask
    binary_mask = squeeze(sum(grid_weights_4d, 1)) ~= 0;
    
    % Number of time points in the signal
    Nt = size(cw_signals_all{1}, 2);
    
    mask_ind = find(binary_mask);
    num_source_points = sum(binary_mask(:));
     
    % Initialize the source signal
    distributed_source_signal = zeros(num_source_points, Nt);
    
    source_labels = zeros(kgrid.Nx, kgrid.Ny, max(kgrid.Nz, 1));
    if canUseGPU
        cw_signals_all{1} = gpuArray(cw_signals_all{1});
        cw_signals_all{2} = gpuArray(cw_signals_all{2});
        distributed_source_signal = gpuArray(distributed_source_signal);
    end

    % Loop through elements for both transducers
    for transducer_idx = 1:2
        for el_i = 1:transducer_pars_all{transducer_idx}.n_elements
            idx = el_i + (transducer_idx - 1) * transducer_pars_all{1}.n_elements;
            source_weights = squeeze(grid_weights_4d(idx,:,:,:));
            el_binary_mask = source_weights ~= 0;
            if canUseGPU
                source_weights = gpuArray(source_weights);
            end
            element_mask_ind = find(el_binary_mask);
            local_ind = ismember(mask_ind, element_mask_ind);
            distributed_source_signal(local_ind, :) = distributed_source_signal(local_ind, :) + source_weights(element_mask_ind) * cw_signals_all{transducer_idx}(el_i, :);
            source_labels = source_labels + idx * el_binary_mask;
        end
    end

    % Assign binary mask and source signals
    source.p_mask = binary_mask;
    source.p = distributed_source_signal;
end % TODO non-kWaveArray option should actually be another case distinction (see original function "setup_source")
end
