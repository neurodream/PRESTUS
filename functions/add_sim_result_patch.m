function data = add_sim_result_patch(data_folder, varargin)

    % TODO consider a clustering approach (DBSCAN?), and option to remove the near
    % field
    
    % TODO add option to add transducer
    
    % Create an input parser
    p = inputParser;
    
    % Add parameters with default values
    % addRequired(p, 'data_folder', @ischar);
    addParameter(p, 'Color', 'cyan', @(x) ischar(x) || (isnumeric(x) && numel(x) == 3));
    addParameter(p, 'CutoffPerc', 0.999, @isnumeric);
    addParameter(p, 'AddTransducer', false, @islogical);
    addParameter(p, 'PatchName', 'simPatch', @ischar);
    addParameter(p, 'HeadData', [], @isnumeric); % if provided, then limit to brain
    
    % Parse input arguments
    parse(p, varargin{:});
    
    % Retrieve values
    % data_folder     = p.Results.data_folder;
    color           = p.Results.Color;
    cutoff_perc     = p.Results.CutoffPerc;
    add_transducer  = p.Results.AddTransducer;
    patch_name      = p.Results.PatchName;
    head_data       = p.Results.HeadData;
    
    % fileInfo = dir(fullfile('D:\', data_folder, '*final_isppa_orig_coord.nii.gz'));
    fileInfo = dir(fullfile('M:\Documents\repos\PRESTUS_forked\data\sims', data_folder, '*final_isppa_orig_coord.nii.gz'));
    % fileInfo = dir(fullfile('M:\Documents\repos\PRESTUS_forked\data\sims\CTX500-024-010_77.0mm', data_folder, '*final_isppa_orig_coord.nii.gz')); 
    
    if numel(fileInfo) == 0
        disp('WARNING: could not find a nifti file with isppa values, generating them from pressure values now;')
        disp('- note that missing file likely means that simulation pipeline aborted early, likely meaning erroneous data -')
        disp('(expected load time of around 30 seconds')
        disp('loading data...')
        fileInfo = dir(fullfile('D:\', data_folder, '*_layered_results.mat'));
        load(fullfile(fileInfo.folder, fileInfo.name));
        disp('... loading data done.')
        p_max_CPU = gather(sensor_data.p_max_all);

        data = p_max_CPU.^2./(2*(kwave_medium.sound_speed.*kwave_medium.density)).*1e-4;

        % TODO seems like there is some transformation missing; check the
        % output var inv_final_transformation_matrix from running function
        % preprocess_brain
        % because the patch looks out of the brain (weird because it should
        % already exclude outside-brain voxels) - needs another sanity
        % check?
        sbj_id_tokens = regexp(data_folder, 'sbj_(\d+)', 'tokens');
        subject_id = str2num(sbj_id_tokens{1}{1});
        parameters.data_path = 'M:\Documents\repos\PRESTUS_forked\..\..\scans';
        [~, ~, ~, ~, ~, t1_image_orig, ~, ~, inv_final_transformation_matrix] = preprocess_brain(parameters, subject_id, false);

        data = tformarray(data, inv_final_transformation_matrix, ...
                                            makeresampler('cubic', 'fill'), [1 2 3], [1 2 3], size(t1_image_orig), [], 0) ;

        % delete unneeded vars from workspace again
        clear p_max_CPU source sensor sensor_data parameters kwave_medium kwave_input_args kgrid

    else
        data = niftiread(fullfile(fileInfo.folder, fileInfo.name));
    end

    % obtain target and transducer pos from folder label
    targetStr   = regexp(data_folder, 'targ_(\d+)-(\d+)-(\d+)', 'tokens');
    target = str2double(targetStr{1});
    posStr      = regexp(data_folder, 'tpos_(\d+)-(\d+)-(\d+)', 'tokens');
    pos = str2double(posStr{1});
    
    if numel(head_data) > 0
        % get the brain (assuming "head" exists globally)
        within_brain = ismember(head_data, [1 2 3]); % brain

        % shrink the shape: conservative brain estimate, so to prevent
        % "near field" in estimation of max pos
        se = strel('sphere', 1); % A spherical structuring element with a radius of 1
        within_brain = imerode(within_brain, se);

        % within = ismember(head, [1 2 3 4 7 8 9]); % skull
        within_brain = ~within_brain;
        within_brain = imfill(within_brain, 'holes');
        within_brain = ~within_brain;
    
        % limit to brain
        data(~within_brain) = 0;
    end

    % plot position of max:
    % Find the linear index of the maximum value
    [max_val, linear_index] = max(data(:));
    % Convert the linear index to [x, y, z] coordinates
    % [x_max, y_max, z_max] = ind2sub(size(data), linear_index);
    [x_max, y_max, z_max] = ind2sub(size(data), linear_index);
    max_pos = [x_max, y_max, z_max];
    disp(max_pos)
    % add cross
    add_3d_cross(max_pos([2 1 3]), 15);

    % Flatten the 3D matrix
    flattened_data = data(:);
    
    % Sort the flattened values
    sorted_data = sort(flattened_data);
    
    cutoff = sorted_data(floor(cutoff_perc*numel(sorted_data)));
    
    Ds_sim = smooth3(double(data > cutoff));
    sim_isosurface = isosurface(Ds_sim,0.5);
    sim_patch = patch(sim_isosurface,'FaceColor',color,'EdgeColor','none','facealpha',0.3);
    sim_patch.Tag = patch_name;
    
    % print diagnostic info to console

    % pressureStr = regexp(data_folder, 'amp_(\d+)', 'tokens');
    % pressureStr = pressureStr{1}{1};
    % pressure = str2double(pressureStr);
    % disp('free water intensity: ');
    % disp(['  ' num2str(get_intensity_from_source_amp(pressure)) ' W/cm^2']);

    disp('max Isppa: ');
    disp(['  ' num2str(data(x_max, y_max, z_max)) ' W/cm^2']);

    disp('Focus depth: (TODO check unit - likely needs scaling by voxel size)')
    disp(['  ' num2str(norm(pos - target)) ' mm (?)']);

    if add_transducer

        add_transducer_shape(...
            pos([2 1 3]), target([2 1 3]), ...
            'Shrink', false, ...
            'NoCasing', true, ...
            'TargetType', 'line', ...
            'PatchName', patch_name);
    end

end