function data = add_sim_result_patch(parameters, sbj_ID, varargin)

    % TODO consider a clustering approach (DBSCAN?), and option to remove the near
    % field
    
    % TODO add option to add transducer
    
    % Create an input parser
    p = inputParser;
    
    % Add parameters with default values
    % addRequired(p, 'data_folder', @ischar);
    addParameter(p, 'Color', '#7E2F8E', @(x) ischar(x) || (isnumeric(x) && numel(x) == 3));
    addParameter(p, 'CutoffPerc', 0.999, @isnumeric);
    addParameter(p, 'PatchName', 'simPatch', @ischar);
    addParameter(p, 'HeadData', [], @isnumeric); % if provided, then limit to brain
    
    % Parse input arguments
    parse(p, varargin{:});
    
    % Retrieve values
    color           = p.Results.Color;
    cutoff_perc     = p.Results.CutoffPerc;
    patch_name      = p.Results.PatchName;
    head_data       = p.Results.HeadData;

    % sbj_ID = parameters.subject_subfolder;
    data_folder = fullfile(parameters.data_path, 'sim_outputs', sprintf('sub-%03d', sbj_ID));
    data_file = fullfile(data_folder, sprintf('sub-%03d_final_intensity%s.nii.gz', sbj_ID, parameters.results_filename_affix));
    data = niftiread(data_file);

    % fileInfo = dir(fullfile('M:\Documents\repos\PRESTUS_forked\data\sims', data_folder, '*final_isppa_orig_coord.nii.gz'));
    
    % data = niftiread(fullfile(fileInfo.folder, fileInfo.name));
    
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
    add_3d_cross(max_pos([2 1 3]), 50, 'black', 3);

    % 
    % Flatten the 3D matrix
    flattened_data = data(:);
    % Sort the flattened values
    sorted_data = sort(flattened_data);
    % based on this, compute cutoff
    cutoff = sorted_data(floor(cutoff_perc*numel(sorted_data)));
    
    thresholds = [0.1 0.2 0.4 0.8 1.6 3.2 6.4 12.8 1000] + 4;

    colors = {
        '#FFFF00'  % Yellow
        '#FFD700'  % Golden Yellow
        '#FFA500'  % Orange
        '#FF8C00'  % Dark Orange
        '#FF4500'  % Orange-Red
        '#FF0000'  % Red
        '#FF6347'  % Tomato
        '#FFFFFF'  % White
    };
    % colors = {
    %     '#4B0000'  % Very Dark Red
    %     '#800000'  % Dark Red
    %     '#B22222'  % Firebrick
    %     '#FF4500'  % Orange-Red
    %     '#FF6347'  % Tomato
    %     '#FF8C00'  % Dark Orange
    %     '#FFA500'  % Orange
    %     '#FFFFFF'  % White
    % };
    alphas = linspace(0.1, 0.1, 8);
    % for threshold_i = 1:numel(thresholds) - 1
    %     lower_threshold = thresholds(threshold_i);
    %     upper_threshold = thresholds(threshold_i + 1);
        Ds_sim = smooth3(double(data > cutoff)); %  & data < upper_threshold
        sim_isosurface = isosurface(Ds_sim,0.5);
        color_values = interp3(data, sim_isosurface.vertices(:,1), sim_isosurface.vertices(:,2), sim_isosurface.vertices(:,3));
        sim_patch = patch(sim_isosurface,'FaceColor','interp','FaceVertexCData', color_values, 'EdgeColor','none', 'facealpha',0.5);
        sim_patch.Tag = patch_name;
    % end
    colormap hot; % set desired colormap
    ax = gca; % get current axis
    ax.CLim = [0.5*cutoff, 0.75*max(data(:))];
    colorbar;
    
    % print diagnostic info to console

    % pressureStr = regexp(data_folder, 'amp_(\d+)', 'tokens');
    % pressureStr = pressureStr{1}{1};
    % pressure = str2double(pressureStr);
    % disp('free water intensity: ');
    % disp(['  ' num2str(get_intensity_from_source_amp(pressure)) ' W/cm^2']);


    disp('max Isppa: ');
    disp(['  ' num2str(data(x_max, y_max, z_max)) ' W/cm^2']);

    % TODO not possible with 2 transducers
    % disp('Focus depth: (TODO check unit - likely needs scaling by voxel size)')
    % disp(['  ' num2str(norm(pos - target)) ' mm (?)']);


end