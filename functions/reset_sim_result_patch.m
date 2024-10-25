function reset_sim_result_patch(varargin)

    % not really needed. for code maintenance, easier to remove via
    % remove_transducer_patches and then set again

    % TODO cannot handle multiple patches yet

    global data data_folder;

    remove_transducer_patches('simPatch');

    p = inputParser;

    addParameter(p, 'AddTransducer', true, @islogical);
    addParameter(p, 'Color', 'cyan', @(x) ischar(x) || (isnumeric(x) && numel(x) == 3)); % TODO derive color from existing patch
    addParameter(p, 'CutoffPerc', 0.999, @isnumeric);

    parse(p, varargin{:});

    add_transducer  = p.Results.AddTransducer;
    color           = p.Results.Color;
    cutoff_perc     = p.Results.CutoffPerc;

    % Flatten the 3D matrix & Sort the flattened values
    flattened_data = data(:);
    sorted_data = sort(flattened_data);
    
    cutoff = sorted_data(floor(cutoff_perc*numel(sorted_data)));


    Ds_sim = smooth3(double(data > cutoff));
    sim_isosurface = isosurface(Ds_sim,0.5);
    sim_patch = patch(sim_isosurface,'FaceColor',color,'EdgeColor','none','facealpha',0.5);
    sim_patch.Tag = 'simPatch';
    
    if add_transducer
        % obtain target and transducer pos from folder label
        targetStr   = regexp(data_folder, 'targ_(\d+)-(\d+)-(\d+)', 'tokens');
        posStr      = regexp(data_folder, 'tpos_(\d+)-(\d+)-(\d+)', 'tokens');
        target = str2double(targetStr{1});
        pos = str2double(posStr{1});

        add_transducer_shape(...
            pos([2 1 3]), target([2 1 3]), ...
            'Shrink', false, ...
            'NoCasing', true, ...
            'TargetType', 'line', ...
            'PatchName', 'simPatch' ...
            );
    end

end