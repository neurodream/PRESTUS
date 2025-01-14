function [measure_data] = postprocessing_quantification(sub_id, parameters, medium_masks, containers, output_pressure_file, data)

ROI_r = round(parameters.focus_area_radius/parameters.grid_step_mm);

% create output containers
measure_data = []; % measures = {'pressure', 'intensity', 'mechanicalindex', 'heating', 'maxCEM43'};

% % targets
T = readtable('data/transducer_pos/position_LUT.xlsx');
% targetL = [T.y_l(T.sbj_ID == sub_id) T.x_l(T.sbj_ID == sub_id) T.z_l(T.sbj_ID == sub_id)];
targetR = [T.x_r(T.sbj_ID == sub_id) T.y_r(T.sbj_ID == sub_id) T.z_r(T.sbj_ID == sub_id)];

% Equation of the sphere: (x - px)^2 + (y - py)^2 + (z - pz)^2 <= r^2
[x, y, z] = ndgrid(1:size(medium_masks,1), 1:size(medium_masks,2), 1:size(medium_masks,3));
% ROItarget_L = (x - targetL(1)).^2 + (y - targetL(2)).^2 + (z - targetL(3)).^2 <= r^2;
ROItarget_R = (x - targetR(1)).^2 + (y - targetR(2)).^2 + (z - targetR(3)).^2 <= ROI_r^2;

% segmenting tissues

% Creates a logical skull mask and register skull_ids
labels = fieldnames(parameters.layer_labels);
skull_i = find(strcmp(labels, 'skull_cortical'));
trabecular_i = find(strcmp(labels, 'skull_trabecular'));
all_skull_ids = [skull_i, trabecular_i];
value_masks.skull = ismember(medium_masks,all_skull_ids);
brain_i = find(strcmp(labels, 'brain'));
value_masks.brain = ismember(medium_masks,brain_i);
skin_i = find(strcmp(labels, 'skin'));
value_masks.scalp = ismember(medium_masks,skin_i);
value_masks.whole = true(size(medium_masks));

% calculations for inner brain
SE = strel('sphere', 6);
value_masks.skulldilated = imdilate(value_masks.skull, SE);
% braineroded  = imerode(brain, SE); % not optimal: erodes around the gyri
value_masks.braininner = value_masks.brain & ~value_masks.skulldilated;
% manually defined targets
value_masks.target = ROItarget_R;
value_masks.offtarget = value_masks.brain & ~value_masks.target;

DVs = fieldnames(data);
for i = 1:numel(fieldnames(data))
    fieldname = DVs{i};
    if contains(fieldname, 'pressure')
        data.(fieldname)    = data.(fieldname)/1000000; % convert to MPa
    end
end

%% fill output containers

for outmeasure_ID = 1:length(containers)
    outmeasure = containers{outmeasure_ID};
    outmeasure_parts = cellstr(split(outmeasure, '_'));

    if strcmp(outmeasure_parts{1}, 'FWHM')

        % calculate FWHM overlap
        FWHM = get_FWHM(data.pressure, 0.25, false, false); % -6 dB, hence 0.25
        overlap = FWHM & value_masks.target;
        measure_data.(outmeasure) = (nnz(overlap) / nnz(FWHM))*100;

    else
        
        if strcmp(outmeasure_parts{3}, 'FW')
            values = data.([outmeasure_parts{1} '_FW']);
            mask_name = 'whole';
        else
            values = data.(outmeasure_parts{1});
            mask_name = outmeasure_parts{3};
        end
        
        if strcmp(outmeasure_parts{2}, 'max')
            output_function = @max;
        elseif strcmp(outmeasure_parts{2}, 'avg')
            output_function = @mean;
        elseif strcmp(outmeasure_parts{2}, '95')
            output_function = @(x) prctile(x, 95);
        else
            continue
        end

        measure_data.(outmeasure) = output_function(values(value_masks.(mask_name)));

    end

end

%% Save values to file
out_table = rows2vars(struct2table(measure_data));

% Check if file exists
if isfile(output_pressure_file)
    % Read existing data
    existingTable = readtable(output_pressure_file);
    
    % Append new data
    combinedTable = [existingTable; out_table];
else
    % If file doesn't exist, just use the new table
    combinedTable = out_table;
end

% Write combined data to the file
writetable(combinedTable, output_pressure_file);

end