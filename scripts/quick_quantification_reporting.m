close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(currentFile);
cd(rootpath); % repos/PRESTUS_forked/
cd ..

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

sham = false;
ROI_mm = 6;

% base config ("hard" params)
parameters = load_parameters('nico_test_double_acoustic_100mm_config.yaml');

% create output containers
measure_data = []; % measures = {'pressure', 'intensity', 'mechanicalindex', 'heating', 'maxCEM43'};
containers = {
    'pressure_max_FW', ...
    'pressure_max_brain', ...
    'pressure_max_scalp', ...
    'pressure_max_target', ...
    'pressure_max_offtarget', ...
    'pressure_avg_target', ...
    'pressure_avg_scalp', ...
    'pressure_avg_offtarget', ...
    'pressure_95_brain', ...
    'pressure_95_scalp', ...
    'intensity_max_FW', ...
    'mechanicalindex_max_FW', ...
    'mechanicalindex_max_whole', ...
    'mechanicalindex_max_brain', ...
    'mechanicalindex_max_scalp', ...
    'heating_max_brain', ...
    'heating_avg_target', ...
    'maxCEM43_max_whole', ...
    'maxCEM43_max_brain', ...
    'maxCEM43_max_scalp', ...
    'maxCEM43_95_whole', ...
    'maxCEM43_95_brain', ...
    'maxCEM43_95_scalp', ...
    'FWHM_in_ROI_perc', ...
};
containers = {
    'pressure_max_brain', ...
    'pressure_max_scalp', ...
    'pressure_max_target', ...
    'pressure_max_offtarget', ...
    'pressure_avg_target', ...
    'pressure_avg_scalp', ...
    'pressure_avg_offtarget', ...
    'pressure_95_brain', ...
    'pressure_95_scalp', ...
    'mechanicalindex_max_whole', ...
    'mechanicalindex_max_brain', ...
    'mechanicalindex_max_scalp', ...
    'heating_max_brain', ...
    'heating_avg_target', ...
    'maxCEM43_max_whole', ...
    'maxCEM43_max_brain', ...
    'maxCEM43_max_scalp', ...
    'maxCEM43_95_whole', ...
    'maxCEM43_95_brain', ...
    'maxCEM43_95_scalp', ...
    'FWHM_in_ROI_perc', ...
};
for container_ID = 1:length(containers)
    container = containers{container_ID};
    measure_data.(container) = zeros(1,6);
end

if sham
    targeting = '-';
else
    targeting = '--';
end

for sub_id = 8%1:6

    disp(['calculating subject ' num2str(sub_id) '...']);

    % load the MRI and segment
    nifti_path = fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', sub_id), 'final_tissues.nii.gz');
    
    layers = niftiread(nifti_path);
    info   = niftiinfo(nifti_path);
    
    % % targets
    T = readtable('data/transducer_pos/position_LUT.xlsx');
    % targetL = [T.y_l(T.sbj_ID == sub_id) T.x_l(T.sbj_ID == sub_id) T.z_l(T.sbj_ID == sub_id)];
    targetR = [T.x_r(T.sbj_ID == sub_id) T.y_r(T.sbj_ID == sub_id) T.z_r(T.sbj_ID == sub_id)];
    
    r = round(ROI_mm/mean(info.PixelDimensions)); % 6 mm around target
    % Equation of the sphere: (x - px)^2 + (y - py)^2 + (z - pz)^2 <= r^2
    [x, y, z] = ndgrid(1:size(layers,1), 1:size(layers,2), 1:size(layers,3));
    % ROItarget_L = (x - targetL(1)).^2 + (y - targetL(2)).^2 + (z - targetL(3)).^2 <= r^2;
    ROItarget_R = (x - targetR(1)).^2 + (y - targetR(2)).^2 + (z - targetR(3)).^2 <= r^2;
    
    % segmenting tissues
    
    value_masks = [];
    value_masks.whole = true(size(layers));
    value_masks.skull = layers == 7 | layers == 8;
    value_masks.brain = layers == 1 | layers == 2;
    value_masks.scalp = layers == 5;
    % calculations for inner brain
    SE = strel('sphere', 6);
    value_masks.skulldilated = imdilate(value_masks.skull, SE);
    % braineroded  = imerode(brain, SE); % not optimal: erodes around the gyri
    value_masks.braininner = value_masks.brain & ~value_masks.skulldilated;
    % manually defined targets
    value_masks.target = ROItarget_R;
    value_masks.offtarget = value_masks.brain & ~value_masks.target;
    
    % load the data
    data = [];
    measures = {'pressure', 'intensity', 'mechanicalindex', 'heating', 'maxCEM43'};
    
    for measure_ID = 1:length(measures)
    
        measure = measures{measure_ID};
        if strcmp(measure, 'heating') |  strcmp(measure, 'maxCEM43')
            basename = 'sub-%03d/sub-%03d_layered_%sL*%sr_R*%sr_it1_%s_imprecisionnone.nii.gz';
        else
            basename = 'sub-%03d/sub-%03d_final_%sL*%sr_R*%sr_it1_%s_imprecisionnone.nii.gz';
        end
    
        filename = sprintf(basename, sub_id, sub_id, measure, targeting, targeting, 'heatingtimeline');
        filename = dir(fullfile(parameters.data_path, 'sim_outputs', filename));
        % necessary because placeholder also finds '--' instead of only '-'
        if length(filename) > 1 & sham
            file_names = {filename.name};
            filename = filename(~contains(file_names, '--'));
        end
        data.(measure) = niftiread(fullfile(filename.folder, filename.name));
        
        try
            filename = sprintf(basename, sub_id, sub_id, measure, targeting, targeting, 'FW');
            filename = dir(fullfile(parameters.data_path, 'sim_outputs', filename));
            % necessary because placeholder also finds '--' instead of only '-'
            if length(filename) > 1 & sham
                file_names = {filename.name};
                filename = filename(~contains(file_names, '--'));
            end
            data.([measure '_FW']) = niftiread(fullfile(filename.folder, filename.name));
        catch
            % no FW for heating
        end
    
    end
    
    data.pressure    = data.pressure/1000000; % convert to MPa
    % data.pressure_FW = data.pressure_FW/1000000; % convert to MPa
    
    % fill output containers
    for outmeasure_ID = 1:length(containers)
        outmeasure = containers{outmeasure_ID};
        outmeasure_parts = cellstr(split(outmeasure, '_'));

        if strcmp(outmeasure_parts{1}, 'FWHM')

            % calculate FWHM overlap
            FWHM = get_FWHM(data.pressure, 0.25, false, false); % -6 dB, hence 0.25
            overlap = FWHM & value_masks.target;
            measure_data.(outmeasure)(1, sub_id) = (nnz(overlap) / nnz(FWHM))*100;

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

            measure_data.(outmeasure)(1, sub_id) = output_function(values(value_masks.(mask_name)));

        end

    end

end

disp_format = @(x) fprintf('%.2f < %.2f < %.2f\n', min(measure_data.(x)), median(measure_data.(x)), max(measure_data.(x)));
disp_format = @(x) fprintf('%s %.2f\n', x, max(measure_data.(x)));
disp_format_debug = @(x) fprintf('%f < %f < %f\n', min(measure_data.(x)), median(measure_data.(x)), max(measure_data.(x)));

% disp_format('pressure_max_FW');
disp_format('pressure_max_brain');
disp_format('pressure_max_target');
disp_format('pressure_max_offtarget');

fprintf('\n');

disp_format('pressure_avg_scalp');
disp_format('pressure_avg_target');
disp_format('pressure_avg_offtarget');

fprintf('\n');

% disp_format('intensity_max_FW');

fprintf('\n');

% disp_format('mechanicalindex_max_FW');
disp_format('mechanicalindex_max_whole');
disp_format('mechanicalindex_max_brain');
disp_format('mechanicalindex_max_scalp');

fprintf('\n');

disp_format('heating_max_brain');

fprintf('\n');

disp_format('heating_avg_target');

fprintf('\n');

disp_format('maxCEM43_max_whole');
disp_format('maxCEM43_max_brain');
disp_format('maxCEM43_max_scalp');

fprintf('\n');

disp_format('maxCEM43_95_whole');
disp_format('maxCEM43_95_brain');
disp_format('maxCEM43_95_scalp');

fprintf('\n');

disp_format('FWHM_in_ROI_perc');