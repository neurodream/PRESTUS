function [] = descriptive_stats_postprocessing()
    
    % big TODO
    % read measure_data from files
    % think if this is really a function or rather a script

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

end