function [data_isppa, data_p, data_mi] = acoustic_mat_to_nifti(sbj_id, identifier, filepath, store_nifti)

fname_in = sprintf('sub-%03d_%s_results%s.mat', sbj_id, 'layered', identifier);
fname_out_isppa = sprintf('sub-%03d_%s_isppa%s', sbj_id, 'layered', identifier);
fname_out_p = sprintf('sub-%03d_%s_pressure%s', sbj_id, 'layered', identifier);

disp('loading file, might take a few moments...');
load(fullfile(filepath, fname_in), 'parameters', 'sensor_data', 'kwave_medium', 't1_header', 'inv_final_transformation_matrix');
disp('loading file done.');

p = gather(sensor_data.p_max_all);
Isppa_map = p.^2 ./ (2 * (kwave_medium.sound_speed .* kwave_medium.density)) * 1e-4;
MI_map = (p/10^6)/sqrt((parameters.transducers(1).source_freq_hz/10^6)); % TODO assumes that source frequency is the same for all transducers

% backtransform the data
data_isppa = tformarray(Isppa_map, inv_final_transformation_matrix, makeresampler('cubic', 'fill'), [1 2 3], [1 2 3], t1_header.ImageSize, [], 0) ;
data_p     = tformarray(p,         inv_final_transformation_matrix, makeresampler('cubic', 'fill'), [1 2 3], [1 2 3], t1_header.ImageSize, [], 0) ;
data_mi    = tformarray(MI_map,    inv_final_transformation_matrix, makeresampler('cubic', 'fill'), [1 2 3], [1 2 3], t1_header.ImageSize, [], 0) ;

if store_nifti
    t1_header.Datatype = 'single';
    niftiwrite(data_isppa, fullfile(filepath, fname_out_isppa), t1_header, 'Compressed', true);
    niftiwrite(data_p,     fullfile(filepath, fname_out_p),     t1_header, 'Compressed', true);
end