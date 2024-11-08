function [error, i_axial_oneil] = phase_optimization_annulus(phase, parameters, velocity, axial_position, desired_focal_dist_mm)
  p_axial_oneil = focusedAnnulusONeil(parameters.transducer.curv_radius_mm/1e3, ...
    [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm]/1e3, repmat(velocity,1,parameters.transducer.n_elements), ...
    [0 phase], parameters.transducer.source_freq_hz, parameters.medium.water.sound_speed, ...
    parameters.medium.water.density, (axial_position-0.5)*1e-3);
  
  % [max_amp, max_pos] = max(p_axial_oneil);
  % halfMax = (min(p_axial_oneil) + max_amp) / 2;
  % % Find where the data first drops below half the max.
  % index1 = find(p_axial_oneil <= halfMax & axial_position<axial_position(max_pos), 1, 'last')+1;
  % % Find where the data last rises above half the max.
  % index2 = find(p_axial_oneil <= halfMax & axial_position>axial_position(max_pos), 1, 'first')-1;
  % 
  % actual_focal_dist_mm = (axial_position(index2) - axial_position(index1))/2+axial_position(index1);

  %actual_focal_dist_mm = axial_position(p_axial_oneil==max(p_axial_oneil));

  i_axial_oneil = p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4;

  % [peaks_list, peak_inds] = findpeaks(i_axial_oneil); 
  % actual_focal_dist_mm = axial_position(peak_inds(end));
  [~, max_ind] = max(i_axial_oneil);
  actual_focal_dist_mm = axial_position(max_ind);

  error = abs(desired_focal_dist_mm-actual_focal_dist_mm);

end