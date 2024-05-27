function kwave_medium = setup_medium_pCT(parameters, medium_mask, pseudoCT_cropped)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                         Setup k-wave medium                       %
    %                                                                   %
    % This function sets up the medium for k-wave simulations to take   %
    % place in. The way this mask is set up is by starting with an      %
    % empty grid and adjust the acoustic and thermal parameters in the  %
    % shape of the provided mask. This can be repeated to incorporate   %
    % different kinds of tissue (see 'layered').                        %
    %                                                                   %
    % Some notes:                                                       %
    % - The alpha value can be tweaked at the end of this script.       %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % DECIDE IF I WANT TO START FROM GRID OF 0s OR GRID OF 1s FOR DENSITY
    % REMOVE THE IMAGES TO CHECK

    % Loads the medium settings from the config file
    medium = parameters.medium;

    % Finds maximum and minumum values
    HU_min = min(pseudoCT_cropped,[],'all');
    HU_max = max(pseudoCT_cropped,[],'all');

    % Creates an empty grid according to the grid dimensions
    grid_of_ones = ones(parameters.grid_dims);
    
    % Fills this grid with a baseline medium
    baseline_medium = medium.('water');
    
    % 'sound_speed' and 'density' are 3D arrays determining the speed of sound and the medium density within the simulation grid
    sound_speed = baseline_medium.sound_speed*grid_of_ones; 
    density = baseline_medium.density*grid_of_ones;    %waterDensity(temp_0) * ones(Nx,Ny,Nz);                     % [kg/m^3]
    %density = baseline_medium.density*zeros(parameters.grid_dims);

    % Assume a homogeneous attenuation across the skull of 0.6 dB/MHz/cm in water and 7.4 dB/MHz/cm in the skull
    alpha_0_true =  baseline_medium.alpha_0_true*grid_of_ones; %0.6 * ones(Nx,Ny,Nz);                              % [dB/MHz/cm] (Kremkau et al., 1981)
    alpha_power_true = baseline_medium.alpha_power_true*grid_of_ones;

    % 'thermal conductivity' and 'specific_heat' are constants that define
    % the transmission of heat and the amount of heat required to change
    % the temperature within a givven mass of the tissue.
    thermal_conductivity = baseline_medium.thermal_conductivity * grid_of_ones;                                    % [W/(m.K)]
    specific_heat = baseline_medium.specific_heat_capacity * grid_of_ones;   % [J/(kg.K)]

    % Changes the values of the acoustic and thermal properties in the
    % baseline_medium in the shape of the labelled mask
    
    labels = fieldnames(parameters.layer_labels);
    % Loops through each labelled layer to create a new mask
        for label_i = 1:length(labels)
            label_name = labels{label_i};
            if strcmp(label_name, 'water')
               continue % Loops to next label since the baseline medium is set to water anyways
            end
            if strcmp(label_name, 'skull_cortical') || strcmp(label_name, 'skull_trabecular')
                thermal_conductivity(medium_mask==label_i) = medium.(label_name).thermal_conductivity;                
                specific_heat(medium_mask==label_i) = medium.(label_name).specific_heat_capacity;

                % Offset of 1000 for CT values to use housfield2density
                % function
                pseudoCT_cropped(medium_mask==label_i) = pseudoCT_cropped(medium_mask==label_i) + 1000;
                % Ensure there are no negative values
                pseudoCT_cropped(medium_mask==label_i) = max(pseudoCT_cropped(medium_mask==label_i),1);
                % Apply the function
                density(medium_mask==label_i) = hounsfield2density(pseudoCT_cropped(medium_mask==label_i));
                % Ensure there are no negative values
                density = max(density,1);
                
                sound_speed(medium_mask==label_i) = 1.33.*density(medium_mask==label_i) + 167; 
                alpha_0_true(medium_mask==label_i) =  4+4.7.*[1-(pseudoCT_cropped(medium_mask==label_i)-1000 -HU_min)./(HU_max-HU_min)].^0.5;
                % added -1000 to compensate for the previous +1000
                %4.7 = 8.7 - 4 = alpha_bone_max-alpha_bone_min
                alpha_power_true(medium_mask==label_i) = medium.(label_name).alpha_power_true;
            else
            % Sets the parameters in the shape of the mask
            thermal_conductivity(medium_mask==label_i) = medium.(label_name).thermal_conductivity;                 % [W/(m.K)]
            specific_heat(medium_mask==label_i) = medium.(label_name).specific_heat_capacity;                      % [J/(kg.K)]

            sound_speed(medium_mask==label_i) = medium.(label_name).sound_speed; 
            density(medium_mask==label_i) = medium.(label_name).density;
            alpha_0_true(medium_mask==label_i) =  medium.(label_name).alpha_power_true; 
            alpha_power_true(medium_mask==label_i) = medium.(label_name).alpha_power_true;
            end
        end
    filename_density = fullfile(parameters.output_dir, sprintf('density'));
    niftiwrite(density, filename_density, 'Compressed',true);
    filename_sound_speed = fullfile(parameters.output_dir, sprintf('sound_speed'));
    niftiwrite(sound_speed, filename_sound_speed, 'Compressed',true);
    filename_alpha_0_true = fullfile(parameters.output_dir, sprintf('alpha_0_true'));
    niftiwrite(alpha_0_true, filename_alpha_0_true, 'Compressed',true);
    filename_alpha_power_true = fullfile(parameters.output_dir, sprintf('alpha_power_true'));
    niftiwrite(alpha_power_true, filename_alpha_power_true, 'Compressed',true);

    % Account for actual absorption behaviour in k-Wave, which varies when high
    % absorption is used (see https://doi.org/10.1121/1.4894790).

    % 'alpha_coeff' is tweaked so any alpha_power could be used, 
    % as long as the alpha_0_true and alpa_power_true are correct for a
    % given frequency.
    % Here, I use 2.
    alpha_power_fixed = 2;
    
    alpha_coeff = fitPowerLawParamsMulti(alpha_0_true, alpha_power_true, sound_speed, parameters.transducer.source_freq_hz, alpha_power_fixed );

    filename_alpha_0_true = fullfile(parameters.output_dir, sprintf('alpha_0_true_fit'));
    niftiwrite(alpha_0_true, filename_alpha_0_true, 'Compressed',true);
    filename_alpha_power_true = fullfile(parameters.output_dir, sprintf('alpha_power_true_fit'));
    niftiwrite(alpha_power_true, filename_alpha_power_true, 'Compressed',true);

    % Outputs the medium as a structure
    kwave_medium = struct('sound_speed', sound_speed, ...
                          'density', density, ...
                          'alpha_coeff',alpha_coeff,...
                          'alpha_power', alpha_power_fixed , ...
                          'thermal_conductivity', thermal_conductivity,...
                          'specific_heat', specific_heat);

end