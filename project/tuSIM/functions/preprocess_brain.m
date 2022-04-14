function [skull_mask, segmented_image_cropped, skull_edge, trans_pos_final, focus_pos_final] = preprocess_brain(parameters, subject_id)
    %% CHECK INPUTS AND TRANSLATE PATTERNS
    disp('Checking inputs...')
    
    filename_t1 = fullfile(parameters.data_path, sprintf(parameters.t1_path_template, subject_id));
    filename_t2 = fullfile(parameters.data_path, sprintf(parameters.t2_path_template, subject_id));
    
    files_to_check = ["filename_t1", "filename_t2"];
    if parameters.transducer_from_localite
        localite_file = fullfile(parameters.data_path, sprintf(parameters.localite_instr_file_template, subject_id));
        files_to_check = [files_to_check, "localite_file"];
    end
    % check if files exist
    for filename_var = files_to_check
        eval(sprintf('filename = %s;',  filename_var ))
        % if there is a wildcard in the string, use dir to find file
        if contains(filename, '*')
            matching_files = dir(filename);
            if length(matching_files)>1
                error('More than 1 file matches the template %s', filename)
            elseif isempty(matching_files)
                error('No files match the template %s', filename)
            else 
                filename = fullfile(matching_files.folder , matching_files.name);
                eval(sprintf('%s = filename;',  filename_var ))
            end
        end

        if ~isfile(filename)
            error(sprintf('File does not exist: \r\n%s', filename));
        end
    end
    %% LOAD AND PROCESS t1 IMAGE (OR CT IMAGE OR CONTRAST IMAGE)
    disp('Loading T1...')

    % load mri image
    t1_image = niftiread(filename_t1);
    t1_header= niftiinfo(filename_t1);

    if parameters.transducer_from_localite
        % position the transducer on original (unrotated) T1 
        [trans_pos_grid, focus_pos_grid, ~, ~] = ...
            position_transducer_localite(localite_file, t1_header, parameters);
    else
        trans_pos_grid = parameters.transducer.pos_t1_grid';
        focus_pos_grid = parameters.focus_pos_t1_grid';
    end

    [t1_with_trans_img, transducer_pars] = plot_t1_with_transducer(t1_image, t1_header.PixelDimensions(1), trans_pos_grid, focus_pos_grid, parameters);

    %% SEGMENTATION
    disp('Starting segmentation...')

    headreco_folder = fullfile(parameters.data_path, sprintf('m2m_sub-%03d', subject_id));
    filename_segmented_headreco = fullfile(headreco_folder, sprintf('sub-%03d_final_contr.nii.gz', subject_id));

    if confirm_overwriting(filename_segmented_headreco, parameters) && (~isfield( parameters,'overwrite_simnibs') || parameters.overwrite_simnibs || ~exist(filename_segmented_headreco, 'file'))
        % double-check since segmentation takes a long time
        if parameters.interactive == 0 || confirmation_dlg('This will run SEGMENTATION WITH SIMNIBS that takes a long time, are you sure?', 'Yes', 'No')
            run_headreco(parameters.data_path, subject_id, filename_t1, filename_t2, parameters.simnibs_env_path);
            disp('\nThe script will continue with other subjects in the meanwhile...')
            skull_mask = [];
            segmented_image_cropped = [];
            skull_edge = [];
            trans_pos_final = [];
            focus_pos_final = [];
            return;
        end
    end

    %% Rotate to match the stimulation trajectory
    disp('Rotating to match the focus axis...')

    assert(exist(filename_segmented_headreco,'file')>0, ...
        'Head segmentation is not completed (%s does not exist), see logs in the batch_logs folder and in %s folder',...
            filename_segmented_headreco, headreco_folder)

    filename_reoriented_scaled_data = fullfile(parameters.data_path, sprintf('sub-%03d_after_rotating_and_scaling%s.mat', subject_id, parameters.results_filename_affix));

    if confirm_overwriting(filename_reoriented_scaled_data, parameters)
        segmented_img_orig = niftiread(filename_segmented_headreco);
        segmented_hdr_orig = niftiinfo(filename_segmented_headreco);

        scale_factor = segmented_hdr_orig.PixelDimensions(1)/parameters.grid_step_mm;

        [segmented_img_rr, trans_pos_upsampled_grid, focus_pos_upsampled_grid, scale_rotate_recenter_matrix, rotation_matrix] = ...
            align_to_focus_axis_and_scale(segmented_img_orig, t1_header, trans_pos_grid, focus_pos_grid, scale_factor, parameters);

        t1_img_rr = align_to_focus_axis_and_scale(t1_image, t1_header, trans_pos_grid, focus_pos_grid, scale_factor, parameters);

        assert(isequal(size(trans_pos_upsampled_grid,1:2),size(focus_pos_upsampled_grid, 1:2)),...
            "After reorientation, the first two coordinates of the focus and the transducer should be the same")

        save(filename_reoriented_scaled_data, 'segmented_img_rr', 'trans_pos_upsampled_grid', 'focus_pos_upsampled_grid', 'scale_rotate_recenter_matrix', 'rotation_matrix', 't1_img_rr');
    else 
        load(filename_reoriented_scaled_data);
    end

    %% Plot the skin & skull from simnibs

    % unsmoothed skull & skin masks
    skull_mask_unsmoothed = segmented_img_rr==4;
    skin_mask_unsmoothed = segmented_img_rr==5;

    skin_slice = squeeze(skin_mask_unsmoothed(:,trans_pos_upsampled_grid(2),:));
    skull_slice = squeeze(skull_mask_unsmoothed(:,trans_pos_upsampled_grid(2),:));

    t1_slice =  repmat(mat2gray(squeeze(t1_img_rr(:,trans_pos_upsampled_grid(2),:))), [1 1 3]);
    % skull is blue, skin is green
    skin_skull_img = cat(3, zeros(size(skull_slice)), skin_slice, skull_slice);

    % plot images for comparison
    montage(cat(4, t1_slice*255, skin_skull_img*255, imfuse(mat2gray(t1_slice), skin_skull_img, 'blend')) ,'size',[1 NaN]);

    %% SMOOTH & CROP SKULL
    disp('Smoothing and cropping the skull...')

    filename_cropped_smoothed_skull_data = fullfile(parameters.data_path, sprintf('sub-%03d_after_cropping_and_smoothing%s.mat', subject_id, parameters.results_filename_affix));
    if confirm_overwriting(filename_cropped_smoothed_skull_data, parameters)

        [skull_mask, skull_edge, segmented_image_cropped, trans_pos_final, focus_pos_final, ~, ~, new_grid_size, crop_translation_matrix] = ...
            smooth_and_crop_skull(segmented_img_rr, parameters.grid_step_mm, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters);
        save(filename_cropped_smoothed_skull_data, 'skull_mask', 'skull_edge', 'segmented_image_cropped', 'trans_pos_final', 'focus_pos_final', 'new_grid_size', 'crop_translation_matrix')
    else 
        load(filename_cropped_smoothed_skull_data);
    end    
    parameters.grid_dims = new_grid_size;


    % check that the transformations can be correctly reversed

    final_transformation_matrix = scale_rotate_recenter_matrix*crop_translation_matrix';
    inv_final_transformation_matrix = maketform('affine', inv(final_transformation_matrix')');

    backtransf_coordinates = round(tformfwd([trans_pos_final, focus_pos_final]', inv_final_transformation_matrix));
    if ~all(all(backtransf_coordinates ==[trans_pos_grid focus_pos_grid]'))
        disp('Backtransformed focus and transducer parameters differ from the original ones. Something went wrong (but note that small rounding errors could be possible.')
        disp('Original coordinates')
        disp([trans_pos_final, focus_pos_final]')
        disp('Backtransformed coordinates')
        disp(backtransf_coordinates)
        exit()
    end

    
end