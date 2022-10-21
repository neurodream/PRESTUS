function [smoothed_segmented_img, skull_edge, segmented_image_cropped, trans_pos_final, focus_pos_final, min_dims, max_dims, new_grid_dims, translation_matrix] = smooth_and_crop_layered(segmented_img, voxel_size_mm, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters)
    
    %% split into layers of interest
    
    % Smoothing window size
    windowSize = 4;
    % Creates an empty grid the size of the segmented image
    smoothed_segmented_img = zeros(size(segmented_img));
    
    % Adds a smoothing threshold to bone and other non-water tissue.
    labels = fieldnames(parameters.layer_labels);
    for label_i = 1:length(labels)
        if strcmp(labels{label_i}, 'water')
           continue
        end
        sim_nibs_layers = parameters.layer_labels.(labels{label_i});
        layer_mask = ismember(segmented_img, sim_nibs_layers);
        if strcmp(labels{label_i}, 'skull')
            smooth_threshold = parameters.skull_smooth_threshold;    
        else
            smooth_threshold = parameters.other_smooth_threshold;
        end
        
        layer_mask_smoothed = smooth_img(layer_mask, windowSize, smooth_threshold);
        smoothed_segmented_img(layer_mask_smoothed~=0) = label_i;
    end
    
    %% remove gaps between skull & skin
    if any(strcmp(labels,  'skull')) && any(strcmp(labels,  'skin'))
        skin_i = find(strcmp(labels,  'skin'));
        skull_i = find(strcmp(labels,  'skull'));
        skin = smoothed_segmented_img==skin_i;
        skull = smoothed_segmented_img==skull_i;
        skin_skull = skin+skull;
        skin_skull_filled = imfill(skin_skull);

        [labeledImage, ~] = bwlabeln(~skin_skull);
        blobMeasurements = regionprops(labeledImage, 'area');
        allAreas = [blobMeasurements.Area];

        [~, sortIndexes] = sort(allAreas, 'descend');
        biggestBlob = ismember(labeledImage, sortIndexes(2));

        smoothed_segmented_img((skin_skull_filled-skin_skull -  biggestBlob)>0) = skull_i;
    end
    
    % Shows the segmented figures before and after smoothing side-by-side
    imshowpair(label2rgb(squeeze(segmented_img(:,trans_pos_upsampled_grid(2),:))), label2rgb(squeeze(smoothed_segmented_img(:,trans_pos_upsampled_grid(2),:))), 'montage')
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_segmented_img_before_after_smoothing%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    title('Original (left) and smoothed (right) segmented image')
    export_fig(output_plot, '-native')

    csf_mask = segmented_img==3;
        
    % Creates a transducer for overlay in a figure
    transducer_bowl = transducer_setup(parameters.transducer, trans_pos_upsampled_grid, focus_pos_upsampled_grid, ...
                                                            size(segmented_img), voxel_size_mm);

    %% Cropping the bone tissue
    % Some bone tissue is not close to the neural tissue, and thus not essential for simulations.
    % The csf will be expanded to be used as a guide for what bone tissue should be left in.
    
    % To do this, the csf volume will first be expanded by a lot.
    SE = strel('cube', parameters.csf_mask_expansion_factor/voxel_size_mm);
    csf_mask_expanded = imdilate(csf_mask, SE);

%     imshow(plot_t1_with_transducer(segmented_image, segmented_hdr, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters))
    
    % Plots the result of the CSF expansion
    imshowpair(squeeze(smoothed_segmented_img(:,trans_pos_upsampled_grid(2),:)), squeeze(csf_mask_expanded(:,trans_pos_upsampled_grid(2),:)), 'falsecolor')
    title('CSF mask (pink) and the segmented image (white)')
    
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_skull_csf_mask%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    export_fig(output_plot, '-native')

    % The margin around the edge of the figure where the CSF will be cropped
    crop_margin = parameters.pml_size+1;
    
    % get_crop_dims ensures that the CSF mask is not expanded outside the
    % bounds of the image dimensions
    smoothed_segmented_img(~csf_mask_expanded) = 0;
    %imshow(squeeze(skull_mask_smoothed(trans_pos_upsampled_grid(1),:,:)))
    %montage({label2rgb(squeeze(segmented_img(trans_pos_upsampled_grid(1),:,:)),'jet','k'), label2rgb(squeeze(segmented_img(:,trans_pos_upsampled_grid(2),:)),'jet','k'),label2rgb(squeeze(segmented_img(:,:,3)),'jet','k')},'Size',[1 3])
    [min_dims, max_dims, new_grid_dims] = get_crop_dims(smoothed_segmented_img+transducer_bowl, crop_margin);
    
    % In case the minimum dimensions are too small, 
    if any(min_dims < 1)
        pad_amount = abs(min(min_dims, [1 1 1]));
        
        segmented_img = padarray(segmented_img, pad_amount,0,'pre');
        smoothed_segmented_img = padarray(smoothed_segmented_img, pad_amount,0,'pre');
        min_dims = max(min_dims, [1 1 1]);
        max_dims = max_dims+pad_amount;
        new_grid_dims = max_dims - min_dims + 1;
        
    end

    new_grid_dims(1) = find_min_factor(new_grid_dims(1),new_grid_dims(1)+parameters.prime_factor_max_grid_expansion);
    new_grid_dims(2) = find_min_factor(new_grid_dims(2),new_grid_dims(2)+parameters.prime_factor_max_grid_expansion);
    new_grid_dims(3) = find_min_factor(new_grid_dims(3),new_grid_dims(3)+parameters.prime_factor_max_grid_expansion);
    max_dims = min_dims+new_grid_dims - 1;
    
    % Ensures that the maximum size of the new grid matches the skull mask
    if any( max_dims > size(smoothed_segmented_img) )
        pad_amount = max( max_dims, size(segmented_img))-size(segmented_img);
        segmented_img = padarray(segmented_img, pad_amount,0,'post');
        smoothed_segmented_img = padarray(smoothed_segmented_img, pad_amount,0,'post');        
    end

    % Crops the figure
    smoothed_segmented_img = smoothed_segmented_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3)); % Crop image
    % Creates a mask around the edge of the figure to define the edge of the skull
    skull_edge = edge3(smoothed_segmented_img==find(strcmp(labels,'skull')), 'approxcanny',0.1);

    % Crops the segmented figure
    segmented_image_cropped = segmented_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3)); % Crop image

    % Saves the skull mask in the parameters structure as the new grid dimensions
    parameters.grid_dims = size(smoothed_segmented_img);

    % Moves the transducer and focus positions to the correct location in
    % the cropped grid
    trans_pos_final = trans_pos_upsampled_grid-min_dims';
    focus_pos_final = focus_pos_upsampled_grid-min_dims';
    
    % Creates a translation matrix for later reference
    translation_matrix = makehgtform('translate', -min_dims);
    
    disp('New grid dimensions:')
    disp(new_grid_dims)
    disp('New transducer position:')
    disp(trans_pos_final)
    
     % Shows and saves figures of the smoothed and unsmoothed skull mask with the transducer and focus locations
    imshowpair(plot_t1_with_transducer(smoothed_segmented_img, voxel_size_mm, trans_pos_final, focus_pos_final, parameters),...
        plot_t1_with_transducer(segmented_img, voxel_size_mm, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters),...
    'montage')
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_skull_final%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    title('Cropped, padded, & smoothed (left) and original (right) segmented mask')
    export_fig(output_plot, '-native')

end

% Currently not in use
function smoothed_img = smooth_img(unsmoothed_img, windowSize, threshold)
    arguments
       unsmoothed_img (:,:,:) double
       windowSize (1,1) double = 4
       threshold (1,1) double = 0.5       
    end
    kernel = ones(windowSize,windowSize,windowSize)/ windowSize ^ 3;
    blurryImage = convn(double(unsmoothed_img), kernel, 'same');
    smoothed_img = blurryImage > threshold; % Rethreshold
end
