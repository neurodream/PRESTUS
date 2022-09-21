function [res_image, transducer_pars] = plot_t1_with_transducer(t1_image, voxel_size_mm, trans_pos_grid, focus_pos_grid, parameters, plot_options)
    arguments
        t1_image (:,:,:) double
        voxel_size_mm (1,1) double
        trans_pos_grid (1,3) double
        focus_pos_grid (1,3) double
        parameters struct
        plot_options.slice_dim (1,1) double = 2
        plot_options.slice_ind (1,1) double = 0

    end
    if isempty(plot_options.slice_ind) || plot_options.slice_ind == 0
        plot_options.slice_ind = trans_pos_grid(plot_options.slice_dim);
    end
    
    im_size = max(size(t1_image), trans_pos_grid);
    if ~isequal(im_size, size(t1_image))
        t1_image = padarray(t1_image, im_size-size(t1_image),'post');
    end
    [transducer_bowl, ~, transducer_pars] = transducer_setup(parameters.transducer, trans_pos_grid, focus_pos_grid, ...
                                                            im_size, voxel_size_mm);

    % define cell array:
    slice_pointer = repmat({':'},1, 3);
    slice_pointer{plot_options.slice_dim} = plot_options.slice_ind;
    % comma-separated list to supply indices:

    t1_slice = double(t1_image(slice_pointer{:}));
    res_image = repmat(mat2gray(squeeze(t1_slice)),[1,1,3]);
    
    transducer_bowl = squeeze(transducer_bowl(slice_pointer{:}));
    res_image(:,:,1) = max(res_image(:,:,1), transducer_bowl); % red transducer
    
%     focus_pos_ind = 
    focus_pos_ind = find(~cellfun(@(x) ~strcmp(x,':'),slice_pointer));
    res_image(focus_pos_grid(focus_pos_ind(1))+(-2:2),focus_pos_grid(focus_pos_ind(2))+(-2:2),2) = 1; % blue focal point

end
