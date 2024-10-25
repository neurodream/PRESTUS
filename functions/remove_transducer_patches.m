function remove_transducer_patches(patchname)

    % patch names: transducerPatch, simPatch

    if nargin < 1
        patchname = 'transducerPatch';
    end

    % to be used in command window when Figure from script
    % b2_transducer_plot is opened
    
    % Find and delete all patches with the tag 'transducerPatch'
    patches = findobj(gca, 'Tag', patchname);
    delete(patches);
end