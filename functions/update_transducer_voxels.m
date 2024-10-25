function [transducer_voxels] = update_transducer_voxels(trans_pos, target, parameters, grid_dims)
    
    % make a grid of transducer
    
    add_transducer_shape(trans_pos, target, 'Parameters', parameters, 'TargetType', 'none');
    
    % get the faces and vertices for use in voxelization script
    bowl2case = findobj('Type', 'Patch', 'Tag', 'transducerPatch');
    bowl = findobj('Type', 'Surface', 'Tag', 'transducerPatch');
    bowl = surf2patch(bowl);
    
    FV = [];
    FV.faces = triangulate_faces(bowl2case.Faces);
    FV.vertices = bowl2case.Vertices;
    bowl2case_voxels = polygon2voxel(FV, grid_dims, 'none');
    
    FV.faces = triangulate_faces(bowl.faces);
    FV.vertices = bowl.vertices;
    bowl_voxels = polygon2voxel(FV, grid_dims, 'none');

    delete(findobj('Tag', 'transducerPatch')); % to only show voxel representation, not the patches
    
    transducer_voxels = bowl_voxels | bowl2case_voxels;

end