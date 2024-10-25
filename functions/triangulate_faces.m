function [faces_tri] = triangulate_faces(faces)

    faces_tri = zeros(2*size(faces, 1), 3);
    
    for i = 1:size(faces,1)
        v1 = faces(i, 1);
        v2 = faces(i, 2);
        v3 = faces(i, 3);
        v4 = faces(i, 4);
        
        % Define the two triangles for this quadrilateral
        faces_tri(2*i-1, :) = [v1, v2, v3];  % First triangle
        faces_tri(2*i, :) = [v1, v3, v4];    % Second triangle
    end

end