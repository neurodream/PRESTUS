function [new_vertex, new_face] = move_on_sphere(faces, current_face, current_vertex, direction, step_size)
    % Default step_size to 1 if not provided
    if nargin < 5
        step_size = 1;
    end
    
    % Get the indices of the vertices for the current face
    face_vertices = faces(current_face, :);  % 4 vertices for this face
    
    % Initialize new face and vertex
    new_face = current_face;
    new_vertex = current_vertex;
    
    % Determine the new vertex based on the direction
    switch direction
        case 'd'
            % Move to the next vertex in the quad (wrap around), according to step_size
            for k = 1:step_size
                idx = find(face_vertices == new_vertex);
                new_vertex = face_vertices(mod(idx, 4) + 1);  % Move to next in the cycle
            end
            
        case 'a'
            % Move to the previous vertex in the quad (wrap around), according to step_size
            for k = 1:step_size
                idx = find(face_vertices == new_vertex);
                new_vertex = face_vertices(mod(idx - 2, 4) + 1);  % Move to previous in the cycle
            end
            
        case 'w'
            % Move up by step_size adjacent faces
            for k = 1:step_size
                temp_face = find_adjacent_face(faces, new_face, new_vertex, 'up');
                if ~isempty(temp_face)
                    new_face = temp_face;
                    new_vertex = find_new_vertex_on_face(faces, new_face, new_vertex);
                else
                    break;  % Stop if no more adjacent face found
                end
            end
            
        case 's'
            % Move down by step_size adjacent faces
            for k = 1:step_size
                temp_face = find_adjacent_face(faces, new_face, new_vertex, 'down');
                if ~isempty(temp_face)
                    new_face = temp_face;
                    new_vertex = find_new_vertex_on_face(faces, new_face, new_vertex);
                else
                    break;  % Stop if no more adjacent face found
                end
            end
            
        otherwise
            error('Invalid direction. Use "a" for left, "d" for right, "w" for up, or "s" for down.');
    end
end

function adjacent_face = find_adjacent_face(faces, current_face, current_vertex, direction)
    % Get the vertices of the current face
    face_vertices = faces(current_face, :);
    
    % Find the index of the current vertex in the face
    idx = find(face_vertices == current_vertex);
    
    % Determine the edge based on the direction
    switch direction
        case 'up'
            % Move along the top edge (vertices 1 and 2, or corresponding pair in the quad)
            edge_vertices = [face_vertices(mod(idx, 4) + 1), current_vertex];
        case 'down'
            % Move along the bottom edge (vertices 3 and 4, or corresponding pair in the quad)
            edge_vertices = [current_vertex, face_vertices(mod(idx - 2, 4) + 1)];
        case 'left'
            % Move along the left edge (between two vertices that form a vertical edge)
            edge_vertices = [current_vertex, face_vertices(mod(idx - 2, 4) + 1)];
        case 'right'
            % Move along the right edge (between two vertices that form a vertical edge)
            edge_vertices = [face_vertices(mod(idx, 4) + 1), current_vertex];
        otherwise
            error('Unsupported direction. Use "up", "down", "left", or "right".');
    end
    
    % Search for an adjacent face that shares the same edge
    adjacent_face = [];
    for f = 1:size(faces, 1)
        if f == current_face
            continue;  % Skip the current face
        end
        
        % Check if the face contains both vertices of the edge
        if all(ismember(edge_vertices, faces(f, :)))
            adjacent_face = f;
            return;  % Exit once the adjacent face is found
        end
    end
end

function new_vertex = find_new_vertex_on_face(faces, new_face, current_vertex)
    % Get the vertices of the new face
    new_face_vertices = faces(new_face, :);
    
    % Find the position of the current vertex in the new face
    idx = find(new_face_vertices == current_vertex);
    
    % If the current vertex is in the new face, find the next vertex
    if ~isempty(idx)
        new_vertex = new_face_vertices(mod(idx, 4) + 1);  % Move to next vertex in the cycle
    else
        % If the current vertex is not in the new face, return the first vertex
        new_vertex = new_face_vertices(1);
    end
end