function array = floodFill3D(array)
    % Ensure input is a 3D logical array
    assert(ndims(array) == 3, 'Input must be a 3D logical array.');
    array = logical(array);

    gap = 1; % this is just for try & error. function does not work yet.
    
    % Flood fill each boundary side using imfill
    % Left boundary
    leftBoundary = squeeze(array(1+gap, :, :));
    leftBoundary_filled = imfill(leftBoundary, 'holes');
    array(1+gap, :, :) = leftBoundary_filled;
    
    % Right boundary
    rightBoundary = squeeze(array(end-gap, :, :));
    rightBoundary_filled = imfill(rightBoundary, 'holes');
    array(end-gap, :, :) = rightBoundary_filled;
    
    % Front boundary
    frontBoundary = squeeze(array(:, 1+gap, :));
    frontBoundary_filled = imfill(frontBoundary, 'holes');
    array(:, 1+gap, :) = frontBoundary_filled;
    
    % Back boundary
    backBoundary = squeeze(array(:, end-gap, :));
    backBoundary_filled = imfill(backBoundary, 'holes');
    array(:, end-gap, :) = backBoundary_filled;
    
    % Top boundary
    topBoundary = squeeze(array(:, :, 1+gap));
    topBoundary_filled = imfill(topBoundary, 'holes');
    array(:, :, 1+gap) = topBoundary_filled;
    
    % Bottom boundary
    bottomBoundary = squeeze(array(:, :, end-gap));
    bottomBoundary_filled = imfill(bottomBoundary, 'holes');
    array(:, :, end-gap) = bottomBoundary_filled;


    % now fill the holes in 3D space
    array = ~array;
    array = imfill(array, 'holes');
    array = ~array;
end