function [head] = fill_head(head)

% input: head: 3D logical matrix

% try to fill the holes
% (not perfect approach here, because it is blind to parallel concavities
% and could erroneously fill some outside surface concavities)
for i = 1:size(head, 1)
    slice = head(i, :, :); % Extract the i-th slice (2D)
    slice = squeeze(slice);  % Remove singleton dimensions (make it 2D)
    filledSlice = imfill(slice, 'holes'); % Fill the holes in the 2D slice
    head(i, :, :) = filledSlice; % Replace the original slice with the filled slice
end
for i = 1:size(head, 2)
    slice = head(:, i, :); % Extract the i-th slice (2D)
    slice = squeeze(slice);  % Remove singleton dimensions (make it 2D)
    filledSlice = imfill(slice, 'holes'); % Fill the holes in the 2D slice
    head(:, i, :) = filledSlice; % Replace the original slice with the filled slice
end
for i = 1:size(head, 3)
    slice = head(:, :, i); % Extract the i-th slice (2D)
    slice = squeeze(slice);  % Remove singleton dimensions (make it 2D)
    filledSlice = imfill(slice, 'holes'); % Fill the holes in the 2D slice
    head(:, :, i) = filledSlice; % Replace the original slice with the filled slice
end

end