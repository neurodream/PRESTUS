function [swappedData, swappedInfo] = swapNiftiXY(data, info)
    swappedData = permute(data, [2,1,3]);  % swap x<->y
    swappedInfo = info;
    T = swappedInfo.Transform.T;

    % Swap rows 1 and 2
    T([1,2], :) = T([2,1], :);

    % Swap columns 1 and 2
    T(:, [1,2]) = T(:, [2,1]);

    swappedInfo.Transform.T = T;
end