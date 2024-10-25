function add_3d_cross(center, length, color)
    % add_3d_cross(center, length, color)
    % center: [x, y, z] coordinates of the center of the cross
    % length: Length of each arm of the cross
    % color: Color of the cross (optional, default is 'k' (black))

    if nargin < 3
        color = 'k'; % default color is black
    end

    % Extract center coordinates
    x0 = center(1);
    y0 = center(2);
    z0 = center(3);

    % Calculate end points of the cross lines
    x = [x0 - length/2, x0 + length/2];
    y = [y0 - length/2, y0 + length/2];
    z = [z0 - length/2, z0 + length/2];

    % Plot the cross
    hold on;
    l1 = plot3([x0 x0], [y0 y0], [z(1) z(2)], color, 'LineWidth', 2); % z-axis line
    l2 = plot3([x(1) x(2)], [y0 y0], [z0 z0], color, 'LineWidth', 2); % x-axis line
    l3 = plot3([x0 x0], [y(1) y(2)], [z0 z0], color, 'LineWidth', 2); % y-axis line
    set(l1, 'Tag', 'cross'); % Tag for the first line
    set(l2, 'Tag', 'cross'); % Tag for the first line
    set(l3, 'Tag', 'cross'); % Tag for the first line
    hold off;
end