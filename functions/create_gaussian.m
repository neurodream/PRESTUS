function y_interpolated = create_gaussian(center_x, width, x_given, min_y, max_y)
% CREATE_GAUSSIAN Generates a Gaussian function and interpolates it to given x values.
%
%   y_interpolated = CREATE_GAUSSIAN(center_x, width, x_given, min_y, max_y)
%
%   This function creates a Gaussian function centered at 'center_x' with a
%   standard deviation derived from 'width' for scaling. The function returns the
%   interpolated values of the Gaussian at the specified x positions in 'x_given'.
%   The y-values are scaled between 'min_y' and 'max_y'.
%
%   Inputs:
%       center_x - The center of the Gaussian along the x-axis.
%       width - A parameter controlling the width of the Gaussian.
%       x_given - A vector of x-values where the Gaussian function will be interpolated.
%       min_y (optional) - The minimum y-value for scaling (default: 0).
%       max_y (optional) - The maximum y-value for scaling (default: 1).
%
%   Outputs:
%       y_interpolated - The interpolated y-values of the Gaussian function at the specified
%                        x-values given in 'x_given'.
%
%   Example:
%       center_x = 5;
%       width = 2;
%       x_given = linspace(0, 10, 100);
%       y = create_gaussian(center_x, width, x_given);

    % Set default values for min_y and max_y if not provided
    if nargin < 4
        min_y = 0;
    end
    if nargin < 5
        max_y = 1;
    end
    
    % Derive min_x and max_x from x_given
    min_x = min(x_given);
    max_x = max(x_given);
    
    % Create x range for the Gaussian function
    x = linspace(min_x, max_x, 1000);
    
    % Calculate the Gaussian standard deviation from the width
    sigma = width / 2;  % Adjust width to control spread of the Gaussian
    
    % Define the Gaussian function
    y = exp(-((x - center_x).^2) / (2 * sigma^2));
    
    % Scale y to fit within [min_y, max_y]
    y = y / max(y) * (max_y - min_y) + min_y;
    
    % Interpolate the Gaussian function to the given x values
    y_interpolated = interp1(x, y, x_given, 'linear');
end
