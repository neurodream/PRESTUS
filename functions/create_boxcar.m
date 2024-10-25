function y_interpolated = create_boxcar(center_x, width, x_given, min_y, max_y)

% (note: docstring created by ChatGPT)
% CREATE_BOXCAR Generates a boxcar (rectangular) function and interpolates it to given x values.
%
%   y_interpolated = CREATE_BOXCAR(center_x, width, x_given, min_y, max_y)
%
%   This function creates a boxcar (rectangular) function centered at 'center_x' with a
%   width of 'width'. The function returns the interpolated values of the boxcar at the
%   specified x positions in 'x_given'. The y-values of the boxcar outside of the defined
%   width are set to 'min_y', while the y-values inside the boxcar are set to 'max_y'.
%
%   Inputs:
%       center_x - The center of the boxcar along the x-axis.
%       width - The width of the boxcar.
%       x_given - A vector of x-values where the boxcar function will be interpolated.
%       min_y (optional) - The y-value outside of the boxcar (default: 0).
%       max_y (optional) - The y-value inside the boxcar (default: 1).
%
%   Outputs:
%       y_interpolated - The interpolated y-values of the boxcar function at the specified
%                        x-values given in 'x_given'.
%
%   Example:
%       center_x = 5;
%       width = 2;
%       x_given = linspace(0, 10, 100);
%       y = create_boxcar(center_x, width, x_given);


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
    
    % Calculate the boundaries of the boxcar
    left_edge = center_x - width/2;
    right_edge = center_x + width/2;
    
    % Create x range for the boxcar function
    x = linspace(min_x, max_x, 1000);
    
    % Initialize the boxcar function
    y = min_y * ones(size(x));
    y(x >= left_edge & x <= right_edge) = max_y;
    
    % Interpolate the boxcar function to the given x values
    y_interpolated = interp1(x, y, x_given, 'linear');
end
