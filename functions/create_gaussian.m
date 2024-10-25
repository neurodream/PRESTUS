function y_interpolated = create_gaussian(center_x, width, x_given, min_y, max_y)
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
    
    % Calculate the Gaussian function
    sigma = width / (2 * sqrt(2 * log(2)));  % Convert width to standard deviation (FWHM to sigma)
    y = min_y + (max_y - min_y) * exp(-((x - center_x).^2) / (2 * sigma^2));
    
    % Interpolate the Gaussian function to the given x values
    y_interpolated = interp1(x, y, x_given, 'linear');
end
