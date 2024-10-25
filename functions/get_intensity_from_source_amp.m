

function result = get_intensity_from_source_amp(source_amp)
    % TODO figure out why correction factor needed
    corr_factor = 59.4;
    % TODO read out from config rather than hardcoding
    result = ((source_amp^2)/(2*994*1500))/corr_factor;
end



