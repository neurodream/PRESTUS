function reset_transducer_patches(varargin)

    global trans_pos target

    p = inputParser;

    addParameter(p, 'PosChange', [0, 0, 0], @(x) (isnumeric(x) && numel(x) == 3));
    addParameter(p, 'TargetChange', [0, 0, 0], @(x) (isnumeric(x) && numel(x) == 3));
    addParameter(p, 'TransDistChange', 0, @isnumeric); % TODO misnomer; other way around!
    addParameter(p, 'PosDistChange', 0, @isnumeric); % TODO misnomer; other way around!

    parse(p, varargin{:});

    % Default values
    posChange = p.Results.PosChange;
    targetChange = p.Results.TargetChange;
    transDistChange = p.Results.TransDistChange;
    posDistChange = p.Results.PosDistChange;

    % % Check for keyword arguments
    % for k = 1:2:length(varargin)
    %     switch varargin{k}
    %         case 'PosChange'
    %             posChange = varargin{k+1};
    %         case 'TargetChange'
    %             targetChange = varargin{k+1};
    %     end
    % end

    % TODO figure out inhowfar order matters here

    trans_pos_add = posChange([2 1 3])*2;
    trans_pos = trans_pos + trans_pos_add;
    target_add = targetChange([2 1 3])*2;
    target = target + target_add;

    direction = target - trans_pos;
    unit_direction = direction / norm(direction);
    new_target = target + transDistChange * unit_direction;
    new_pos = trans_pos + posDistChange * unit_direction;
    % careful here: rounding means destructive operations (TODO figure out
    % if necessary)
    target = [round(new_target(1)) round(new_target(2)) round(new_target(3))];
    trans_pos = [round(new_pos(1)) round(new_pos(2)) round(new_pos(3))];

    % Display the values (or use them in your function)
    disp(['new trans_pos: ', mat2str(trans_pos)]);
    disp(['new target: ', mat2str(target)]);

    % Example of usage: remove existing transducer patches and re-add with new parameters
    remove_transducer_patches();
    add_transducer_shape(trans_pos([2 1 3])/2, target([2 1 3])/2, 'TargetType', 'line'); % TODO cosmetic: read out the target type
end