var = [...
    1,2,3; ...
    2,3,4; ...
    3,4,5 ...
    ];

for i = 1:size(var, 1)
    var_part = var(i, :);
end

var_comp = [3 4 5];

disp(num2str(var_comp))