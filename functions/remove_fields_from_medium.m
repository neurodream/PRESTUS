function s = remove_fields_from_medium(s, fieldNames)
    if ischar(fieldNames)
        fieldNames = {fieldNames};
    end
    
    if isstruct(s)
        % Get the fields of the main struct
        fields = fieldnames(s);
        % Iterate over each field of the main struct
        for i = 1:numel(fields)
            if isstruct(s.(fields{i}))
                % Remove specified fields from the second layer (sub-structs)
                for j = 1:numel(fieldNames)
                    if isfield(s.(fields{i}), fieldNames{j})
                        s.(fields{i}) = rmfield(s.(fields{i}), fieldNames{j});
                    end
                end
            end
        end
    end
end