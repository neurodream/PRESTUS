% close all;
clear; clc;

cd(fileparts(mfilename('fullpath')));
cd ..;

% add paths
addpath('functions')
addpath(genpath('toolboxes'))

% % hide confirmation dlgs
% set(groot, 'ShowHiddenHandles', 'on');
% delete(findall(groot, 'Type', 'figure', 'Name', 'Confirm'));

parameters = load_parameters('default_paths_windows.yaml');

plot_scalp = true;
plot_skull = false;
plot_intensity = true;
save = false;

data_path = parameters.data_path;
seg_path = parameters.seg_path;

sub_ids = [1];
iterations = [3 4];
found_iterations = [];

folders = dir(fullfile(data_path, 'sim_outputs', 'sub-*')); % change back
folders = folders(~ismember({folders.name}, {'.', '..'})); % TODO ? necessary?
for i = 1:numel(folders)
    sub_id = sscanf(folders(i).name, 'sub-%d');
    if ismember(sub_id, sub_ids)
        files = dir(fullfile(data_path, 'sim_outputs', folders(i).name, '*parameters*it*'));

        % if numel(files) > 2
        %     files = files(~contains({files.name}, 'parallel'));
        % end
        for j = 1:numel(files)
            file = fullfile(files(j).folder, files(j).name);
            % file = fullfile(files(1).folder, files(1).name);
            tokens = regexp(file, '_it(\d+)_', 'tokens');
            if ~isempty(tokens)
                iteration = str2double(tokens{1}{1}); % Extract number
                if ismember(iteration, iterations) && ~ismember(iteration, found_iterations)
                    parameters = load(file);
                    parameters = parameters.parameters;
                    parameters.data_path = data_path;
                    parameters.seg_path = seg_path;
                    plot_transducer_pos(parameters, sub_id, save, 'Functional', 'intensity', 'LowCutoff', 3);
                    title(['sub' num2str(sub_id) ', it' num2str(iteration)]);
                    found_iterations(end+1) = iteration;
                end
            % break; % only plot the first match
            end
        end
    end

end