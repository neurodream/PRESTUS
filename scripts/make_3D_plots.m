close all; clear; clc;

cd(fileparts(mfilename('fullpath')));
cd ..;


% add paths
addpath('functions')
addpath(genpath('toolboxes'))

parameters = load_parameters();

plot_scalp = true;
plot_skull = false;
plot_intensity = true;
save = true;

folders = dir(fullfile(parameters.data_path, 'sim_outputs'));
folders = folders(~ismember({folders.name}, {'.', '..'}));
for i = 6:numel(folders)
    files = dir(fullfile(parameters.data_path, 'sim_outputs', folders(i).name, '*parameters_L*--*_R*--*'));
    if numel(files) > 2
        files = files(~contains({files.name}, 'parallel'));
    end
    for j = 1:numel(files)
        file = fullfile(files(j).folder, files(j).name);
        parameters = load(file);
        parameters = parameters.parameters;
        plot_transducer_pos(parameters, i, plot_scalp, plot_skull, plot_intensity, save)
    end

end