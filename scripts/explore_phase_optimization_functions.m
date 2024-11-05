
close all; clear; clc;

cd(fileparts(mfilename('fullpath')));
cd ..

% add paths
addpath('functions')
addpath(genpath('toolboxes'))

%% adjust params here:

% base config ("hard" params)
parameters = load_parameters('nico_test_double_acoustic_100mm_config.yaml');

expected_focal_distances_mm = [70,75,80,85,90,95,100,105];
ROI_widths_mm = [5 10 15 20 25 30 35];

% Define SLURM job parameters
job_name = 'simulation_job';
output_file = '/home/sleep/nicade/Documents/scans/sim_outputs/output.log';
error_file = '/home/sleep/nicade/Documents/scans/sim_outputs/error.log';
partition = 'gpu';
time = '00:05:00';
cpus_per_task = 1;
memory = '2G';
num_gpus = 1; % Number of GPUs requested

% Create the SLURM batch script as a temporary file
slurm_script = tempname; % Generates a unique temporary filename

% Write the SLURM batch script
fid = fopen(slurm_script, 'w');
fprintf(fid, '#!/bin/bash\n');
fprintf(fid, '#SBATCH --job-name=%s\n', job_name);
fprintf(fid, '#SBATCH --output=%s\n', output_file);
fprintf(fid, '#SBATCH --error=%s\n', error_file);
fprintf(fid, '#SBATCH --partition=%s\n', partition);
fprintf(fid, '#SBATCH --time=%s\n', time);
fprintf(fid, '#SBATCH --cpus-per-task=%d\n', cpus_per_task);
fprintf(fid, '#SBATCH --mem=%s\n', memory);
fprintf(fid, '#SBATCH --gres=gpu:%d\n', num_gpus); % Request GPU
fprintf(fid, 'module load matlab\n'); % Load MATLAB module

filepath = fullfile(pwd, 'data'); % 'home/nicade/';
filename = 'simulation_grid';

% matlab_command = 'matlab -nodisplay -nosplash -batch "store_simulated_profiles_grid(XXX, filepath, filename, [1,2], [1,2]); exit;"';
depths_arg = regexprep(mat2str(expected_focal_distances_mm), ' ', ',');
widths_arg = regexprep(mat2str(ROI_widths_mm), ' ', ',');
parameters_fname = 'nico_test_double_acoustic_100mm_config.yaml';
save('parameters.mat', 'parameters');
matlab_command = ['matlab -nodisplay -nosplash -batch "addpath(''functions'');addpath(genpath(''toolboxes''));store_simulated_profiles_grid(''' parameters_fname ''',''' filepath ''',''' filename ''',' depths_arg ',' widths_arg ');exit;"'];
disp(matlab_command);

% store_simulated_profiles_grid(parameters_fname, filepath, filename, expected_focal_distances_mm, ROI_widths_mm);

fprintf(fid, '%s\n', matlab_command);

% Close the script file
fclose(fid);

% Submit the job to SLURM
system(['sbatch ' slurm_script]);

% Optional: Clean up the temporary SLURM script
delete(slurm_script);