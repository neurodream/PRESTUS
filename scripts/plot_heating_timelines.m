close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(currentFile);
cd(rootpath); % repos/PRESTUS_forked/
cd ..

% TODO allow for sham computation

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

% base config ("hard" params)
path_params = load_parameters('nico_test_double_acoustic_same_temp0_config.yaml');
parameters = load('../../scans/sim_outputs/sub-009/sub-009_parametersL--r_R--r_pre-pilot_it19_imprecisionnone.mat', 'parameters');
parameters = parameters.parameters;
parameters.data_path = path_params.data_path;
%load_parameters('nico_test_double_acoustic_same_temp0_config.yaml');

filepath = fullfile(parameters.data_path, 'sim_outputs');

sbj_ID = 9;

disp(['subject ' num2str(sbj_ID)]);

filename_glob = sprintf('sub-%03d/sub-%03d_layered_heating_resL--r_R--r_pre-pilot_it19_imprecisionnone.mat', sbj_ID, sbj_ID);
f = dir(fullfile(filepath, filename_glob));
load(fullfile(f.folder, f.name), 'time_status_seq', 'tissue_CEM43', 'tissue_heat');

% filename_glob = sprintf('sub-%03d/sub-%03d_layered_heating_resL*--r_R*--r_it8_heatingtimeline_imprecisionnone.mat', sbj_ID, sbj_ID);
% filename_glob = sprintf('sub-%03d/sub-%03d_layered_heating_resL+z--r_R+z--r_pre_pilot_imprecisionnone.mat', sbj_ID, sbj_ID);
% f = dir(fullfile(filepath, filename_glob));
% filename = 'sub-008/sub-008_layered_heating_resL-r_R-r_it4_heatingtimeline_debug_imprecisionnone.mat';
% load(fullfile(f.folder, f.name), 'time_status_seq', 'tissue_CEM43', 'tissue_heat');
% load('M:\Documents\scans\sim_outputs\sub-008\sub-008_layered_heating_resL-r_R-r_it4_heatingtimeline_debug_imprecisionnone.mat');

time = [time_status_seq.time];
figure; hold on;
for i = 1:6
    tissue_CEM43_max = tissue_CEM43(i,:);
    to_remove = numel(tissue_CEM43_max) - numel(time);
    plot(time*2, tissue_CEM43_max(1:end-to_remove)); % TODO: check why times 2 needed!!!
end
xlabel('time (seconds) since stim onset');
ylabel('CEM43');
legend(fieldnames(parameters.medium));

data(i).time = [time_status_seq.time];

tissue_CEM43 = tissue_CEM43(:,1:end-4);
tissue_CEM43(tissue_CEM43 <= 0) = NaN;
data(i).tissue_heat = tissue_heat(:,1:end-4);
data(i).tissue_heat = data(i).tissue_heat - 37;
data(i).tissue_heat(data(i).tissue_heat <= 0) = NaN;
data(i).tissue_CEM43_max = max(tissue_CEM43, [], 1);

figure; hold on;
for i = 1:3%length(sbj_IDs)
    plot(data(i).time*2, data(i).tissue_CEM43_max); % TODO: check why times 2 needed!!!
end

% TR
figure; hold on;
for i = 1:3%length(sbj_IDs)
    plot(data(i).time*2, data(i).tissue_heat(3,:)); % TODO: check why times 2 needed!!!
end
ylim([37 39]);