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
parameters = load_parameters('nico_test_double_acoustic_100mm_same_temp0_config.yaml');

filepath = fullfile(parameters.data_path, 'sim_outputs');

sbj_IDs = [1 2 8];
data = repmat(struct(), 1, length(sbj_IDs));



i = 1;
for sbj_ID = sbj_IDs%1:6

    disp(['subject ' num2str(sbj_ID)]);

    % try
    %     filename = sprintf('sub-%03d/sub-%03d_layered_heating_resL+z--r_R+z--r_it1_heatingtimeline_imprecisionnone.mat', sbj_ID, sbj_ID);
    %     load(fullfile(filepath, filename), 'time_status_seq', 'tissue_CEM43', 'tissue_heat');
    % catch
    %     filename = sprintf('sub-%03d/sub-%03d_layered_heating_resL--r_R--r_it1_heatingtimeline_imprecisionnone.mat', sbj_ID, sbj_ID);
    %     load(fullfile(filepath, filename), 'time_status_seq', 'tissue_CEM43', 'tissue_heat');
    % end
    filename_glob = sprintf('sub-%03d/sub-%03d_layered_heating_resL*--r_R*--r_it5_heatingtimeline_imprecisionnone.mat', sbj_ID, sbj_ID);
    f = dir(fullfile(filepath, filename_glob));
    % filename = 'sub-008/sub-008_layered_heating_resL-r_R-r_it4_heatingtimeline_debug_imprecisionnone.mat';
    load(fullfile(f.folder, f.name), 'time_status_seq', 'tissue_CEM43', 'tissue_heat');
    % load('M:\Documents\scans\sim_outputs\sub-008\sub-008_layered_heating_resL-r_R-r_it4_heatingtimeline_debug_imprecisionnone.mat');
    
    data(i).time = [time_status_seq.time];
    
    tissue_CEM43 = tissue_CEM43(:,1:end-4);
    tissue_CEM43(tissue_CEM43 <= 0) = NaN;
    data(i).tissue_heat = tissue_heat(:,1:end-4);
    % data(i).tissue_heat = data(i).tissue_heat - 37;
    % data(i).tissue_heat(data(i).tissue_heat <= 0) = NaN;
    data(i).tissue_CEM43_max = max(tissue_CEM43, [], 1);

    i = i + 1;
    
end

figure; hold on;
for i = 1:length(sbj_IDs)
    % plot(data(i).time*2, data(i).tissue_heat(3,:)); % TODO: check why times 2 needed!!!
    plot(data(i).time*2, data(i).tissue_CEM43_max); % TODO: check why times 2 needed!!!
end