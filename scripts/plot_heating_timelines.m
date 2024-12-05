close all; clear; clc;

currentFile = matlab.desktop.editor.getActiveFilename;
rootpath = fileparts(currentFile);
cd(rootpath); % repos/PRESTUS_forked/
cd ..

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

% base config ("hard" params)
parameters = load_parameters('nico_test_double_acoustic_100mm_same_temp0_config.yaml');

filepath = fullfile(parameters.data_path, 'sim_outputs');
figure; hold on;
for sbj_ID = 8%1:6

    disp(['subject ' num2str(sbj_ID)]);

    % try
    %     filename = sprintf('sub-%03d/sub-%03d_layered_heating_resL+z--r_R+z--r_it1_heatingtimeline_imprecisionnone.mat', sbj_ID, sbj_ID);
    %     load(fullfile(filepath, filename), 'time_status_seq', 'tissue_CEM43', 'tissue_heat');
    % catch
    %     filename = sprintf('sub-%03d/sub-%03d_layered_heating_resL--r_R--r_it1_heatingtimeline_imprecisionnone.mat', sbj_ID, sbj_ID);
    %     load(fullfile(filepath, filename), 'time_status_seq', 'tissue_CEM43', 'tissue_heat');
    % end
    filename = sprintf('sub-%03d/sub-%03d_layered_heating_resL--r_R--r_it4_heatingtimeline_imprecisionnone.mat', sbj_ID, sbj_ID);
    filename = 'sub-008/sub-008_layered_heating_resL-r_R-r_it4_heatingtimeline_debug_imprecisionnone.mat';
    % load(fullfile(filepath, filename), 'time_status_seq', 'tissue_CEM43', 'tissue_heat');
    load('M:\Documents\scans\sim_outputs\sub-008\sub-008_layered_heating_resL-r_R-r_it4_heatingtimeline_debug_imprecisionnone.mat');
    
    time = [time_status_seq.time];
    
    tissue_CEM43 = tissue_CEM43(:,1:end-4);
    tissue_CEM43(tissue_CEM43 <= 0) = NaN;
    tissue_heat = tissue_heat(:,1:end-4);
    tissue_heat = tissue_heat - 37;
    tissue_heat(tissue_heat <= 0) = NaN;
    tissue_CEM43_max = max(tissue_CEM43, [], 1);

    % plot(time*2, tissue_heat(3,:)); % TODO: check why times 2 needed!!!
    % plot(time*2, tissue_CEM43_max); % TODO: check why times 2 needed!!!
    
end