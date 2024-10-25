clear all; clc; close all;

foldername = 'sbj_3_duty_45_ntrl_5_amp_500000_targ_128-142-152_tpos_212-152-164';
rootpath = 'M:\Documents\repos\PRESTUS_forked';

disp('loading data... (takes around 30 seconds)')
load(fullfile('D:\', foldername, '\sub-003_layered_results.mat'));
disp('... loading data done.')

%% add relevant paths (copy pasted here; don't need most of it actually)
% TODO why is pn != parameters? where are params loaded in from yaml?
% pn.tuSIM = fullfile(rootpath, 'tools', 'PRESTUS');                addpath(pn.tuSIM); % why PRESTUS within PRESTUS??
pn.tuSIM = rootpath;
pn.scripts = fullfile(pn.tuSIM, 'scripts');                       addpath(pn.scripts);
pn.tuSIM_fun = fullfile(pn.tuSIM, 'functions');                   addpath(pn.tuSIM_fun);
pn.tuSIM_tools = fullfile(pn.tuSIM, 'toolboxes');                 addpath(genpath(pn.tuSIM_tools));
pn.kwave = fullfile(rootpath, 'toolboxes', 'k-wave-toolbox-version-1.4', 'k-Wave');
                                                                  addpath(pn.kwave);
pn.minimize = fullfile(rootpath, 'toolboxes', 'FEX-minimize');    addpath(pn.minimize);
pn.configs = fullfile(rootpath, 'configs');                       % why no add path??
pn.profiles = fullfile(rootpath, 'data', 'transducer_profiles');  % why no add path??
pn.data_path = fullfile(rootpath, '..', '..', 'scans');           % why no add path?? % seems like it is "only" output path? or also loading tpars etc.?
% previously: "data_higher_angle"
% it seems like data are saved within the "data_sims" folder
pn.seg_path = fullfile(rootpath, '..', '..', 'scans', 'segmentation_results');         % why no add path??
pn.nifti = (fullfile(rootpath, 'tools', 'nifti_toolbox'));        addpath(pn.nifti); % does not exist
pn.sim_path = fullfile(rootpath, 'data', 'transducer_pos'); % TODO: careful: transducer name not specified in here yet!

cd(fullfile(rootpath)) % needs to contain a folder called "configs" in this setup

%% run heating sims

% obtain target
posStr = regexp(foldername, 'tpos_(\d+)-(\d+)-(\d+)', 'tokens');
pos = str2double(posStr{1});

% convert the relevant arrays to CPU arrays
kwave_medium.density        = gather(kwave_medium.density);
kwave_medium.sound_speed    = gather(kwave_medium.sound_speed);
sensor_data.p_max_all       = gather(sensor_data.p_max_all);

% ensure no GPU is used
parameters.code_type = 'matlab_cpu';

% set trials to very small number to speed things up a bit
parameters.thermal.n_trials = 2;

run_heating_simulations(sensor_data, kgrid, kwave_medium, sensor, source, parameters, pos);