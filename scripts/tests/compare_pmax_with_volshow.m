
clear all;
filepath = 'D:\';
foldername = 'sbj_3_duty_45_ntrl_10_amp_500000_targ_128-142-150_tpos_212-152-170'; % change this parameter
filename = 'sub-003_layered_results.mat';
disp('loading data... (takes around half a minute)')
load(fullfile(filepath, foldername, filename));
disp('...loading data done.')
p_max_CPU = gather(sensor_data.p_max_all);

% convert to ISPPA
Isppa_map = p_max_CPU.^2./(2*(kwave_medium.sound_speed.*kwave_medium.density)).*1e-4;
volshow(Isppa_map);