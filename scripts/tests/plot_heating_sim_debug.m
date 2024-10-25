% TODO load all the relevant paths
% TODO ensure that windows paths are active, not linux (assuming debugging
% on windows)

% takes a few moments to load
path_to_file = 'D:\CTX500-024-010_77.0mm\sbj_3_duty_23_ntrl_7_amp_105000_targ_128-142-160_tpos_209-146-180';
load(fullfile(path_to_file, 'sub-003_layered_heating_res.mat')); % assuming it is subj 3
parameters = load_parameters('config_IMASONIC_IGT_300-15473_100.0mm.yaml');

parameters.output_dir = 'M:\Documents\repos\PRESTUS_forked\scripts\tests';
parameters.subject_id = 3;
plot_heating_sims(focal_planeT, time_status_seq, parameters, [212 144 172], []);