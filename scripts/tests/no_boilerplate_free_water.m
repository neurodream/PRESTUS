% basically replicating the single subject pipeline

%% load kWave functions

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..','..'))
rootpath = pwd;

pn.tuSIM_tools = fullfile(rootpath, 'toolboxes');
addpath(genpath(pn.tuSIM_tools));
pn.kwave = fullfile(rootpath, 'toolboxes', 'k-wave-toolbox-version-1.4', 'k-Wave');
addpath(pn.kwave);

%% 