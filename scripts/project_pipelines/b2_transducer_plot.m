%% Pipeline to obtain sensible transducer positions.
% It assumes that you have already did the segmentation for your subjects in SimNIBS (can be done through the main pipeline)
% The script uses gpuArray so run it on a CUDA-enabled machine. 

clear all;
% close all;
clc;

global trans_pos target

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..','..'))
rootpath = pwd;

%%
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

%%
transducer_list = {
    %'CTX250-001-010_60.9mm';...
    % 'CTX500-024-010_77.0mm';...
    'IMASONIC_IGT_300-15473_100.0mm';...
    };
% the specification should actually not vary by depth

% mni_targets = struct('BF_right_test1',[107 127 60]);
all_targets = {...
    'NBM_left', ...
    % 'NBM_right'...
    };

all_subjs = 3;

for i_transducer = 1:length(transducer_list)
    transducer_name = transducer_list{i_transducer};
    parameters = load_parameters(['config_',transducer_name,'.yaml']); % '_windows.yaml'
    % Note: load_parameters always requires a default config that will then be overwritten
    if ismac
        % reset paths to local mac install; NOTE: horrible not to have this be dynamic
        parameters.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
    else
        parameters.simnibs_bin_path = fullfile('/home', 'neuromod', 'julkos', 'SimNIBS', 'bin');
    end
    parameters.ld_library_path ="/opt/gcc/7.2.0/lib64";
    parameters.paths_to_add = {pn.kwave};
    parameters.dist_close = 85; % set distance in mm that is close enough
    
    for subject_id = all_subjs

        fprintf('Current subject: %03i\n', subject_id)

        % Make subfolder (if enabled) and check if directory exists
        if isfield(parameters,'subject_subfolder') && parameters.subject_subfolder == 1
            parameters.output_dir = fullfile(pn.sim_path, sprintf('sub-%03d', subject_id));
        else 
            parameters.output_dir = pn.sim_path;
        end

        parameters.output_dir = pn.sim_path;
        
        if ~isfolder(parameters.output_dir)
            mkdir(parameters.output_dir);
        end  

        % TODO why was this previously in inner loop? it should be only
        % subject-specific?
        headreco_folder = fullfile(pn.seg_path, sprintf('m2m_sub-%03d', subject_id));
        filename_segmented_headreco = fullfile(headreco_folder,'final_tissues.nii.gz');
        segmented_img_orig = niftiread(filename_segmented_headreco);
        segmented_img_head = niftiinfo(filename_segmented_headreco);

        for target_idx = 1:length(all_targets)
            target_name = all_targets{target_idx};
            
            %% load results

            locfile = fullfile(parameters.output_dir, ['tpars',sprintf('_sub-%03d_', subject_id),target_name,'.csv']);
            locs = readtable(locfile);

            pixel_size = segmented_img_head.PixelDimensions(1); % only for the case of ernie; otherwise: segmented_img_head.PixelDimensions(1)

            tppf = locs(locs.prop_intersect<0.05,:);
            % % the following 2 lines were commented out:
            % tppf = tppf(tppf.mean_dist_skull <= quantile(tppf.mean_dist_skull, 0.5) & tppf.var_dist_skull <= quantile(tppf.var_dist_skull, 0.1),:);
            % tppf = tppf(tppf.var_dist_skin==min(tppf.var_dist_skin),:);

            % instead of hard selection of optimum: sort, so to allow
            % selecting second, third, ... best options
            % tppf = tppf(tppf.dist_to_target==min(tppf.dist_to_target),:);
            % i = find(locs.idx==tppf.idx(1));
            % best_trans_pos = locs(i,:);

            % sort by minimum
            tppf = sortrows(tppf, 'dist_to_target');

            % find the highest z position (to avoid skull heterogeneity:
            % temporal window and thicker part) within the 100 nearest
            tppf_nearest = tppf(1:100,:);
            tppf_nearest = sortrows(tppf_nearest, 'trans_z', 'descend');
            best_trans_pos = tppf_nearest(1,:);
            



            trans_pos = [18  137  146]; % [37 56 101]; % [48 47 101]; %[18  137  146]; %floor(table2array(best_trans_pos(1,["trans_x","trans_y","trans_z"])));
            target = floor(table2array(best_trans_pos(1,["targ_x","targ_y","targ_z"])));




            % % TODO remove overwriting
            % trans_pos = [212 152 130];
            % target = [128 142 150];

            [t1_x, t1_y, t1_z] = ndgrid(1:size(segmented_img_orig, 1),1:size(segmented_img_orig, 2),1:size(segmented_img_orig, 3));
            coord_mesh.xyz = [reshape(t1_x,[],1) reshape(t1_y,[],1) reshape(t1_z,[],1)];

            % % TODO clear debug overwrite:
            % trans_pos = [200 200 0];
            % target = [150 210 50];

            figure;
            rotation_vector = show_3d_head(segmented_img_orig, ...
                target, ... % TODO before hardcoded to [128 142 150] - from where?
                trans_pos, ...
                parameters, ...
                pixel_size, ...
                coord_mesh.xyz, ...
                [0 0 0],...
                [0,0],...
                0)
            view([90,-90])
            title(target_name)
            hold on;

            % % plot also the skull
            % tissue = segmented_img_orig;
            % tissue(tissue ~= 7) = 0; % see LUT: 1/2: white/gray matter; 7/8: bone; 5: scalp
            % NA_plot_mri_patch(tissue, 0.5, [0.6 0.6 0.6]);

            % % plot also the spongy skull
            % tissue = segmented_img_orig;
            % tissue(tissue ~= 8) = 0; % see LUT: 1/2: white/gray matter; 7/8: bone; 5: scalp
            % NA_plot_mri_patch(tissue, 0.5, [0.5 0.5 0.5]);

            % % for both skull types together:
            % tissue(tissue ~= 7 & tissue ~= 8) = 0;

            % TODO figure out why these transforms are necessary
            add_transducer_shape(trans_pos([2 1 3])/2, target([2 1 3])/2, 'TargetType', 'line');

            % add ncl basalis
            % if strcmp(target_name, 'NBM_right')
            %     NBM = niftiread(sprintf('M:/Documents/scans/segmentation_results/m2m_sub-%03d/NBM_R.nii', subject_id));
            % else
            %     NBM = niftiread(sprintf('M:/Documents/scans/segmentation_results/m2m_sub-%03d/NBM_L.nii', subject_id));
            % end
            % Ds_NBM = smooth3(NBM);
            % NBM_isosurface = isosurface(Ds_NBM,0.5);
            % NBM_isosurface.vertices = NBM_isosurface.vertices(:, [2 1 3]);
            % NBM_isosurface.vertices = NBM_isosurface.vertices / 2;
            % hiso_NBM = patch(NBM_isosurface,...
            %    'FaceColor',[1.0 0.0 0.0],...
            %    'EdgeColor','none');

            % for c = {'L', 'R'}
            %     NBM = niftiread(sprintf(['M:/Documents/scans/segmentation_results/m2m_sub-%03d/NBM_' c{1} '.nii'], subject_id));
            %     Ds_NBM = smooth3(NBM);
            %     NBM_isosurface = isosurface(Ds_NBM,0.5);
            %     % NBM_isosurface.vertices = NBM_isosurface.vertices(:, [2 1 3]);
            %     NBM_isosurface.vertices = NBM_isosurface.vertices / 2;
            %     patch(NBM_isosurface,...
            %        'FaceColor',[1.0 0.0 0.0],...
            %        'EdgeColor','none');
            % end

            view(90, 90);

        end
    end
end