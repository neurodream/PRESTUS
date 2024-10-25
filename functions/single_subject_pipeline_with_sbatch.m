function single_subject_pipeline_with_sbatch(subject_id, parameters, timelimit, memorylimit)
    arguments
        subject_id double
        parameters struct
        timelimit (1,1) double = 60*60*4 % time limit for a job in seconds (4 hours by default)
        memorylimit (1,1) double = 40 % memory limit for a job in Gb (40 Gb by default)
    end
    if parameters.interactive
        warning('Processing is set to interactive mode, this is not supported when running jobs with sbatch, switching off interactive mode.')
        parameters.interactive = 0;
    end
    assert(matches(parameters.overwrite_files, ["always", "never"]), "When running jobs with sbatch, it is not possible to create dialog windows to ask for a confirmation when a file already exists. Set parameters.overwrite_files to 'always' or 'never'");

    % Make subfolder (if enabled) and check if directory exists
    if isfield(parameters, 'subject_subfolder') && parameters.subject_subfolder == 1
        output_dir = fullfile(parameters.temp_output_dir, sprintf('sub-%03d', subject_id));
    else
        output_dir = fullfile(parameters.temp_output_dir);
    end
    
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end    

    log_dir = fullfile(output_dir, 'batch_job_logs');
    if ~exist(log_dir, 'dir')
        mkdir(log_dir)
    end

    [path_to_pipeline, ~, ~] = fileparts(which('scripts/single_subject_pipeline'));
    
    subj_id_string = sprintf('sub-%03d', subject_id);
    
    % Save inputs in the temp file
    temp_data_path = tempname(log_dir);
    [tempdir, tempfile] = fileparts(temp_data_path);
    tempfile = [tempfile '.mat'];
    temp_data_path = [temp_data_path '.mat'];
    save(temp_data_path, "subject_id", "parameters");

    temp_m_file = tempname(log_dir);
    fid = fopen([temp_m_file '.m'], 'w+');
    fprintf(fid, "load %s; cd %s; single_subject_pipeline(subject_id, parameters); delete %s; delete %s;", temp_data_path, path_to_pipeline, temp_data_path, [temp_m_file '.m']);
    fclose(fid);
    [~, temp_m_file_name, ~] = fileparts(temp_m_file);
    
    matlab_cmd = sprintf('matlab -batch "%s"', temp_m_file_name);

    if ~isfield(parameters, 'sbatch_job_prefix')
        parameters.sbatch_job_prefix = 'PRESTUS';
    end
    job_name = [parameters.sbatch_job_prefix '_' subj_id_string];

    % Convert time limit from seconds to the format needed by SLURM
    hours = floor(timelimit / 3600);
    minutes = floor((timelimit - hours * 3600) / 60);
    seconds = timelimit - hours * 3600 - minutes * 60;
    timelimit_str = sprintf('%d:%02d:%02d', hours, minutes, seconds);

    sbatch_call = sprintf('bash -c "/usr/bin/sbatch" --job-name=%s --gres=gpu:1 --constraint=cudacap>=5.0 --mem=%iG --time=%s --output=%s --error=%s --chdir=%s', ...
        job_name, memorylimit, timelimit_str, ...
        sprintf('%s_sbatch_pipeline_output_%%j.log', subj_id_string), ...
        sprintf('%s_sbatch_pipeline_error_%%j.log', subj_id_string), ...
        log_dir);

    % Include loading SLURM module or setting PATH
    slurm_setup = 'module load slurm;'; % Modify this if your cluster uses a different command or path setup

    % Construct the full command
    full_cmd = sprintf('%s cd %s; echo ''%s'' | %s', slurm_setup, log_dir, matlab_cmd, sbatch_call);
    
    % Print and submit the job
    fprintf('Submitted the job to the cluster.\nSee logs in %s in case there are errors. \n', log_dir)
    [res, out] = system(full_cmd);
    
    fprintf('Job name: %s; job ID: %s\n', job_name, out)
end
