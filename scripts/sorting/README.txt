===================
SORTING CONFIG DATA
===================

"create_total_yaml" not really needed; almost all vars are already in there
(from the docstring of the script single_subject_pipeline:
"The file 'default_config' contains all possible to-be altered parameters, but not all of them need to be used to succesfully run the pipeline.")
a few though are not (printed by "find_nonstandard_yaml_keys"):
	run_simulations_with_qsub
        transducer.source_phase_rad
        transducer.pos_grid
        output_location
        localite_instr_file_template
hence, either make sure they are implemented everywhere or drop them

also, TODO make sure there are no var corpses (includes overwriting of params in the scripts)