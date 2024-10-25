import os
import yaml
from collections import defaultdict

def merge_dictionaries(original, new):
    """Merge two dictionaries, updating values intelligently."""
    for key, value in new.items():
        if key in original:
            if isinstance(original[key], dict) and isinstance(value, dict):
                merge_dictionaries(original[key], value)
            else:
                # If the key exists but is not a dictionary, or the new value is not a dictionary, replace or ignore.
                original[key] = value
        else:
            original[key] = value

def load_yaml_files(directory):
    """ Load all YAML files recursively from the specified directory """
    all_data = {}
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith('.yaml'):
                filepath = os.path.join(root, filename)
                with open(filepath, 'r') as file:
                    data = yaml.safe_load(file)
                    if data is not None:
                        merge_dictionaries(all_data, data)

    return all_data

def save_to_yaml(data, output_file):
    """ Save the combined data to a YAML file """
    with open(output_file, 'w') as file:
        yaml.safe_dump(data, file, default_flow_style=None, sort_keys=False)

# Specify the directory containing YAML files and the output file
directory = '../../configs'
output_file = '../../configs/all_configs.yaml'

# Load, combine, and save YAML files
combined_data = load_yaml_files(directory)
save_to_yaml(combined_data, output_file)
