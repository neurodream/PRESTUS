import os
import yaml

def extract_keys(data, prefix=''):
    """ Recursively extract all keys from nested dictionaries with path-like keys """
    keys = set()
    for key, value in data.items():
        current_path = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            keys.update(extract_keys(value, current_path))
        else:
            keys.add(current_path)
    return keys

def load_yaml_file(filepath):
    """ Load a single YAML file and return its data """
    with open(filepath, 'r') as file:
        return yaml.safe_load(file)
    
def print_keys(additional_keys):
    output = "\n"
    for k in additional_keys:
        output += "\t"
        output += k
        output += '\n'
    return output

def find_additional_keys(default_keys, directory):
    """ Recursively find and print keys in YAML files not present in default_keys """
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith('.yaml') and not filename.endswith('default_config.yaml'):
                filepath = os.path.join(root, filename)
                data = load_yaml_file(filepath)
                if data is not None:
                    file_keys = extract_keys(data)
                    additional_keys = file_keys - default_keys
                    if additional_keys:
                        print(f'{filepath}: {print_keys(additional_keys)}')

# Load the default config YAML
default_config_path = '../../configs/default_config.yaml'
default_config_data = load_yaml_file(default_config_path)

# Assuming the default configuration data is a dictionary
default_keys = extract_keys(default_config_data)

# Specify the directory containing the other YAML files
directory = '../../configs'

# Find and print additional keys
find_additional_keys(default_keys, directory)