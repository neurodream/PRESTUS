import os
import pandas as pd
import glob
from datetime import datetime

# Define source base directory
src_base = r"M:\Documents\scans\sim_outputs"

# Find all subdirectories matching the pattern "sub-###"
sub_dirs = [d for d in os.listdir(src_base) if os.path.isdir(os.path.join(src_base, d)) and d.startswith("sub-") and d[4:].isdigit()]

# Dictionary to store data from all files
merged_data = {}

# Process each subdirectory
for sub_dir in sub_dirs:
    src_dir = os.path.join(src_base, sub_dir)
    
    # Find matching CSV files
    csv_files = glob.glob(os.path.join(src_dir, "sub*layered_output_table*.csv"))
    
    for file in csv_files:
        df = pd.read_csv(file, header=None, skiprows=1)  # Skip first row
        var_names = df.iloc[:, 0].tolist()
        values = df.iloc[:, 1].tolist()
        
        # Determine category based on filename
        if any(pattern in file for pattern in ["L-r", "L+z-r"]):
            category = "sham"
        elif any(pattern in file for pattern in ["L--r", "L+z--r"]):
            category = "stim"
        else:
            continue  # Skip files that don't match expected patterns
        
        for var, val in zip(var_names, values):
            key = f"{category}_{var}"
            if key not in merged_data:
                merged_data[key] = {}
            merged_data[key][sub_dir] = val

# Convert to DataFrame
merged_df = pd.DataFrame.from_dict(merged_data, orient='index')

# Drop the first row in the merged file
merged_df = merged_df.iloc[1:]

# Generate output file name with date
date_str = datetime.now().strftime("%y-%m-%d")
output_file = fr"M:\Documents\repos\PRESTUS_forked\data\merged_output_table_{date_str}.csv"

# Save merged data to CSV
merged_df.to_csv(output_file)
print(f"Merged data saved to {output_file}")
