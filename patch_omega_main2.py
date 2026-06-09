import sys

with open('subworkflows/local/omega/main.nf', 'r') as f:
    lines = f.readlines()

new_lines = []
for line in lines:
    if "site_comparison_results_flattened }" in line:
        pass # Handle this differently
    
    new_lines.append(line)

# Let's just do it directly using a simple sed-like replacement in python

