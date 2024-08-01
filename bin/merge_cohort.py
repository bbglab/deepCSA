#!/usr/local/bin/python

import os, sys
import pandas as pd
from read_utils import custom_na_values


samp_name = sys.argv[1]

# Directory containing the .tsv.gz files
directory = '.'

# Get a list of file names in the directory ending with .tsv.gz
file_names = [file for file in os.listdir(directory) if file.endswith('.tsv.gz')]

# Initialize an empty list to store DataFrames
dfs = []

# Read each .tsv.gz file into a DataFrame and append to the list
for file_name in file_names:
    file_path = os.path.join(directory, file_name)
    df = pd.read_csv(file_path, compression='gzip', header = 0, sep='\t', na_values = custom_na_values)  # Read gzipped TSV
    dfs.append(df)

# Concatenate the list of DataFrames into a single DataFrame
concatenated_df = pd.concat(dfs, ignore_index=True)

# 'concatenated_df' now contains the concatenated DataFrame


concatenated_df.to_csv(f"{samp_name}.cohort.tsv.gz",
                        sep = "\t",
                        header = True,
                        index = False)


