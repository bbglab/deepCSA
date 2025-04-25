#!/usr/bin/env python

import sys
import pandas as pd
import json

## TODO reimplement with click
filename_of_matrices = sys.argv[1]
samples_json_file = sys.argv[2] # either tab or comma


with open(samples_json_file, 'r') as file:
    samples_info = json.load(file)

samples_only_matrix = None
groups_matrix = None

with open(filename_of_matrices, 'r') as file:
    for line in file.readlines():
        filename = line.strip()
        sample_name = filename.split('.')[0]
        sample_data = pd.read_table(filename, sep = '\t', header = 0)
        sample_data.columns = [ "CONTEXT_MUT", sample_name ]
        sample_data = sample_data.set_index("CONTEXT_MUT")

        if sample_name in samples_info.keys():
            if samples_only_matrix is None:
                samples_only_matrix = sample_data
            else:
                # TODO ensure that this concat preserves the original order of the 96 channels
                # and that it is the same one as the one that SigProfilerTools need
                samples_only_matrix = pd.concat((samples_only_matrix, sample_data), axis = 1)
        else:
            if groups_matrix is None:
                groups_matrix = sample_data
            else:
                # TODO ensure that this concat preserves the original order of the 96 channels
                # and that it is the same one as the one that SigProfilerTools need
                groups_matrix = pd.concat((groups_matrix, sample_data), axis = 1)


if groups_matrix.shape[0] > 0:
    groups_matrix.reset_index().to_csv("groups_matrix.tsv", sep = '\t', header = True, index = False)


if samples_only_matrix.shape[0] > 0:
    samples_only_matrix.reset_index().to_csv("samples_matrix.tsv", sep = '\t', header = True, index = False)
