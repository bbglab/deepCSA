#!/usr/local/bin/python

import sys
import pandas as pd
import json
from read_utils import custom_na_values

table_filename = sys.argv[1]
separator_input = sys.argv[2] # either tab or comma
json_information = sys.argv[3]

separator2character = {
    'tab' : '\t',
    'comma' : ','
}

def reformatnames(name):
    return name.lower().replace(' ','').replace('_','').replace('.','').capitalize()


separator = separator2character.get(separator_input, separator_input)

features_table = pd.read_table(table_filename, header = 0, sep=separator, na_values = custom_na_values)


with open(json_information, 'r') as file:
    dict_of_columns = json.load(file)

if len(dict_of_columns) > 0:
    uniq_name = dict_of_columns['unique_identifier']
    groups_of_interest = dict_of_columns['groups_of_interest']
else:
    uniq_name = "sample"
    groups_of_interest = []

samples_json = { str(name) : [str(name)] for name in features_table[uniq_name] }

json_groups = {}



for col_names2group in groups_of_interest:

    dict_sample_groups = features_table.groupby(by = col_names2group).agg({uniq_name : list}).reset_index().to_dict()
    for category in list(dict_sample_groups[uniq_name].keys()):
        dict_sample_groups[uniq_name][category] = [ str(x) for x in dict_sample_groups[uniq_name][category] ]

    for i in dict_sample_groups[col_names2group[0]].keys():
        sample_batch_name = ""
        for cat in col_names2group:
            formatted_category = reformatnames(cat)
            value_of_variable = dict_sample_groups[cat][i]
            if type(value_of_variable) == bool:
                sample_batch_name += f"{formatted_category}{ 'Yes' if value_of_variable else 'No'}_"

            # TODO consider adding an option to categorize numerical variables

            else:
                sample_batch_name += f"{formatted_category}{reformatnames(dict_sample_groups[cat][i])}_"

        sample_batch_name = sample_batch_name.strip("_")
        # print(sample_batch_name, len(dict_sample_groups[uniq_name][i]), dict_sample_groups[uniq_name][i])
        json_groups[sample_batch_name] = dict_sample_groups[uniq_name][i]
    # print()

# Write samples json
with open("samples.json", "w") as outfile:
    json.dump(samples_json, outfile)


if len(json_groups) > 0:
    # Write groups json
    with open("groups.json", "w") as outfile:
        json.dump(json_groups, outfile)


# Write all samples and groups json
json_groups.update(samples_json)

json_groups["all_samples"] = [str(x) for x in samples_json.keys()]
with open("all_groups.json", "w") as outfile:
    json.dump(json_groups, outfile)


