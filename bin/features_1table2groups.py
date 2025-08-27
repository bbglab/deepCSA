#!/usr/bin/env python

import click
import pandas as pd
import json
from read_utils import custom_na_values

separator2character = {
    'tab' : '\t',
    'comma' : ','
}

def reformatnames(name):
    return str(name).lower().replace(' ','').replace('_','').replace('.','').capitalize()


@click.command()
@click.option('--table-filename', required=True, type=click.Path(exists=True), help='Input features table file')
@click.option('--separator', required=True, type=click.Choice(['tab', 'comma']), help='Separator: tab or comma')
@click.option('--unique-identifier', default=None, type=str, help='Unique identifier column name')
@click.option('--groups', default=None, type=str, help='List of columns with grouping information')
def main(table_filename, separator, unique_identifier, groups):

    sep_char = separator2character[separator]
    features_table = pd.read_table(table_filename, header=0, sep=sep_char, na_values=custom_na_values)

    uniq_name = unique_identifier if unique_identifier else "sample"

    # groups may contain lists of lists, but all formatted into a string
    groups_of_interest_init = [group.strip().strip(",").split(",") for group in groups.replace("[", ";;;").replace("]", "").split(";;;")] if groups else []

    groups_of_interest = []
    for comparison in groups_of_interest_init:
        comparison_group = [item for item in comparison if item != '']
        if len(comparison_group) > 0:
            groups_of_interest.append(comparison_group)

    print(f"Unique identifier: {uniq_name}")
    print(f"Groups of interest: ")
    for group in groups_of_interest:
        print(f"{group}", type(group))

    samples_json = { str(name) : [str(name)] for name in features_table[uniq_name] }
    json_groups = {}

    for col_names2group in groups_of_interest:
        dict_sample_groups = features_table.groupby(by=col_names2group).agg({uniq_name : list}).reset_index().to_dict()
        for category in list(dict_sample_groups[uniq_name].keys()):
            dict_sample_groups[uniq_name][category] = [ str(x) for x in dict_sample_groups[uniq_name][category] ]
        for i in dict_sample_groups[col_names2group[0]].keys():
            sample_batch_name = ""
            for cat in col_names2group:
                formatted_category = reformatnames(cat)
                value_of_variable = dict_sample_groups[cat][i]
                if type(value_of_variable) == bool:
                    sample_batch_name += f"{formatted_category}{ 'Yes' if value_of_variable else 'No'}_"
                else:
                    sample_batch_name += f"{formatted_category}{reformatnames(dict_sample_groups[cat][i])}_"
            sample_batch_name = sample_batch_name.strip("_")
            json_groups[sample_batch_name] = dict_sample_groups[uniq_name][i]

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


if __name__ == '__main__':
    main()