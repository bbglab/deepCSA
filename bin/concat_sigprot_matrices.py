#!/usr/bin/env python

import click
import pandas as pd
import json


# identify samples without mutations
def remove_zero_mutation_samples(matrix):
    zero_col = matrix.columns[matrix.sum(axis=0) == 0]
    if len(zero_col) > 0:
        print(f"Excluding samples without mutations: {list(zero_col)}")
        matrix = matrix.drop(columns=zero_col)
    return matrix


def save_all_matrices(matrix, type_of_profile, group_meaning):
    samples_only_matrix_init = remove_zero_mutation_samples(matrix)

    # SigProfilerTools
    samples_only_matrix = samples_only_matrix_init.copy()
    samples_only_matrix.reset_index().to_csv(f"{group_meaning}_matrix.{type_of_profile}.sp.tsv", sep='\t', header=True, index=False)
    rounded_samples_only_matrix = remove_zero_mutation_samples(samples_only_matrix.round().astype(int))
    rounded_samples_only_matrix.reset_index().to_csv(f"{group_meaning}_matrix.{type_of_profile}.sp.round.tsv", sep='\t', header=True, index=False)

    # HDP
    samples_only_matrix = samples_only_matrix_init.transpose()
    samples_only_matrix.index.name = None
    hdp_sorted_contexts = sorted(samples_only_matrix.columns, key = lambda x : x[0] + x[2] + x[-1] + x[1:])
    samples_only_matrix = samples_only_matrix[hdp_sorted_contexts]
    samples_only_matrix.to_csv(f"{group_meaning}_matrix.{type_of_profile}.hdp.tsv", sep = '\t', header = True, index = True, quoting=1)

    # MSA
    mutation_matrix = samples_only_matrix_init.reset_index().copy()
    mutation_matrix["Type"] = mutation_matrix["CONTEXT_MUT"].apply(lambda x: x[2:-2])
    mutation_matrix["SubType"] = mutation_matrix["CONTEXT_MUT"].apply(lambda x: x[0] + x[2] + x[-1] )
    mutation_matrix = mutation_matrix.drop(columns = ["CONTEXT_MUT"])
    mutation_matrix = mutation_matrix[
            ["Type", "SubType"] + [col for col in mutation_matrix.columns if col not in ["Type", "SubType"]]
        ].sort_values(by = ["Type", "SubType"])

    mutation_matrix.to_csv(f"{group_meaning}_matrix.{type_of_profile}.msa.csv", sep = ',', index = False)


def concat_sigprot_matrices(filename_of_matrices, samples_json_file, type_of_profile):
    """
    Concatenate signature profile matrices for samples and groups.
    """
    # Load the samples information from the JSON file
    with open(samples_json_file, 'r') as file:
        samples_info = json.load(file)

    samples_only_matrix = None
    groups_matrix = None

    # Read the list of matrix filenames
    with open(filename_of_matrices, 'r') as file:
        for line in file.readlines():
            filename = line.strip()
            sample_name = filename.split('.')[0]
            sample_data = pd.read_table(filename, sep='\t', header=0)
            sample_data.columns = ["CONTEXT_MUT", sample_name]
            sample_data = sample_data.set_index("CONTEXT_MUT")

            # Separate into samples and groups
            if sample_name in samples_info.keys():
                if samples_only_matrix is None:
                    samples_only_matrix = sample_data
                else:
                    # Ensure the concatenation preserves the original order of the 96 channels
                    samples_only_matrix = pd.concat((samples_only_matrix, sample_data), axis=1)
            else:
                if groups_matrix is None:
                    groups_matrix = sample_data
                else:
                    # Ensure the concatenation preserves the original order of the 96 channels
                    groups_matrix = pd.concat((groups_matrix, sample_data), axis=1)


    # Save the groups matrix if it exists
    if groups_matrix is not None and groups_matrix.shape[0] > 0:
        save_all_matrices(groups_matrix, type_of_profile, 'groups')


    # Save the samples matrix if it exists
    if samples_only_matrix is not None and samples_only_matrix.shape[0] > 0:
        save_all_matrices(samples_only_matrix, type_of_profile, 'samples')


@click.command()
@click.option('--filename_of_matrices', type=click.Path(exists=True), required=True, help='File containing the list of matrix filenames.')
@click.option('--samples_json_file', type=click.Path(exists=True), required=True, help='JSON file containing sample information.')
@click.option('--type_of_profile', type=str, required=True, help='Type of profile (e.g.: all, non_protein_affecting, exons, introns, etc.).')
def main(filename_of_matrices, samples_json_file, type_of_profile):
    """
    CLI entry point for concatenating signature profile matrices.
    """
    click.echo(f"Processing matrices from: {filename_of_matrices}")
    click.echo(f"Using sample information from: {samples_json_file}")
    click.echo(f"Profile type: {type_of_profile.strip('.')}")
    concat_sigprot_matrices(filename_of_matrices, samples_json_file, type_of_profile.strip('.'))
    click.echo("Matrices concatenation completed.")


if __name__ == '__main__':
    main()
