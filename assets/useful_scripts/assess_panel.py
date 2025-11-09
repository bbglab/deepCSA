#!/usr/bin/env python


import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import click


sys.path.append("../../bin")

from utils_context import triplet_context_iterator

# we need a stable background mutation density estimation,
# for this we need at least 3 mutations

def mutation_based_assessment(panel_tsv, expected_mutation_density_in_mb, samples, depth, mutational_profile = None):
    if mutational_profile is None:
        mutational_profile = pd.DataFrame(np.ones(96) / 96)
        mutational_profile.index = triplet_context_iterator()
    mutational_profile = mutational_profile / mutational_profile.mean()
    mutational_profile.columns = ["Probability"]

    size_per_consequence = panel_tsv.groupby(by = ["GENE", "IMPACT", "CONTEXT_MUT"]).size().reset_index(name = "DEPTH")
    size_per_consequence["DEPTH"] *= depth/3

    size_per_consequence = size_per_consequence.merge(mutational_profile,
                                                      left_on = "CONTEXT_MUT",
                                                      right_index = True,
                                                      how = "left")
    size_per_consequence["ADJUSTED_DEPTH"] = size_per_consequence["DEPTH"] * size_per_consequence["Probability"]

    total_adjusted_size = size_per_consequence.groupby(by = ["GENE", "IMPACT"])["ADJUSTED_DEPTH"].sum() / 3
    expected_mutations_per_cnsq = (expected_mutation_density_in_mb * total_adjusted_size * samples) / 1e6

    return expected_mutations_per_cnsq



@click.command()
@click.option('--panel', '-p', 'panel_path', required=True, type=click.Path(exists=True, dir_okay=False),
              help='Path to panel TSV file')
@click.option('--expected-mutation-density', '-e', 'expected_mutation_density', default=1.2, show_default=True,
              type=float, help='Expected mutation density (mutations per Mb)')
@click.option('--samples', '-s', default=100, show_default=True, type=int,
              help='Number of samples to assume')
@click.option('--depth', '-d', default=5000, show_default=True, type=int,
              help='Sequencing depth (used in size adjustment)')
@click.option('--mutational-profile', '-m', 'mutational_profile_path', default=None, type=click.Path(exists=True, dir_okay=False),
              help='Optional mutational profile TSV with CONTEXT_MUT and Probability (or column named like all_samples.all)')
@click.option('--output-prefix', '-o', default=None, help='Optional prefix to write TSV outputs (will create <prefix>*.tsv)')
def main(panel_path, expected_mutation_density, samples, depth, mutational_profile_path, output_prefix):
    """CLI entrypoint for assess_panel.

    Loads the panel TSV and optional mutational profile, runs the mutation-based assessment
    (with and without the provided profile) and optionally writes results to TSV files.
    """
    panel = pd.read_table(panel_path)

    mutational_prof = None
    if mutational_profile_path:
        mutational_prof = pd.read_table(mutational_profile_path)
        # Try to normalize common formats: expect CONTEXT_MUT as index and a column named Probability
        if 'CONTEXT_MUT' in mutational_prof.columns:
            # If the profile already has a Probability column, keep it; otherwise try to rename
            if 'Probability' in mutational_prof.columns:
                mutational_prof = mutational_prof.set_index('CONTEXT_MUT')[['Probability']]
            elif 'all_samples.all' in mutational_prof.columns:
                mutational_prof = mutational_prof.set_index('CONTEXT_MUT').rename(columns={'all_samples.all': 'Probability'})[['Probability']]
            else:
                # Fallback: use all other numeric column(s) and take the first as Probability
                numeric_cols = mutational_prof.select_dtypes(include=[np.number]).columns.tolist()
                if numeric_cols:
                    mutational_prof = mutational_prof.set_index('CONTEXT_MUT')[[numeric_cols[0]]].rename(columns={numeric_cols[0]: 'Probability'})
                else:
                    # If no numeric columns, create uniform profile (will be normalized in function)
                    mutational_prof = None

    # Run assessments
    mutations_per_gene_consequence = mutation_based_assessment(panel, expected_mutation_density, samples, depth)
    mutations_per_gene_consequence_mut_prof = mutation_based_assessment(panel, expected_mutation_density, samples, depth, mutational_prof)


    signatures_based_assessment(panel_tsv, expected_mutation_density_in_mb, samples, depth, mutational_profile)


    # Print a short summary to stdout
    click.echo("Mutation-based assessment (uniform profile) - top rows:")
    click.echo(mutations_per_gene_consequence.head().to_string())
    click.echo("\nMutation-based assessment (provided mutational profile) - top rows:")
    click.echo(mutations_per_gene_consequence_mut_prof.head().to_string())

    # Optionally write outputs
    if output_prefix:
        out1 = f"{output_prefix}.mutations_per_gene_consequence.tsv"
        out2 = f"{output_prefix}.mutations_per_gene_consequence_mutprof.tsv"
        mutations_per_gene_consequence.to_csv(out1, sep='\t')
        mutations_per_gene_consequence_mut_prof.to_csv(out2, sep='\t')
        click.echo(f"Wrote results to: {out1} and {out2}")