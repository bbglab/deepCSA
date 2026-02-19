#!/usr/bin/env python


import click
import pandas as pd
from read_utils import custom_na_values


def compute_dnds_proxy(mutdensity_file, cohort_syn_mutdensities_file, output_file, mode):
    """
    TODO: explain what this function does
    TODO 2: store a log file that is also outputted and can be used to check some basic statistics

    right now the use of mode is not implemented,
    since we only compute one type of synonymous mutation densities.
    """

    mutdensity_df_init = pd.read_csv(mutdensity_file, sep = "\t", header = 0, na_values = custom_na_values)
    all_possible_genes = list(mutdensity_df_init["GENE"].unique())

    cohort_syn_mutdensity_df = pd.read_csv(cohort_syn_mutdensities_file, sep = "\t", header = 0, na_values = custom_na_values)
    cohort_syn_mutdensity_df.columns = ['GENE', 'cohort_synonymous']
    cohort_syn_mutdensity_df = cohort_syn_mutdensity_df.set_index("GENE")

    init_cohort_syn_df = pd.DataFrame(index = all_possible_genes)
    cohort_syn_df = pd.concat((init_cohort_syn_df, cohort_syn_mutdensity_df), axis = 0)
    
    # filling the null mutation densities with the value of the 1st decile
    cohort_syn_df = cohort_syn_df.fillna(cohort_syn_df[~(cohort_syn_df.isna())].quantile(.1))

    mutdensity_df = mutdensity_df_init.merge(cohort_syn_df, on = "GENE")
    for impact in ["missense", "truncating", "nonsynonymous_splice"]:
        mutdensity_df[f"d_{impact}/d_synonymous"] = mutdensity_df[impact] / mutdensity_df["synonymous"]
        mutdensity_df[f"d_{impact}/d_cohort_synonymous"] = mutdensity_df[impact] / mutdensity_df["cohort_synonymous"]

    # summary at all_samples level
    subset_mutdensities = mutdensity_df[(mutdensity_df["SAMPLE"] == 'all_samples')]
    for impact in ["missense", "truncating"]:
        print(subset_mutdensities.sort_values(by=f"d_{impact}/d_synonymous", ascending=False)[
            ["GENE", "SAMPLE", impact, "synonymous", f"d_{impact}/d_synonymous"]
            ].head(10))


    # # summary at sample-level
    # subset_mutdensities = mutdensity_df[(mutdensity_df["SAMPLE"] != 'all_samples')]
    # for impact in ["missense", "truncating"]:
    #     print(subset_mutdensities.sort_values(by=f"d_{impact}/d_synonymous", ascending=False)[
    #         ["GENE", "SAMPLE", impact, "synonymous", f"d_{impact}/d_synonymous"]
    #         ].head(10))

    # TODO implement these different modes if appropriate
    # if mode == 'mutations':
    #     synonymous_mutdensities_genes = synonymous_mutdensities_all_samples[['GENE', 'synonymous']]
    # elif mode == 'mutated_reads':
    #     synonymous_mutdensities_genes = synonymous_mutdensities_all_samples[['GENE', 'synonymous']]

    mutdensity_df.to_csv(f"{output_file}",
                            header=True,
                            index=False,
                            sep="\t")


@click.command()
@click.option('--mutdensities', type=click.Path(exists=True), help='Input mutation density file')
@click.option('--cohort-syn-mutdensities', type=click.Path(exists=True), help='Input cohort synonymous mutation densities')
@click.option('--output', type=click.Path(), help='Output file')
@click.option('--mode', type=click.Choice(['mutations', 'mutated_reads']), default='mutations')
def main(mutdensities, cohort_syn_mutdensities, output, mode):
    click.echo("Selecting the gene synonymous mutation densities...")
    compute_dnds_proxy(mutdensities, cohort_syn_mutdensities, output, mode)

if __name__ == '__main__':
    main()

