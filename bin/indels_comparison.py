#!/opt/conda/bin/python

import sys
import click
from scipy.stats import chi2
import pandas as pd
import numpy as np


@click.command()
@click.option('--sample', type=str)
@click.option('--filename', type=click.Path())
def cli(sample, filename):
    '''
    Computes the proportion of mutated genomes per sample, gene-by-gene and also panel-wide.
    It uses the VAF estimates from each unique mutation, then applies the probabilistic principle of inclusion-exclusion.
    The main assumption is that mutations at different sites are statistically independent.
    Output:
    - Lower bound proportion of mutated genomes gene-by-gene and panel-wide.
    - Expected proportion of mutated genomes gene-by-gene and panel-wide.
    - Upper bound proportion of mutated genomes gene-by-gene and panel-wide.
    '''

    maf = load_maf(filename)


    # filter and somatic only
    indel_maf_df = maf.loc[
                            (maf["TYPE"].isin(["INSERTION", "DELETION"]))
                            ].reset_index(drop = True)

    # count total number of mutations (this includes all genes sequenced, not only panel genes)
    indels_panel_df = indel_maf_df.groupby(["SYMBOL", "INDEL_INFRAME", 'Protein_affecting'], dropna = False).size().to_frame('number_mutations').reset_index()
    indels_panel_df = indels_panel_df.pivot(index='SYMBOL', columns=["INDEL_INFRAME", "Protein_affecting"], values='number_mutations')
    indels_panel_df.columns = [f"{x}.{y}" for x, y in indels_panel_df.columns]

    # Ensure all expected columns are present
    expected_columns = ['False.non_protein_affecting', 'True.non_protein_affecting',
                        'False.protein_affecting', 'True.protein_affecting']
    for col in expected_columns:
        if col not in indels_panel_df.columns:
            indels_panel_df[col] = 0

    indels_panel_df = indels_panel_df[expected_columns].copy().fillna(0).astype(int)
    indels_panel_df.columns = ['NOT_MULTIPLE_3.non_protein_affecting', 'MULTIPLE_3.non_protein_affecting',
                        'TRUNCATING.protein_affecting', 'NON_TRUNCATING.protein_affecting']
    indels_panel_df = indels_panel_df.reset_index()


    # Define "new gene" as the combination of all genes
    # This could also be done for certain gene groups
    indels_panel_df_all = pd.DataFrame(indels_panel_df[['NOT_MULTIPLE_3.non_protein_affecting', 'MULTIPLE_3.non_protein_affecting',
                                            'NON_TRUNCATING.protein_affecting', 'TRUNCATING.protein_affecting']].sum()).T
    indels_panel_df_all["SYMBOL"] = 'ALL_GENES'
    indels_panel_df = pd.concat((indels_panel_df, indels_panel_df_all), axis = 0)


    indels_panel_df["Npa_NM3/M3"] = indels_panel_df['NOT_MULTIPLE_3.non_protein_affecting'] / indels_panel_df['MULTIPLE_3.non_protein_affecting']
    indels_panel_df["pa_TRUNC/NOTTRUNC"] = indels_panel_df['TRUNCATING.protein_affecting'] / indels_panel_df['NON_TRUNCATING.protein_affecting']
    indels_panel_df["pa/Npa"] = indels_panel_df['pa_TRUNC/NOTTRUNC'] / indels_panel_df['Npa_NM3/M3']

    indels_panel_df["Npa_NM3/M3"] = indels_panel_df["Npa_NM3/M3"].fillna(-1)
    indels_panel_df["pa_TRUNC/NOTTRUNC"] = indels_panel_df["pa_TRUNC/NOTTRUNC"].fillna(-1)

    # if there is anything missing fill the ratio with 1
    indels_panel_df["pa/Npa"] = indels_panel_df["pa/Npa"].fillna(1)


    indels_panel_df[["g_test_score", "pvalue"]] = [gtest(np.array([w, x]), np.array([y, z])) for w, x, y, z in indels_panel_df[['MULTIPLE_3.non_protein_affecting', 'NOT_MULTIPLE_3.non_protein_affecting',
                                            'NON_TRUNCATING.protein_affecting', 'TRUNCATING.protein_affecting']].values]

    indels_panel_df.to_csv(f"{sample}.sample.indels.tsv",
                            sep="\t",
                            header=True,
                            index=False
                            )



def gtest(null_vector, alt_vector):

    """
    Likelihood-ratio test comparing the proportion point estimate inferred
    from the frequencies of observed species (maximum-likelihood) versus
    the proportion inferred from a null reference.

    This implementation does not require that the count in the null
    reference is the same as in the observed.
    """

    null_prob = null_vector / null_vector.sum()
    expected = sum(alt_vector) * null_prob

    # G-statistic must skip zero cells
    # because lim x->0 xlog(x) = 0
    args = np.where(alt_vector > 0)[0]

    # G-test statistic
    g = 2 * np.sum(alt_vector[args] * np.log(alt_vector[args] / expected[args]))

    # G-test statistic follows chi-squared distribution with dof = {number of categories - 1}
    dof = len(alt_vector) - 1

    # one sided p-val = 1 - CDF(x > g)
    pval = chi2.sf(g, df=dof)

    return g, pval



def load_maf(maf_file):
    '''
    Loads MAF and applies filters
    '''

    maf = pd.read_csv(maf_file, sep='\t', header=0, na_filter = False)

    return maf



if __name__ == '__main__':

    cli()

