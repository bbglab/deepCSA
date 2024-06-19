#!/opt/conda/bin/python

import sys
import click
from scipy.stats import beta
import pandas as pd


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

    res = {
        'SYMBOL': [],    # gene symbol or 'all' for panel-wide
        'LOWER': [],     # lower bound proportion of mutated genomes
        'EXPECTED': [],  # expected proportion of mutated genomes
        'UPPER': [],     # upper bound proportion of mutated genomes
        'SAMPLE': []     # sample ID
        }

    genes = maf['SYMBOL'].unique()
    df = maf.groupby(by=['CHROM', 'POS', 'REF']).agg({'VAF': 'sum', 'VAF_lower': 'sum', 'VAF_upper': 'sum', 'SYMBOL': 'first'}).reset_index()

    # proportion by gene

    for g in genes:

        dh = df[df['SYMBOL'] == g]
        v = dh['VAF'].tolist()
        v_low = dh['VAF_lower'].tolist()
        v_high = dh['VAF_upper'].tolist()

        res['SYMBOL'].append(g)
        res['LOWER'].append(pie(v_low) if len(v) > 0 else 0)
        res['EXPECTED'].append(pie(v) if len(v) > 0 else 0)
        res['UPPER'].append(pie(v_high) if len(v) > 0 else 0)
        res['SAMPLE'].append(sample)

    # global proportion

    v = df['VAF'].tolist()
    v_low = df['VAF_lower'].tolist()
    v_high = df['VAF_upper'].tolist()

    res['SYMBOL'].append('all')
    res['LOWER'].append(pie(v_low) if len(v) > 0 else 0)
    res['EXPECTED'].append(pie(v) if len(v) > 0 else 0)
    res['UPPER'].append(pie(v_high) if len(v) > 0 else 0)
    res['SAMPLE'].append(sample)

    # save as a dataframe

    res = pd.DataFrame(res)
    res.to_csv(f"{sample}.sample.mutated_genomes_from_vaf.tsv",
               sep="\t",
               header=True,
               index=False
               )


def load_maf(maf_file):
    '''
    Loads MAF and applies filters
    '''

    maf = pd.read_csv(maf_file, sep='\t', header=0)

    # Clopper-Pearson 95% CI
    maf.loc[:, 'VAF_lower'] = get_binomial_low(maf['ALT_DEPTH'].astype(float).values, maf['DEPTH'].astype(float).values)
    maf.loc[:, 'VAF_upper'] = get_binomial_high(maf['ALT_DEPTH'].astype(float).values, maf['DEPTH'].astype(float).values)

    return maf


def pie(plist):
    '''
    Probabilistic version of the principle of inclusion-exclusion
    '''

    plist = list(plist)
    n = len(plist)

    if n > 2:
        p1 = pie(plist[: n // 2])
        p2 = pie(plist[n // 2: ])
        return pie([p1, p2])
    if n == 2:
        return plist[0] + plist[1] - plist[0] * plist[1]
    if n == 1:
        return plist[0]


'''
Compute the confidence interval about phat = x / n
following the Clopper-Pearson method: https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Clopper%E2%80%93Pearson_interval
'''

def get_binomial_low(x, n, alpha=0.05):
    return beta(x, n - x + 1).ppf(alpha/2)


def get_binomial_high(x, n, alpha=0.05):
    return beta(x + 1, n - x).ppf(1 - alpha/2)


if __name__ == '__main__':

    cli()

