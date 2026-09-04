import warnings
warnings.filterwarnings("ignore")

import os
import tqdm
import functools
import sys
import click

import numpy as np
import pandas as pd
from pathlib import Path

PROTEIN_AFFECTING_SET = {
    'nonsense',
    'missense',
    'essential_splice',
    'protein_altering_variant',
    'transcript_amplification',
}

def prob_min_uniform_sample_below_cut(N, n, cut):
    """
    Probability that the minimum of a sampling of n elements from [N] is lower or equal than cut
    """
    arr = np.array([np.log(N - cut - k) - np.log(N - k) for k in range(n)])
    return 1 - np.exp(np.sum(arr))

def prob_min_uniform_sample_below_cut_vec(N, n, cut):
    """
    Probability that the minimum of a sampling of n elements from [N] is lower or equal than cut
    Vectorized version for speed up.
    Computed for arrays of (N, n, cut) triples.
    """
    N = np.asarray(N, dtype=int)
    n = np.asarray(n, dtype=int)
    cut = np.asarray(cut, dtype=int)

    out = np.zeros_like(N, dtype=float)
    out[cut >= N] = 1.0
    valid = (N > 0) & (n > 0) & (cut > 0) & (cut < N)

    for j in np.where(valid)[0]:
        k = np.arange(n[j])
        log_terms = np.log(N[j] - cut[j] - k) - np.log(N[j] - k)
        out[j] = 1 - np.exp(log_terms.sum())

    return out

def load_mutations(somatic_mutations_file, protein_affecting_type):
    somatic_mutations = pd.read_csv(somatic_mutations_file, sep='\t', low_memory=False)
    mutations = somatic_mutations[
        ~(somatic_mutations['FILTER'].str.contains("not_in_exons"))
        & (somatic_mutations['canonical_Protein_affecting'] == protein_affecting_type)
        & (somatic_mutations['TYPE'] == 'SNV')
    ]

    mutations_lite = mutations[['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID', 'ALT_DEPTH', 'ALT_DEPTH_AM', 'SYMBOL', 'canonical_Consequence_broader']]
    mutations_lite = mutations_lite.groupby(by=['CHROM', 'POS', 'REF', 'ALT']).agg({'ALT_DEPTH_AM': 'sum', 'ALT_DEPTH': 'sum'}).reset_index()
    return mutations_lite


def get_aachange_format(r):

    if r['Protein_position'] == '-':
        return '-'
    elif len(r['Amino_acids']) == 1:
        return r['Amino_acids'] + r['Protein_position'] + r['Amino_acids']
    else:
        return r['Amino_acids'].split('/')[0] + r['Protein_position'] + r['Amino_acids'].split('/')[1]

        
def collect_vep(vep_filename):
    vep_panel = pd.read_csv(vep_filename, sep='\t')
    dg = vep_panel[['CHROM', 'POS', 'REF', 'ALT', 'Protein_position', 'Amino_acids', 'GENE']].copy()
    dg['POS'] = dg['POS'].astype(int)
    dg['AACHANGE'] = dg.apply(get_aachange_format, axis=1)
    return dg


def load_panel(consensus_panel_file, depths_file, samples_group, vep_annotations):

    df_panel = pd.read_csv(consensus_panel_file, sep='\t')
    df_depth = pd.read_csv(depths_file, sep='\t')

    df_panel = pd.merge(df_panel, df_depth[['CHROM', 'POS', samples_group]], on=['CHROM', 'POS'], how='left')
    df_panel.rename(columns={samples_group: 'DEPTH'}, inplace=True)
    df_panel = pd.merge(df_panel, vep_annotations[['CHROM', 'POS', 'REF', 'ALT', 'AACHANGE', 'GENE']], 
                              left_on=['CHROM', 'POS', 'REF', 'ALT', 'GENE'], 
                              right_on=['CHROM', 'POS', 'REF', 'ALT', 'GENE'],
                              how='left')
    df_panel = df_panel[df_panel['AACHANGE'] != '-']
    df_panel.dropna(inplace=True)
    df_panel['RESIDUE'] = df_panel['AACHANGE'].apply(lambda s: s[:-1])
    return df_panel


@click.command()
@click.option("--somatic-mutations-file", required=True, type=click.Path(exists=True),help="Path to the somatic mutations file")
@click.option("--vep-file", required=True, type=click.Path(exists=True),help="Path to the VEP annotation file")
@click.option("--consensus-panel-file", required=True, type=click.Path(exists=True),help="Path to the exons consensus panel file")
@click.option("--depths-file", required=True, type=click.Path(exists=True),help="Path to the depths file")
@click.option("--residue", is_flag=True, show_default=True, default=False, help="either genomic or residue based sites")
@click.option("--group-name", type=str, default="all_samples", show_default=True, help="Name of the group/sample to be used in the code")
@click.option("--impact", type=str, default="protein_affecting", show_default=True, help="either protein_affecting or non_protein_affecting positions")
def cli(somatic_mutations_file, vep_filename, consensus_panel_file, depths_file, residue, group_name, impact):

    print(f"Analyzing {group_name}")
    mutations = load_mutations(somatic_mutations_file, impact)
    print("Mutations loaded")
    vep = collect_vep(vep_filename)
    df_panel = load_panel(consensus_panel_file, depths_file, group_name, vep)

    print("Panel loaded")
    mutations_lite = pd.merge(mutations, 
                          df_panel[['CHROM', 'POS', 'REF', 'ALT', 'DEPTH', 'GENE', 'AACHANGE']],
                          on=['CHROM', 'POS', 'REF', 'ALT'], 
                          how='left')
    mutations_lite.dropna(axis=0, inplace=True)  # drop nan sites
    mutations_lite = mutations_lite[mutations_lite['AACHANGE'] != '-']
    mutations_lite['RESIDUE'] = mutations_lite['AACHANGE'].apply(lambda s: s[:-1])

    if residue:
        mutations_lite = mutations_lite.groupby(['GENE', 'RESIDUE']).agg({'ALT_DEPTH': 'sum', 'DEPTH': 'mean'}).reset_index()    
    else:
        mutations_lite = mutations_lite.groupby(['GENE', 'POS']).agg({'ALT_DEPTH': 'sum', 'DEPTH': 'mean'}).reset_index()
    
    subsampling_rates = np.logspace(-2, np.log10(0.9), num=20)
    print("Starting subsampling")

    depth = mutations_lite["DEPTH"].to_numpy()
    alt_depth = mutations_lite["ALT_DEPTH"].to_numpy()
    print(mutations_lite.shape)
    print(depth.shape, alt_depth.shape)

    for i, p in tqdm.tqdm(enumerate(subsampling_rates), total=len(subsampling_rates)):
        n = np.floor(p * depth).astype(int)
        mutations_lite[f"UNIQUE_RATE_{i}"] = prob_min_uniform_sample_below_cut_vec(depth, n, alt_depth)
    
    if residue:
        out_file = f"{group_name}_mutations_residue_rates.{impact}.tsv"
    else:
        out_file = f"{group_name}_mutations_genomic_rates.{impact}.tsv"

    mutations_lite.to_csv(out_file, sep="\t", index=False)        

if __name__ == '__main__':
    
    cli()