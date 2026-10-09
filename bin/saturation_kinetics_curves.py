#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os

import warnings
# warnings.filterwarnings("ignore")

# import tqdm
import click

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions
tfb = tfp.bijectors


CONSEQUENCES_CATEGORIES = {
    'protein_affecting' : { 'nonsense', 'missense', 'essential_splice', 'protein_altering_variant'},
    'nonsense': {'nonsense'},
    'truncating': {'nonsense', 'essential_splice'},
    'missense' : { 'missense' },
    'synonymous' : { 'synonymous' }
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
        k = np.arange(int(n[j]))
        log_terms = np.log(N[j] - cut[j] - k) - np.log(N[j] - k)
        out[j] = 1 - np.exp(log_terms.sum())

    return out

def load_mutations(somatic_mutations_file):
    somatic_mutations = pd.read_csv(somatic_mutations_file, sep='\t', low_memory=False)
    mutations = somatic_mutations[
        ~(somatic_mutations['FILTER'].str.contains("not_in_exons"))
        & (somatic_mutations['TYPE'] == 'SNV')
    ]

    mutations_lite = mutations[['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID', 'ALT_DEPTH', 'ALT_DEPTH_AM']]
    mutations_lite = mutations_lite.groupby(by=['CHROM', 'POS', 'REF', 'ALT']
                                            ).agg({'ALT_DEPTH_AM': 'sum', 'ALT_DEPTH': 'sum'}).reset_index()
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




def empirical_discovery_index_curve(gene, mutations_dict, df_panel_dict,
                                    subsampling_rates,
                                    replicates=100, sites='genomic'):

    df = mutations_dict[sites]
    df = df[df['GENE'] == gene]

    dg = df_panel_dict[sites]
    dg = dg[dg['GENE'] == gene]

    size = dg.shape[0]
    mean_depth = df['DEPTH'].mean()
    
    x, mean, err_low, err_high = [], [], [], []
    
    unique_dict = {}
    # for i, p in tqdm.tqdm(enumerate(subsampling_rates)):
    for i, p in enumerate(subsampling_rates):
        dist_bernoulli = tfd.Bernoulli(probs=df[f'UNIQUE_RATE_{i}'].values)
        unique_mutations = np.sum(dist_bernoulli.sample(sample_shape=(100,)), axis=1)
        y = list(unique_mutations / size)
        mean += [np.mean(y)]
        err_low += [np.percentile(y, 2.5)]
        err_high += [np.percentile(y, 97.5)]
        x += [mean_depth * p]
    mean += [df.shape[0] / size]
    err_low += [df.shape[0] / size]
    err_high += [df.shape[0] / size]
    x += [mean_depth]

    return x, mean, err_low, err_high


def compute_interval_slopes(x, y):
    """
    Rate of change of y across each consecutive interval of x, in log-space:

        slope_i = (y[i+1] - y[i]) / (log10(x[i+1]) - log10(x[i]))

    Returns (midpoints, slopes), where midpoints is the geometric mean of the
    two depth bounds of each interval. Intervals with non-increasing or
    non-finite depths get a NaN slope.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    d_log = np.diff(np.log10(x))
    d_y = np.diff(y)
    slopes = np.full(d_y.shape, np.nan)
    valid = np.isfinite(d_log) & (d_log > 0) & np.isfinite(d_y)
    slopes[valid] = d_y[valid] / d_log[valid]
    midpoints = np.sqrt(x[:-1] * x[1:])
    return midpoints, slopes


def plot_empirical_discovery(gene, mutations_dict, df_panel_dict, subsampling_rates, sites='genomic',
                             impact="protein_affecting",
                             pdf=None):
    click.echo("Plotting empirical discovery")
    x, mean, err_low, err_high = empirical_discovery_index_curve(gene, mutations_dict, df_panel_dict, subsampling_rates, sites = sites)
    fig, ax1 = plt.subplots(figsize=(2,2))

    ax1.scatter(x, mean, s=50)
    for i, m in enumerate(x):
        ax1.vlines(m, err_low[i], err_high[i])
        
    ax1.set_xscale('log')
    if sites == 'residue':
        ax1.set_ylabel('proportion of\nmutated residues')
    elif sites == 'genomic':
        ax1.set_ylabel('proportion of\nmutated nucleotides')
  
    ax1.set_xlabel('depth per residue')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.set_ylim(0, max(err_high) * 1.1)
    plt.title(f"{gene} ({impact}, {sites})")
    if pdf is not None:
        pdf.savefig(fig, bbox_inches='tight', dpi=300)
    plt.show()
    plt.close()


def plot_slope_comparison(gene, x_empirical, slopes_empirical, slopes_theoretical,
                          sites='genomic', impact="protein_affecting", pdf=None):
    """
    Compare the rate of change (log-space slope) of the empirical curve against
    the theoretical neutral curve, computed over the same intervals.
    """
    click.echo(f"Plotting slope comparison for {gene}")
    fig, ax1 = plt.subplots(figsize=(3, 2.5))

    midpoints = np.sqrt(np.asarray(x_empirical[:-1]) * np.asarray(x_empirical[1:]))

    ax1.scatter(midpoints, slopes_empirical, color='brown', s=25, label='empirical')
    ax1.plot(midpoints, slopes_theoretical, color='grey', lw=2, alpha=0.5, label='neutral theoretical')

    ax1.set_xscale('log')
    ax1.set_xlabel('depth per residue')
    ax1.set_ylabel('rate of change\n(Δ proportion / Δ log10 depth)')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.legend(loc='best', fontsize=6)

    plt.title(f"{gene} ({impact}, {sites})")
    if pdf is not None:
        pdf.savefig(fig, bbox_inches='tight', dpi=300)
    plt.show()
    plt.close()



# theoretical neutral vs empirical discovery curves
# Create a PDF to save the plots
def main_empirical(sample,
                   mutations_dict, panel_df,
                   df_panel_dict,
                   omega_mutability_file, relative_mutability_file,
                   subsampling_rates,
                   sites='genomic', impact = "protein_affecting", logscale=False, genes_list = None,
                   empirical_pdf=None, combined_pdf=None, slope_pdf=None):

    # retrieve relative mutability
    mutability_raw = pd.read_csv(relative_mutability_file, sep='\t',
                                 header=None, names=['CHROM', 'POS', 'REF', 'ALT', 'MUTABILITY'])

    mutability_raw = pd.merge(mutability_raw, panel_df, on=['CHROM', 'POS', 'REF', 'ALT'], how='left')
    mutability_raw = mutability_raw.rename(columns={sample: "MUTABILITY"})

    if genes_list is None:
        genes_list = panel_df['GENE'].unique()

    mutations_lite = mutations_dict[sites]
    mutations_lite['VAF'] = mutations_lite.apply(lambda r: r['ALT_DEPTH']/r['DEPTH'], axis=1)

    all_slope_records = []

    # for gene in tqdm.tqdm(genes_list):

    for gene in genes_list:
        try:
            plot_empirical_discovery(gene, mutations_dict, df_panel_dict, subsampling_rates, sites=sites,
                                     impact=impact, pdf=empirical_pdf)

            synonymous_mutation_rate = pd.read_csv(omega_mutability_file, sep='\t')
            synonymous_mutation_rate = synonymous_mutation_rate[synonymous_mutation_rate['GENE'] == gene]
            mutability_gene = mutability_raw[mutability_raw['GENE'] == gene]
            mutability_gene = pd.merge(mutability_gene, synonymous_mutation_rate[['CONTEXT_MUT', sample]], on=['CONTEXT_MUT'], how='left')
            mutability_gene.rename(columns={sample: 'MUTRATE'}, inplace=True)

            # discard positions in non-CDS regions, probably splicing and intronic
            mutability_gene = mutability_gene[(mutability_gene['AACHANGE'] != '-') & (~mutability_gene['AACHANGE'].isnull())]

            # keep only sites of interest
            mutability_gene = mutability_gene[mutability_gene["IMPACT"].isin(CONSEQUENCES_CATEGORIES[impact])]

            mutability_gene['RESIDUE'] = mutability_gene['AACHANGE'].apply(lambda s: s[:-1])
            print(mutability_gene.head())

            if sites == 'residue':
                mutability_gene = mutability_gene.groupby(['GENE', 'RESIDUE']).agg({'MUTRATE': 'sum', 'DEPTH': 'mean'}).reset_index()
            elif sites == 'genomic':
                mutability_gene = mutability_gene.groupby(['GENE', 'POS']).agg({'MUTRATE': 'sum', 'DEPTH': 'mean'}).reset_index()
            
            mutability_gene['MUTABILITY'] = mutability_gene.apply(lambda s: (s['MUTRATE'] / s['DEPTH']), axis=1)

            if sites == 'residue':
                mutability_gene = pd.merge(mutability_gene, mutations_lite[['GENE', 'RESIDUE', 'VAF']], on=['GENE', 'RESIDUE'], how='left')
            elif sites == 'genomic':
                mutability_gene = pd.merge(mutability_gene, mutations_lite[['GENE', 'POS', 'VAF']], on=['GENE', 'POS'], how='left')


            # neutral rate
            mutability_gene['RATE_NEUTRAL'] = mutability_gene['MUTABILITY']
            total_neutral_rate = mutability_gene['RATE_NEUTRAL'].sum()

            # compute saturation theoretical
            y_unique_neutral = []
            if sites == 'genomic':
                x_theoretical = np.logspace(3, 8, num=100)
            elif sites == 'residue':
                x_theoretical = np.logspace(3, 7, num=100)

            for depth in x_theoretical:
                unique_neutral = np.sum(1 - np.exp(-mutability_gene['RATE_NEUTRAL'].to_numpy(dtype=float) * depth)) / mutability_gene.shape[0]
                y_unique_neutral.append(unique_neutral)

            # compute empirical discovery index curve
            x_empirical, mean, err_low, err_high = empirical_discovery_index_curve(gene, mutations_dict, df_panel_dict,
                                                                                   subsampling_rates, sites=sites)

            # rate of change per interval: empirical vs theoretical
            # (theoretical curve interpolated at the empirical depths so both
            # slopes are computed over identical intervals)
            y_theoretical_at_empirical = np.interp(np.log10(x_empirical), np.log10(x_theoretical), y_unique_neutral)
            _, slopes_empirical = compute_interval_slopes(x_empirical, mean)
            _, slopes_theoretical = compute_interval_slopes(x_empirical, y_theoretical_at_empirical)

            slope_records = [
                {
                    'GENE': gene, 'SITES': sites, 'IMPACT': impact,
                    'DEPTH_LOW': x_empirical[i], 'DEPTH_HIGH': x_empirical[i + 1],
                    'DEPTH_MID': np.sqrt(x_empirical[i] * x_empirical[i + 1]),
                    'SLOPE_EMPIRICAL': slopes_empirical[i],
                    'SLOPE_THEORETICAL': slopes_theoretical[i],
                    'SLOPE_RATIO': (slopes_empirical[i] / slopes_theoretical[i]
                                    if np.isfinite(slopes_theoretical[i]) and slopes_theoretical[i] != 0
                                    else np.nan),
                }
                for i in range(len(x_empirical) - 1)
            ]
            all_slope_records.extend(slope_records)

            # plot
            fig, ax1 = plt.subplots(figsize=(2,2))
            ax1.set_xscale('log')
            if logscale:
                ax1.set_yscale('log')

            # empirical
            ax1.scatter(x_empirical[:-1], mean[:-1], label='downsampling', color='brown', s=5)
            ax1.scatter(x_empirical[-1], mean[-1], label='observed', color='white', edgecolors='brown', alpha=1, s=100)
            for i, m in enumerate(x_empirical[:-2]):
                plt.vlines(m, err_low[i], err_high[i], color='brown', lw=1)

            # theoretical

            ax1.plot(x_theoretical, y_unique_neutral, color='grey', lw=2, label='neutral theoretical', alpha=0.5)  # neutral
            # ax1.tick_params(axis='y', labelcolor=color)
            if sites == 'residue':
                ax1.set_ylabel('proportion of\nmutated residues')
            elif sites == 'genomic':
                ax1.set_ylabel('proportion of\nmutated nucleotides')
            
            ax1.set_xlabel('depth per residue')
            # ax1.vlines(5e5, 0, 1., linestyles='dashed', color='maroon', label='cohort', alpha=0.3)

            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)

            ax1.set_ylim(0,1)
            ax1.set_xlim(x_theoretical[0], x_theoretical[-1])

            # ax1.legend(loc=(1,0))

            plt.title(f"{gene} ({impact}, {sites})")

            if combined_pdf is not None:
                combined_pdf.savefig(fig, bbox_inches='tight', dpi=300)

            plt.show()
            plt.close(fig)

            # slope comparison plot
            plot_slope_comparison(gene, x_empirical, slopes_empirical, slopes_theoretical,
                                  sites=sites, impact=impact, pdf=slope_pdf)

        except Exception as e:
            print(f"Error occurred while processing {gene}: {e}")
            continue

    # save the slope table for all genes
    if all_slope_records:
        df_slopes = pd.DataFrame(all_slope_records)
        slope_table_file = f"{sample}_slopes_{sites}.{impact}.tsv"
        df_slopes.to_csv(slope_table_file, sep="\t", index=False)
        click.echo(f"Slope table saved to {slope_table_file}")

    return df_slopes if all_slope_records else None

def compute_mutation_rates(mutations, name, impact, subsampling_rates_list, residue=False):
    """Compute mutation rates for a given set of mutations."""

    grp_by = ['GENE', 'RESIDUE'] if residue else ['GENE', 'POS']
    mutations_lite = mutations.groupby(
        grp_by
        ).agg({'ALT_DEPTH': 'sum', 'DEPTH': 'mean'}).reset_index()    

    depth = mutations_lite["DEPTH"].to_numpy()
    alt_depth = mutations_lite["ALT_DEPTH"].to_numpy()

    click.echo(f"Computing mutation rates for {name}")
    print(mutations_lite.shape)
    print(depth.shape, alt_depth.shape)

    # for i, p in tqdm.tqdm(enumerate(subsampling_rates_list), total=len(subsampling_rates_list)):
    for i, p in enumerate(subsampling_rates_list):
        n = np.floor(p * depth).astype(int)
        mutations_lite[f"UNIQUE_RATE_{i}"] = prob_min_uniform_sample_below_cut_vec(depth, n, alt_depth)
    
    out_file = f"{name}_mutations_{'residue' if residue else 'genomic'}_rates.{impact}.tsv"

    click.echo(f"Saving mutation rates to {out_file}")
    mutations_lite.to_csv(out_file, sep="\t", index=False)

    return mutations_lite




@click.command()
@click.option("--somatic-mutations-file", required=True, type=click.Path(exists=True),help="Path to the somatic mutations file")
@click.option("--vep-file", required=True, type=click.Path(exists=True),help="Path to the VEP annotation file")
@click.option("--consensus-panel-file", required=True, type=click.Path(exists=True),help="Path to the exons consensus panel file")
@click.option("--depths-file", required=True, type=click.Path(exists=True),help="Path to the depths file")
@click.option("--omega-mutability-file", required=True, type=click.Path(exists=True),help="Path to the omega mutability file")
@click.option("--relative-mutability-file", required=True, type=click.Path(exists=True),help="Path to the relative mutability file")
@click.option("--resolution", type=click.Choice(['genomic', 'residue', 'genomic,residue']), show_default=True, default='genomic,residue', help="either genomic or residue based sites")
@click.option("--group-name", type=str, default="all_samples", show_default=True, help="Name of the group/sample to be used in the code")
def cli(somatic_mutations_file, vep_file, consensus_panel_file,
        omega_mutability_file, relative_mutability_file,
        depths_file, resolution, group_name,
        ):
    subsampling_rates = np.logspace(-2, np.log10(0.9), num=20)

    click.echo(f"Analyzing {group_name}")

    mutations = load_mutations(somatic_mutations_file)
    click.echo("Mutations loaded")

    vep = collect_vep(vep_file)
    click.echo("VEP data collected")

    df_panel_orig = load_panel(consensus_panel_file, depths_file, group_name, vep)
    click.echo(f"Panel loaded with {df_panel_orig.shape[0]} sites")

    curves_folder = f'{group_name}.curves'
    os.makedirs(curves_folder, exist_ok=True)

    # df_panel represents the total number of mutable sites,
    # either genomic or residue sites
    for impact in CONSEQUENCES_CATEGORIES.keys():
        df_panel = df_panel_orig[df_panel_orig["IMPACT"].isin(CONSEQUENCES_CATEGORIES[impact])]
    
        click.echo(f"Panel filtered with {df_panel.shape[0]} sites for {impact} mutations")
        mutations_lite = pd.merge(mutations, 
                            df_panel[['CHROM', 'POS', 'REF', 'ALT', 'DEPTH', 'GENE', 'AACHANGE', 'RESIDUE']],
                            on=['CHROM', 'POS', 'REF', 'ALT'], 
                            how='left')
        mutations_lite.dropna(axis=0, inplace=True)  # drop nan sites

        if 'residue' in resolution:
            mutations_residue = compute_mutation_rates(mutations_lite, group_name, impact, subsampling_rates, residue=True)
        if 'genomic' in resolution:
            mutations_genomic = compute_mutation_rates(mutations_lite, group_name, impact, subsampling_rates, residue=False)

        click.echo("Mutation rates computed")

        # load mutations
        mutations_dict = {
            'genomic': mutations_genomic if 'genomic' in resolution else None,
            'residue': mutations_residue if 'residue' in resolution else None
        }

        df_panel_genomic = df_panel.groupby(['POS', 'GENE']).agg({'DEPTH': 'mean'}).reset_index()
        df_panel_residue = df_panel.groupby(['RESIDUE', 'GENE']).agg({'DEPTH': 'mean'}).reset_index()
        df_panel_dict = {
            'genomic': df_panel_genomic,
            'residue': df_panel_residue
            }


        if 'residue' in resolution:
            click.echo("Plotting empirical discovery for residue sites")
            with PdfPages(f'{curves_folder}/residue_{impact}_empirical.pdf') as empirical_pdf, \
                 PdfPages(f'{curves_folder}/residue_{impact}_theoretical_empirical.pdf') as combined_pdf, \
                 PdfPages(f'{curves_folder}/residue_{impact}_slopes.pdf') as slope_pdf:
                main_empirical(group_name, mutations_dict, df_panel, df_panel_dict,
                            omega_mutability_file, relative_mutability_file,
                            subsampling_rates,
                            sites='residue',
                            impact = impact,
                            logscale=False,
                            empirical_pdf=empirical_pdf,
                            combined_pdf=combined_pdf,
                            slope_pdf=slope_pdf,
                            # genes_list = ["TP53","RBM10"]
                            )

        if 'genomic' in resolution:
            click.echo("Plotting empirical discovery for genomic sites")
            with PdfPages(f'{curves_folder}/genomic_{impact}_empirical.pdf') as empirical_pdf, \
                 PdfPages(f'{curves_folder}/genomic_{impact}_theoretical_empirical.pdf') as combined_pdf, \
                 PdfPages(f'{curves_folder}/genomic_{impact}_slopes.pdf') as slope_pdf:
                main_empirical(group_name, mutations_dict, df_panel, df_panel_dict,
                            omega_mutability_file, relative_mutability_file,
                            subsampling_rates,
                            sites='genomic',
                            impact = impact,
                            logscale=False,
                            empirical_pdf=empirical_pdf,
                            combined_pdf=combined_pdf,
                            slope_pdf=slope_pdf,
                            # genes_list = ["TP53","RBM10"]
                            )

    # no per-gene files are written anymore; PDFs are closed by their context managers


if __name__ == '__main__':
    cli()