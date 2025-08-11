#!/usr/bin/env python

import click
import pandas as pd
import numpy as np


# Consequence type ontology

broadimpact_grouping_dict = { 
    "missense": ["missense"], 
    "nonsense": ["nonsense"], 
    "essential_splice": ["essential_splice"], 
    "splice_region_variant": ["splice_region_variant"], 
    "truncating": ["nonsense", "essential_splice"], 
    "essential_splice_plus": ["essential_splice", "splice_region_variant"], 
    "truncating_plus": ["nonsense", "essential_splice", "splice_region_variant"], 
    "nonsynonymous_splice": ["missense", "nonsense", "essential_splice"],
    "synonymous": ["synonymous"]
    } 


# Utils

def triplet_context_iterator():
    subs = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    flanks = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
              'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    for sub in subs:
        for flank in flanks:
            yield flank[0] + sub[0] + flank[1] + '>' + sub[-1]


# define canonical order for trinucleotide contexts

triplet_contexts = list(triplet_context_iterator())

# The following function computes the normalization factor so that the mutation density
# matches the equivalent to 1 mut/(Mb-genome) on average across the mappable genome.

def get_correction_factor(sample_name, trinucleotide_counts_df, mutability_df, flat=False):

    # vector of triplet count in 96-channel canonical sorting

    triplet_counts_dict = dict(zip(trinucleotide_counts_df['CONTEXT'].values, trinucleotide_counts_df['COUNT'].values))
    l = [triplet_counts_dict[c[:3]] for c in triplet_contexts]
    triplet_counts = np.array(l)

    # genome length in Mb

    genome_length = sum(triplet_counts) / (3 * 1e6)

    # vector of relative mutabilities in 96-channel canonical sorting
    
    if flat:
        relative_mutability = (1 / 96) * np.ones(96)
    else:
        d = dict(zip(mutability_df['CONTEXT_MUT'], mutability_df[f'{sample_name}.all']))
        relative_mutability = np.array([d.get(k, 0) for k in triplet_contexts])

    # return correction factor a.k.a. "alpha hat"

    correction_factor = genome_length / np.dot(relative_mutability, triplet_counts)
    
    return correction_factor, relative_mutability


def mutation_density(sample_name, depths_file, somatic_mutations_file, mutability_file, panel_file, trinucleotide_counts_file):

    depths_df = pd.read_csv(depths_file, sep='\t')
    somatic_mutations_df = pd.read_csv(somatic_mutations_file, sep='\t')
    mutability_df = pd.read_csv(mutability_file, sep='\t')
    panel_df = pd.read_csv(panel_file, sep='\t')
    trinucleotide_counts_df = pd.read_csv(trinucleotide_counts_file, sep='\t')

    genes = panel_df['GENE'].unique()
    consequences = broadimpact_grouping_dict.keys()

    res = pd.DataFrame(index=genes, columns=consequences)

    for csqn, csqn_set in broadimpact_grouping_dict.items():
        
        for gene in panel_df['GENE'].unique():
            
            # compute vector of sum of depths per trinucleotide context

            region_df = panel_df[(panel_df['IMPACT'].isin(csqn_set)) & (panel_df['GENE'] == gene)].copy()
            dh = pd.merge(region_df[['CHROM', 'POS']], depths_df[['CHROM', 'POS', 'CONTEXT', sample_name]], on=['CHROM', 'POS'], how='left')
            depth_sum_df = dh.groupby(by='CONTEXT').agg({sample_name: 'sum'}).reset_index()
            depth_region_dict = dict(zip(depth_sum_df['CONTEXT'], depth_sum_df[sample_name]))
            depth_region_vect = np.array([depth_region_dict.get(k[:3], 0) for k in triplet_contexts])

            # compute correction factor "alpha hat"

            correction_factor, relative_mutability = get_correction_factor(sample_name, trinucleotide_counts_df, mutability_df)

            # compute effective length

            effective_length = correction_factor * np.dot(relative_mutability, depth_region_vect)
            
            try:
                assert(effective_length != 0)
            except AssertionError:
                res.loc[gene, csqn] = None
                continue
            
            # observed somatic mutations

            n = somatic_mutations_df[(somatic_mutations_df['canonical_Consequence_broader'].isin(csqn_set)) & (somatic_mutations_df['SYMBOL'] == gene) & ((somatic_mutations_df['TYPE'] == 'SNV'))].shape[0]

            res.loc[gene, csqn] = n / effective_length
    
    return res



@click.command()
@click.option('--sample_name', type=str, help='Name of the sample being processed.')
@click.option('--depths_file', type=str, help='Depth per position in panel')
@click.option('--somatic_mutations_file', type=click.Path(exists=True), help='Observed somatic mutations')
@click.option('--mutability_file', type=click.Path(exists=True), help='Relative mutability per trinucleotide context')
@click.option('--panel_file', type=click.Path(exists=True), help='All possible SNVs with csqn type')
@click.option('--trinucleotide_counts_file', type=click.Path(exists=True), help='How many trinucleotides of each type are in the entire genome')
def main(sample_name, depths_file, somatic_mutations_file, mutability_file, panel_file, trinucleotide_counts_file):

    click.echo(f"Running the mutability adjusted mutation density estimation...")
    res = mutation_density(sample_name, depths_file, somatic_mutations_file, mutability_file, panel_file, trinucleotide_counts_file)
    res.to_csv(f'{sample_name}.mutdensities.tsv', sep='\t')


if __name__ == '__main__':
    main()