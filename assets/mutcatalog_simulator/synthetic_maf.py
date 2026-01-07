import os
import json
from itertools import product

import tqdm
import click
import pandas as pd
import numpy as np

# Usage:
# -----
# python synthetic_maf.py --deepcsa_run_dir (path) -- run_config (json) --output_dir (path) 

# Arguments:
# ---------
# * deepcsa_run_dir (path)
# * run_config (json, dict)
#    - samples  
#    - omegas
#    - genes
#    - depth_corrections
#    - random_replicates
# * output_dir (path)


def get_lambda(expected_alt_depth, impact, omega_dict, depth_correction):
    
    return depth_correction * expected_alt_depth * omega_dict.get(impact, 1)



def poisson_simulate(sample, sample_mutability, gene, omega_dict, depth_correction, possible_muts, random_replicates=100):

    possible_muts_gene = possible_muts[possible_muts['GENE'] == gene].copy()

    simulated_maf = pd.merge(possible_muts_gene,
                             sample_mutability[['CHROM', 'POS', 'CONTEXT_MUT', 'expected_ALT_DEPTH']],
                             on=['CHROM', 'POS', 'CONTEXT_MUT'], how='left')

    # define the vector of lambdas, i.e., mean parameter for Poisson distribution

    simulated_maf['SAMPLE_ID'] = sample
    simulated_maf['lambda'] = simulated_maf.apply(lambda x: get_lambda(x['expected_ALT_DEPTH'], x['IMPACT'], omega_dict, depth_correction), axis=1)

    rng = np.random.default_rng(seed=12345)
    random_counts = rng.poisson(lam=simulated_maf['lambda'].values, size=(random_replicates, simulated_maf.shape[0]))

    return random_counts, simulated_maf



@click.command()
@click.option('--deepcsa_run_dir')
@click.option('--run_config')
@click.option('--output_dir')
def synthetic_maf(deepcsa_run_dir, run_config, output_dir):

    # --- load grid params ---

    with open(run_config, 'rt') as f:
        config = json.load(f)
    
    samples = config['samples']
    omegas = config['omegas']
    genes = config['genes']
    depth_corrections = config['depth_corrections']
    n_replicates = config['n_replicates']
    
    # --- prepare all possible mutation sites ---

    fn = f"{deepcsa_run_dir}/createpanels/capturedpanels/captured_panel.exons_splice_sites.tsv"
    possible_muts = pd.read_table(fn)

    # --- protein mapping filtering ---
    
    all_protein_positions_file = f"{deepcsa_run_dir}/dna2proteinmapping/depths_per_position_exon_gene.tsv"
    all_protein_positions = pd.read_table(all_protein_positions_file)
    all_protein_positions = all_protein_positions[all_protein_positions["COVERED"] == 1].reset_index(drop = True)
    all_protein_positions = all_protein_positions.rename({"CHR": "CHROM", "DNA_POS": "POS"}, axis = 'columns')
    cleaned_possible_muts = possible_muts.merge(all_protein_positions, on = ['CHROM', 'POS', 'GENE', 'CONTEXT'], how = 'left')
    cleaned_possible_muts = cleaned_possible_muts[~(cleaned_possible_muts["PROT_POS"].isna())].reset_index(drop = True)

    # --- loop across grid params ---

    for sample, gene, omega, depth_correction in product(samples, genes, omegas, depth_corrections):
        
        # --- omega dict ---

        omega_dict = {'missense': omega, 'nonsense': omega}

        # --- mutabilities per site ---

        fn = f"{deepcsa_run_dir}/absolutemutabilitiesgloballoc/mutabilities_per_site.{sample}.global_loc.tsv.gz"
        sample_mutability = pd.read_table(fn)
        sample_mutability.rename(columns={sample: 'expected_ALT_DEPTH'}, inplace=True)

        # --- DEPTH column ---

        cleaned_possible_muts = cleaned_possible_muts.drop(columns=['DEPTH'])
        fn = f"{deepcsa_run_dir}/annotatedepths/{sample}.depths.annotated.tsv.gz"
        depths_per_position = pd.read_csv(fn, sep = "\t", usecols = ['CHROM', 'POS', sample])
        depths_per_position = depths_per_position.rename(columns={sample: 'DEPTH'})
        cleaned_possible_muts = cleaned_possible_muts.merge(depths_per_position, on=['CHROM', 'POS'], how ='left')

        counts, simulated_maf = poisson_simulate(
            sample, 
            sample_mutability, 
            gene, omega_dict, depth_correction, 
            cleaned_possible_muts, 
            random_replicates=n_replicates
            )

        for i, s in tqdm.tqdm(enumerate(counts)):
            simulated_maf['ALT_DEPTH'] = s
            df = simulated_maf[simulated_maf['ALT_DEPTH'] > 0]
            df.to_csv(
                f"{output_dir}/{sample}.{gene}.{omega}.{depth_correction}.fake_mutations.{i}.tsv",
                sep = "\t",
                header = True,
                index = False
                )

    # --- convenient printouts ---

    print(f'output_dir: {output_dir}')

if __name__ == '__main__':
    
    synthetic_maf()
