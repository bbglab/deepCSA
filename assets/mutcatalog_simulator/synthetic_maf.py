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

def sci_no_dots(x, sig=2):

    if x % 1 == 0:
        return f"{int(x)}"
    
    mantissa, exp = f"{x:.{sig}e}".split('e')
    exp = int(exp)

    decimals = len(mantissa) - mantissa.index('.') - 1
    mantissa = mantissa.replace('.', '')
    exp -= decimals

    return f"{mantissa}e{exp:+03d}"


def get_lambda(expected_alt_depth, impact, omega_dict, depth_correction):
    
    return depth_correction * expected_alt_depth * omega_dict.get(impact, 1)



def poisson_simulate(sample, sample_mutability, genes, omega_dict, depth_correction, possible_muts, random_replicates=100):

    possible_muts_gene = possible_muts[possible_muts['GENE'].isin(genes)].copy()

    simulated_maf = pd.merge(possible_muts_gene,
                             sample_mutability[['CHROM', 'POS', 'CONTEXT_MUT', 'expected_ALT_DEPTH']],
                             on=['CHROM', 'POS', 'CONTEXT_MUT'], how='left')

    # define the vector of lambdas, i.e., mean parameter for Poisson distribution

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

    cleaned_possible_muts_init = pd.read_table(f"{deepcsa_run_dir}/regions/consensuspanels/consensus.exons_splice_sites.tsv")

    # --- loop across grid params ---
    compiled_output_depths_table = []
    for sample, omega, depth_correction in product(samples, omegas, depth_corrections):
        
        # --- omega dict ---

        omega_dict = {'missense': omega, 'nonsense': omega}

        # --- mutabilities per site ---

        fn = f"{deepcsa_run_dir}/processing_files/absolutemutabilitiesgloballoc/mutabilities_per_site.{sample}.global_loc.tsv.gz"
        sample_mutability = pd.read_table(fn)
        sample_mutability.rename(columns={sample: 'expected_ALT_DEPTH'}, inplace=True)

        # --- DEPTH column ---
        fn = f"{deepcsa_run_dir}/depths/individual/{sample}.depths.annotated.tsv.gz"
        depths_per_position = pd.read_csv(fn, sep = "\t", usecols = ['CHROM', 'POS', sample])
        depths_per_position = depths_per_position.rename(columns={sample: 'DEPTH'})
        cleaned_possible_muts = cleaned_possible_muts_init.merge(depths_per_position, on=['CHROM', 'POS'], how ='left')
        cleaned_possible_muts["DEPTH"] = (cleaned_possible_muts["DEPTH"] * depth_correction // 1).astype(int)

        experiment_name = f"{sample}_{sci_no_dots(omega)}_{sci_no_dots(depth_correction)}"
        print(f"Simulating samples for {experiment_name}...")

        counts, simulated_maf = poisson_simulate(
            sample, 
            sample_mutability, 
            genes, omega_dict, depth_correction, 
            cleaned_possible_muts, 
            random_replicates=n_replicates
            )

        # Generate output depths table
        base = cleaned_possible_muts.loc[cleaned_possible_muts["GENE"].isin(genes),["CHROM", "POS", "DEPTH"]].drop_duplicates().reset_index(drop=True)
        rep_cols = [ f"{experiment_name}_{i:02d}" for i in range(n_replicates) ]
        output_depths_table = base.assign( **{c: base["DEPTH"] for c in rep_cols} ).drop(columns="DEPTH")
        output_depths_table.to_csv( f"{output_dir}/{experiment_name}_depths.tsv", sep="\t", header = True,index=False)
        compiled_output_depths_table.append(output_depths_table.iloc[:,2:])

        # Generate output MAFs
        for i, s in tqdm.tqdm(enumerate(counts)):
            simulated_maf['ALT_DEPTH'] = s
            df = simulated_maf[simulated_maf['ALT_DEPTH'] > 0].copy()
            generated_sample_name = f"{experiment_name}_{i:02d}"
            df['SAMPLE_ID'] = generated_sample_name
            df.to_csv(
                f"{output_dir}/{generated_sample_name}.tsv",
                sep = "\t",
                header = True,
                index = False
                )

    # --- convenient printouts ---
    pd.concat([output_depths_table.iloc[:,:2]] + compiled_output_depths_table,
              axis=1).to_csv( f"{output_dir}/all_experiments_depths.tsv.gz",
                             sep="\t", header = True, index=False)
    print(f'output_dir: {output_dir}')

if __name__ == '__main__':
    
    synthetic_maf()
