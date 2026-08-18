#!/usr/bin/env python

import os
import click
# from joblib import Parallel, delayed

import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

from parse_deepcsa import deepCSAparser
from dev_omega import Background, Omega


# all fitting methods
fitting_methods = ['fit_dndscv', 'fit_map', 'fit_posterior_marginal', 'fit_loc']
alias = dict(zip(fitting_methods, ['dndscv', 'map', 'posmar', 'loc']))

# cols for output dataframe
keys = ['gene_id', 'sample_id', 'n_syn', 'n_mis', 'n_trunc', 'offset_syn', 'offset_mis', 'offset_trunc', 'mean_prior', 'alpha', 'beta']
output_cols = ['omega_mis', 'low_mis', 'high_mis', 'pval_mis', 'pneg_mis', 'pcum_mis', 
                'omega_trunc', 'low_trunc', 'high_trunc', 'pval_trunc', 'pneg_trunc', 'pcum_trunc', 'syn_density', 'syn_density_low', 'syn_density_high']
for name in fitting_methods:
    keys += [col + f'_{alias[name]}' for col in output_cols]


def process_row(row, grouped=False):

    # initialize output dict

    res = {k: np.nan for k in keys}

    gene = row.gene_id
    alpha = row.alpha
    beta = row.beta
    n_syn = row.n_syn
    n_mis = row.n_mis
    n_trunc = row.n_trunc
    offset_syn = row.offset_syn
    offset_mis = row.offset_mis
    offset_trunc = row.offset_trunc
    mu = row.mean
    
    if not grouped:
        sample = row.sample_id
    else:
        del res['sample_id']

    # instantiate omega class

    omega = Omega(n_syn, n_mis, n_trunc, offset_syn, offset_mis, offset_trunc, alpha, beta, mu)

    res['n_syn'] = n_syn
    res['n_mis'] = n_mis
    res['n_trunc'] = n_trunc
    res['offset_syn'] = offset_syn
    res['offset_mis'] = offset_mis
    res['offset_trunc'] = offset_trunc
    res['gene_id'] = gene
    res['mean_prior'] = mu
    res['alpha'] = alpha
    res['beta'] = beta

    if not grouped:
        res['sample_id'] = sample
    
    # run different fitting approaches (omega)

    for name in fitting_methods:
        method = getattr(omega, name)
        omega_mis, lower_mis, upper_mis, pval_mis, pneg_mis, pcum_mis, omega_trunc, lower_trunc, upper_trunc, pval_trunc, pneg_trunc, pcum_trunc, t_hat, t_hat_low, t_hat_high = method()

        res[f'omega_mis_{alias[name]}'] = omega_mis
        res[f'low_mis_{alias[name]}'] = lower_mis
        res[f'high_mis_{alias[name]}'] = upper_mis
        res[f'pval_mis_{alias[name]}'] = pval_mis
        res[f'pneg_mis_{alias[name]}'] = pneg_mis
        res[f'pcum_mis_{alias[name]}'] = pcum_mis
        
        res[f'omega_trunc_{alias[name]}'] = omega_trunc
        res[f'low_trunc_{alias[name]}'] = lower_trunc
        res[f'high_trunc_{alias[name]}'] = upper_trunc
        res[f'pval_trunc_{alias[name]}'] = pval_trunc
        res[f'pneg_trunc_{alias[name]}'] = pneg_trunc
        res[f'pcum_trunc_{alias[name]}'] = pcum_trunc

        res[f'syn_density_{alias[name]}'] = t_hat
        res[f'syn_density_low_{alias[name]}'] = t_hat_low
        res[f'syn_density_high_{alias[name]}'] = t_hat_high
    
    return res


def fit_omega(data, grouped=False, test=False):

    if test:
        limit = 220
    else:
        limit = data.shape[0]

    # results = Parallel(n_jobs=-1, verbose=10)(
    #     delayed(process_row)(row, grouped=grouped) for row in data[: limit+1].itertuples(index=False)
    #     )
    results = [
        process_row(row, grouped=grouped) 
        for row in data[: limit + 1].itertuples(index=False)
    ]

    res = pd.DataFrame(results)
    
    # FDR

    for k, v in alias.items():
        for impact in ['mis', 'trunc']:
            colname = f'pval_{impact}_{v}'
            mask = res[colname].notna()
            res.loc[mask, f'qval_{impact}_{v}'] = multipletests(res.loc[mask, colname], method="fdr_bh")[1]

    return res


@click.command()
@click.option('--config', type=click.Path(), help="configuration file with essential input paths and filters", required = True)
@click.option('--outfolder', type=click.Path(), help="output folder")
@click.option('--roadmap-covariates', is_flag=True, help="flag: use the first 20 PCs from the Roadmap of Epigenomics as covariates in the background model")
@click.option('--samples-as-covariates', is_flag=True, help="flag: use samples as covariates in the background model")
@click.option('--all-samples', is_flag=True, help="flag: conduct only the analysis grouping all samples")
@click.option('--verb', is_flag=True, help="flag: prints out summary output of background model")
@click.option('--test', is_flag=True, help="flag: run in test mode: limited gene-samples, loading cached data")
def cli(config, outfolder, roadmap_covariates, samples_as_covariates, all_samples, verb, test):

    data_fn = os.path.join(outfolder, 'data.tsv')

    if not test:

        # parse granular data
        print('Parsing data...')

        parser = deepCSAparser(config)
        data = parser.parse_dataset()
        data.to_csv(data_fn, sep='\t', index=False)

    else:

        print('Loading cached data...')
        
        data = pd.read_csv(data_fn, sep='\t')

    # run negative-binomial regression with granular data

    print('Fitting background model...')

    background = Background(data)
    augmented_data, overdispersion, gene_coefficients, sample_coefficients = background.fit(
        verb=False, 
        covariates=roadmap_covariates, 
        samples=samples_as_covariates
    )
    
    # remove NaNs before proceeding to next step
    
    augmented_data.dropna(inplace=True)

    if not all_samples:
    
        # run omega with granular data
        
        print('Fitting dN/dS models...')

        output = fit_omega(augmented_data, test=test)

        output.to_csv(os.path.join(outfolder, 'omega.tsv'), sep='\t', index=False)

    # grouped data

    grouped = augmented_data.groupby('gene_id').agg(
        {
        **{
            'offset_syn': np.sum,
            'offset_trunc': np.sum,
            'offset_mis': np.sum,
            'n_syn': np.sum,
            'n_mis': np.sum,
            'n_trunc': np.sum
        },
        **{f'PC{i}': 'first' for i in range(1,21)}
        }
    ).reset_index()

    # run negative-binomial regression with grouped data

    print('Fitting background model on grouped data...')
    
    background = Background(grouped)
    augmented_data, overdispersion, gene_coefficients, sample_coefficients = background.fit(
        verb=False, 
        covariates=roadmap_covariates, 
        samples=False)

    # remove NaNs before proceeding to next step
    
    augmented_data.dropna(inplace=True)
    
    # run omega with grouped data

    print('Fitting dN/dS models on grouped data...')

    output = fit_omega(augmented_data, grouped=True, test=test)

    output.to_csv(os.path.join(outfolder, 'omega.grouped.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    cli()