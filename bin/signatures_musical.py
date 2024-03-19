#!/usr/local/bin/python


# Import some necessary modules
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import pandas as pd
import time
import scipy as sp
import pickle
import click
from matplotlib.backends.backend_pdf import PdfPages


# Import MuSiCal
import musical


def run_musical(matrix_file, cpus, minprocess, maxprocess):
    """
    INFO
    """

    X = pd.read_table(matrix_file, index_col=0)
    X["change"] = [ x[1:6] for x in X.index ]
    X = X.sort_values(by = ["change", "CONTEXT_MUT"]).drop("change", axis = 'columns')
    print(X.head())

    model = musical.DenovoSig(X,
                                min_n_components=minprocess,     # Minimum number of signatures to test
                                max_n_components=maxprocess,     # Maximum number of signatures to test
                                init='random',          # Initialization method
                                method='mvnmf',         # mvnmf or nmf
                                n_replicates=20,        # Number of mvnmf/nmf replicates to run per n_components
                                ncpu=cpus,              # Number of CPUs to use
                                max_iter=100000,        # Maximum number of iterations for each mvnmf/nmf run
                                bootstrap=True,         # Whether or not to bootstrap X for each run
                                tol=1e-8,               # Tolerance for claiming convergence of mvnmf/nmf
                                verbose=1,              # Verbosity of output
                                normalize_X=False       # Whether or not to L1 normalize each sample in X before mvnmf/nmf
                                )
    model.fit()

    # Number of discovered de novo signatures
    print(model.n_components)

    with PdfPages(f'test.pdf') as pdf:
        model.plot_selection()
        pdf.savefig()
        plt.close()

        fig = musical.sigplot_bar(model.W)
        pdf.savefig()
        plt.close()

    return model

def matching_n_refitting(model):
    thresh_grid = np.array([
        # 0.0001, 0.0002, 0.0005,
        0.001, 0.002, 0.005,
        0.01, 0.02, 0.05,
        # 0.1, 0.2, 0.5,
        # 1., 2., 5.
    ])


    catalog = musical.load_catalog('COSMIC-MuSiCal_v3p2_SBS_WGS')
    W_catalog = catalog.W
    print(W_catalog.shape[1])

    model.assign_grid(W_catalog,
                    method_assign='likelihood_bidirectional',   # Method for performing matching and refitting
                    thresh_match_grid=thresh_grid,              # Grid of threshold for matchinng
                    thresh_refit_grid=thresh_grid,              # Grid of threshold for refitting
                    thresh_new_sig=0.0,                         # De novo signatures with reconstructed cosine similarity below this threshold will be considered novel
                    connected_sigs=False,                       # Whether or not to force connected signatures to co-occur
                    clean_W_s=False                             # An optional intermediate step to avoid overfitting to small backgrounds in de novo signatures for 96-channel SBS signatures
                    )

    model.validate_grid(validate_n_replicates=1, # Number of simulation replicates to perform for each grid point
                        grid_selection_method='pvalue', # Method for selecting the best grid point
                        grid_selection_pvalue_thresh=0.05 # Threshold used for selecting the best grid point
                    )

    print(model.best_grid_point)
    print(model.thresh_match)
    print(model.thresh_refit)

    W_s = model.W_s
    H_s = model.H_s

    W_s.to_csv("test.W_s.tsv", header = True, index = True, sep = '\t')
    H_s.to_csv("test.H_s.tsv", header = True, index = True, sep = '\t')


@click.command()
@click.option('--matrixfile', type=click.Path(exists=True), help='File listing decomposed mutation probability files.')
@click.option('--cpus', type=int, help='Number of CPUs to use.')
@click.option('--minprocess', type=int, help='Minimum number of processes to test.')
@click.option('--maxprocess', type=int, help='Maximum number of processes to test.')

def main(matrixfile, cpus, minprocess, maxprocess):
    click.echo(f"Computing signatures with MuSiCal...")
    modd = run_musical(matrixfile, cpus, minprocess, maxprocess)
    matching_n_refitting(modd)

if __name__ == '__main__':
    main()

