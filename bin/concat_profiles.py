#!/usr/bin/env python


import click
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import cosine


def compile_profiles(mutation_profile_files):
    """
    Processes mutation probability signature files and outputs per-sample probability files.
    Plots both a cosine similarity heatmap and a clustermap.
    """
    # Read all profile files listed in the input file
    mut_profile_matrix = pd.concat(
        (pd.read_csv(file.strip(), sep='\t', header=0).set_index("CONTEXT_MUT") for file in open(mutation_profile_files, 'r')),
        axis=0
    )
    mut_profile_matrix = mut_profile_matrix.fillna(0)
    cols_list = list(mut_profile_matrix.columns)
    mode = cols_list[0].split(".")[-1]
    mut_profile_matrix.columns = [ '.'.join(x.split(".")[:-1]) for x in cols_list ]
    mut_profile_matrix.to_csv(f"{mode}.compiled_profiles.tsv", header=True, index=True, sep="\t")

    all_samples_names = mut_profile_matrix.columns
    # Compute cosine similarity matrix
    cosine_similarity_df = pd.DataFrame(index=all_samples_names, columns=all_samples_names, dtype=float)
    for i, sample in enumerate(all_samples_names):
        for j, sample2 in enumerate(all_samples_names):
            if j < i:
                # Copy the already computed value (matrix is symmetric)
                cosine_similarity_df.loc[sample, sample2] = cosine_similarity_df.loc[sample2, sample]
            else:
                similarity = 1 - cosine(mut_profile_matrix[sample], mut_profile_matrix[sample2])
                cosine_similarity_df.loc[sample, sample2] = similarity
    cosine_similarity_df.to_csv(f"{mode}.cosine_similarity.tsv", header=True, index=True, sep="\t")

    # Plot heatmap
    plt.figure(figsize=(12, 8))
    ax = sns.heatmap(cosine_similarity_df, 
                annot=True,
                fmt=".2f",
                cmap='coolwarm',
                # vmin=0, vmax=1,
                cbar_kws={'label': 'Cosine Similarity'})
    plt.title('Cosine Similarity Heatmap')
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(f"{mode}.cosine_similarity_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()


    # Plot clustermap
    clustergrid = sns.clustermap(cosine_similarity_df,
                   annot=True,
                   fmt=".2f",
                   cmap='coolwarm',
                #    vmin=0, vmax=1,
                   cbar_kws={'label': 'Cosine Similarity'},
                   figsize=(12, 8))
    plt.suptitle('Cosine Similarity Clustermap')
    clustergrid.savefig(f"{mode}.cosine_similarity_clustermap.png", dpi=300, bbox_inches='tight')
    plt.close()

@click.command()
@click.option('--mutation-profiles-list', type=click.Path(exists=True), help='File listing mutational profile files.')

def main(mutation_profiles_list):
    click.echo("Combining signature probabilities of all samples and groups...")
    compile_profiles(mutation_profiles_list)

if __name__ == '__main__':
    main()
