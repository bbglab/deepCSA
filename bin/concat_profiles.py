#!/usr/bin/env python



import click
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import cosine
import json


def plot_similarity_heatmaps(mut_profile_matrix, mode, keys, label):
    # Only keep columns matching keys
    selected_cols = [col for col in keys if col in mut_profile_matrix.columns]
    if not selected_cols:
        print(f"No columns found for {label}")
        return
    sub_matrix = mut_profile_matrix[selected_cols]
    # Compute cosine similarity matrix
    cosine_similarity_df = pd.DataFrame(index=selected_cols, columns=selected_cols, dtype=float)
    for i, sample in enumerate(selected_cols):
        for j, sample2 in enumerate(selected_cols):
            if j < i:
                cosine_similarity_df.loc[sample, sample2] = cosine_similarity_df.loc[sample2, sample]
            else:
                similarity = 1 - cosine(sub_matrix[sample], sub_matrix[sample2])
                cosine_similarity_df.loc[sample, sample2] = similarity

    if label == 'all':
        cosine_similarity_df.to_csv(f"{mode}.cosine_similarity_{label}.tsv", header=True, index=True, sep="\t")

    # Plot heatmap
    plt.figure(figsize=(max(8, 0.5*len(selected_cols)), max(6, 0.4*len(selected_cols))))
    ax = sns.heatmap(cosine_similarity_df, 
                annot=True,
                fmt=".2f",
                cmap='coolwarm',
                cbar_kws={'label': 'Cosine Similarity'},
                annot_kws={"fontsize":6})
    plt.title(f'Cosine Similarity Heatmap ({label})')
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(f"{mode}.cosine_similarity_heatmap_{label}.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Plot clustermap
    clustergrid = sns.clustermap(cosine_similarity_df,
                    annot=True,
                    fmt=".2f",
                    cmap='coolwarm',
                    cbar_kws={'label': 'Cosine Similarity'},
                    figsize=(max(8, 0.5*len(selected_cols)), max(6, 0.4*len(selected_cols))),
                    annot_kws={"fontsize":6})
    plt.suptitle(f'Cosine Similarity Clustermap ({label})')
    clustergrid.savefig(f"{mode}.cosine_similarity_clustermap_{label}.png", dpi=300, bbox_inches='tight')
    plt.close()


def compile_profiles(mutation_profile_files, groups_json):
    """
    Processes mutation probability signature files and outputs per-sample probability files.
    Plots cosine similarity heatmaps and clustermaps for three sets:
    - keys with value size 1
    - keys with value size > 1
    - all keys
    """
    # Read all profile files listed in the input file
    mut_profile_matrix = pd.concat(
        (pd.read_csv(file.strip(), sep='\t', header=0).set_index("CONTEXT_MUT") for file in open(mutation_profile_files, 'r')),
        axis=1
    )
    mut_profile_matrix = mut_profile_matrix.fillna(0)
    cols_list = list(mut_profile_matrix.columns)
    mode = cols_list[0].split(".")[-1]
    mut_profile_matrix.columns = [ '.'.join(x.split(".")[:-1]) for x in cols_list ]
    mut_profile_matrix.to_csv(f"{mode}.compiled_profiles.tsv", header=True, index=True, sep="\t")

    # Load groups JSON
    with open(groups_json, 'r') as f:
        groups = json.load(f)

    # Partition keys
    keys_size1 = sorted([k for k, v in groups.items() if isinstance(v, list) and len(v) == 1 and v[0] == k])
    keys_sizegt1 = [k for k, v in groups.items() if isinstance(v, list) and len(v) > 1]

    # Sort each list by: (1) increasing number of '_' in the key, then (2) alphabetical order
    sort_key = lambda s: (s.count('_'), s)
    keys_sizegt1 = sorted(keys_sizegt1, key=sort_key)

    all_keys = list(keys_size1) + keys_sizegt1

    plot_similarity_heatmaps(mut_profile_matrix, mode, keys_size1, "size1")
    plot_similarity_heatmaps(mut_profile_matrix, mode, keys_sizegt1, "sizegt1")
    plot_similarity_heatmaps(mut_profile_matrix, mode, all_keys, "all")


@click.command()
@click.option('--mutation-profiles-list', type=click.Path(exists=True), help='File listing mutational profile files.')
@click.option('--groups-json', type=click.Path(exists=True), help='JSON file defining groups.')
def main(mutation_profiles_list, groups_json):
    click.echo("Combining signature probabilities of all samples and groups...")
    compile_profiles(mutation_profiles_list, groups_json)

if __name__ == '__main__':
    main()
