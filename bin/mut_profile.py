#!/usr/bin/env python


import sys
import click
import pandas as pd
import numpy as np

from utils import contexts_formatted, contexts_formatted_sigprofiler
from utils_plot import plot_profile
from read_utils import custom_na_values

def bayesian_update(observed_counts, prior_profile, limit=200):

    N = np.sum(observed_counts).item()

    if N >= limit:
        return observed_counts / N
    else:
        prior_profile = prior_profile.set_index("CONTEXT_MUT")
        posterior_profile = ((limit - N) * prior_profile).values + observed_counts
        posterior_profile /= np.sum(posterior_profile)
        return posterior_profile


def compute_mutation_matrix(sample_name, mutations_file, mutation_matrix, method, pseudocount,
                            sigprofiler, per_sample):
    """
    Compute mutational profile from the input data
          ***Remember to add some pseudocounts to the computation***

        Required information:
            Annotated mutations observed
        Output:
            Mutational profile per sample, possibility to add pseudocounts to prevent some probabilities from being 0
    """

    if method not in ['unique', 'multiple']:
        print('a')
        exit(1)

    # Load your MAF DataFrame (raw_annotated_maf)
    annotated_maf = pd.read_csv(mutations_file, sep = "\t", header = 0, na_values = custom_na_values)

    # create the matrix in the desired order
    empty_matrix = pd.DataFrame(index = contexts_formatted)

    if method == 'unique':
        # make sure to count each mutation only once (avoid annotation issues)
        annotated_minimal_maf = annotated_maf[["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID"]].drop_duplicates().reset_index(drop = True)
        annotated_minimal_maf.columns = ["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID"]

        # count the mutations per sample and per context
        counts_x_sample_context_long = annotated_minimal_maf.groupby(by = ["SAMPLE_ID", "CONTEXT_MUT"])["MUT_ID"].count().reset_index()
        counts_x_sample_context_long.columns = ["SAMPLE_ID", "CONTEXT_MUT", "COUNT"]

    elif method == 'multiple':
        # make sure to count each mutation only once (avoid annotation issues)
        annotated_minimal_maf = annotated_maf[["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID", "ALT_DEPTH"]].drop_duplicates().reset_index(drop = True)
        annotated_minimal_maf.columns = ["SAMPLE_ID", "CONTEXT_MUT", "MUT_ID", "ALT_DEPTH"]

        # count the mutations per sample and per context
        counts_x_sample_context_long = annotated_minimal_maf.groupby(by = ["SAMPLE_ID", "CONTEXT_MUT"])["ALT_DEPTH"].sum().reset_index()
        counts_x_sample_context_long.columns = ["SAMPLE_ID", "CONTEXT_MUT", "COUNT"]

    # here we group the counts of all the samples
    counts_x_sample_matrix = counts_x_sample_context_long.groupby(by = ["CONTEXT_MUT"])["COUNT"].sum().reset_index()
    counts_x_sample_matrix.columns = ["CONTEXT_MUT", sample_name]
    counts_x_sample_matrix = counts_x_sample_matrix.set_index("CONTEXT_MUT")

    counts_x_sample_matrix = pd.concat( (empty_matrix, counts_x_sample_matrix) , axis = 1)
    counts_x_sample_matrix = counts_x_sample_matrix.fillna(0)
    counts_x_sample_matrix.index.name = "CONTEXT_MUT"

    # add a pseudocount if desired
    counts_x_sample_matrix = counts_x_sample_matrix + pseudocount

    counts_x_sample_matrix.to_csv(mutation_matrix,
                                    header = True,
                                    index = True,
                                    sep = "\t")

    if sigprofiler:
        counts_x_sample_matrix.index = contexts_formatted_sigprofiler
        counts_x_sample_matrix.index.name = "CONTEXT_MUT"
        counts_x_sample_matrix = counts_x_sample_matrix.reset_index().sort_values(by = "CONTEXT_MUT")
        counts_x_sample_matrix.to_csv(f"{mutation_matrix}.single.sigprofiler",
                                        header = True,
                                        index = False,
                                        sep = "\t")

    if per_sample:
        counts_x_sample_matrix = counts_x_sample_context_long.pivot(index = "CONTEXT_MUT", columns = "SAMPLE_ID", values = "COUNT")
        counts_x_sample_matrix = pd.concat( (empty_matrix, counts_x_sample_matrix) , axis = 1)
        counts_x_sample_matrix = counts_x_sample_matrix.fillna(0)
        counts_x_sample_matrix.index.name = "CONTEXT_MUT"

        counts_x_sample_matrix.to_csv(f"{mutation_matrix}.per_sample",
                                        header = True,
                                        index = True,
                                        sep = "\t")

        if sigprofiler:
            counts_x_sample_matrix.index = contexts_formatted_sigprofiler
            counts_x_sample_matrix.index.name = "CONTEXT_MUT"
            counts_x_sample_matrix = counts_x_sample_matrix.reset_index().sort_values(by = "CONTEXT_MUT")
            counts_x_sample_matrix.to_csv(f"{mutation_matrix}.per_sample.sigprofiler",
                                            header = True,
                                            index = False,
                                            sep = "\t")



def profile_stability(counts, denominator):
    """
    Compute stability score of a probability profile by
    adding +1 to each channel individually and measuring
    L1 (Total Variation) deviation.

    Parameters
    ----------
    counts : array-like, shape (k,)
        Raw integer counts per channel.
    denominator : array-like, shape (k,)
        Denominator for each channel (e.g., trinucleotide counts).

    Returns
    -------
    results : dict
        {
            'mean_deviation': float,
            'max_deviation': float,
            'std_deviation': float,
            'all_deviations': np.ndarray shape (k,)
        }
    """

    # Convert to numpy arrays
    counts = np.asarray(counts, dtype=float)
    denominator = np.asarray(denominator, dtype=float)

    if counts.ndim != 1:
        raise ValueError("counts must be a 1D vector")
    if denominator.ndim != 1:
        raise ValueError("denominator must be a 1D vector")
    if len(counts) != len(denominator):
        raise ValueError("counts and denominator must have same length")
    if np.any(denominator == 0):
        raise ValueError("denominator contains zeros")

    k = len(counts)

    # ---- Baseline profile ----
    ref_profile = counts / denominator
    ref_profile = ref_profile / ref_profile.sum()
    p = ref_profile

    deviations = np.zeros(k)

    # ---- Perturb each channel ----
    for j in range(k):
        perturbed = counts.copy()
        perturbed[j] += 1

        p_prime = perturbed / denominator
        p_prime = p_prime / p_prime.sum()  # IMPORTANT: renormalize

        deviations[j] = 0.5 * np.sum(np.abs(p - p_prime))

    return {
        "mean_deviation": deviations.mean(),
        "min_deviation": deviations.min(),
        "max_deviation": deviations.max(),
        "std_deviation": deviations.std(),
        "all_deviations": deviations.tolist()
    }


def compute_mutation_profile(sample_name, mutation_matrix_file, trinucleotide_counts_file, plot,
                                wgs = False, wgs_trinucleotide_counts = False, sigprofiler = False,
                                smoothed = False, prior_profile_file = None):
    """
    Compute mutational profile from the input data

        Required information:
            Mutation matrix
            Trinucleotide content of the sequenced region (depth-aware or non-depth-aware)
        Output:
            Mutational profile per sample
    """

    # Load your mutation matrix
    mutation_matrix = pd.read_csv(mutation_matrix_file, sep = "\t", header = 0)
    mutation_matrix = mutation_matrix.set_index("CONTEXT_MUT")
    total_mutations = np.sum(mutation_matrix[sample_name])

    if smoothed:
        prior_profile = pd.read_table(prior_profile_file)
        minimum_mutations = 200
        mutation_matrix = bayesian_update(mutation_matrix, prior_profile, minimum_mutations)
        total_mutations = max(minimum_mutations, total_mutations)

    # proportion of SBS mutations per trinucleotide in panel
    mutation_matrix_proportions = mutation_matrix.copy()
    mutation_matrix_proportions[sample_name] = mutation_matrix_proportions[sample_name] / total_mutations
    mutation_matrix_proportions.to_csv(f"{sample_name}.proportion_mutations.tsv",
                                header = True,
                                index = True,
                                sep = "\t")


    # Load the trinucleotide counts
    trinucleotide_counts = pd.read_csv(trinucleotide_counts_file, sep = "\t", header = 0)
    trinucleotide_counts = trinucleotide_counts.set_index("CONTEXT")

    trinuc_per_contextmuts = []
    for x in contexts_formatted:
        trinuc_per_contextmuts.append(trinucleotide_counts.loc[x[:3], sample_name])

    contextmut_depth = pd.DataFrame({sample_name : trinuc_per_contextmuts, "CONTEXT_MUT": contexts_formatted})
    contextmut_depth = contextmut_depth.set_index("CONTEXT_MUT")

    # compute profile stability and store the outputs in a file
    stability_output_dict = profile_stability(mutation_matrix[sample_name], contextmut_depth[sample_name])
    stability_outputs = pd.DataFrame(stability_output_dict.items()).set_index(0).T
    stability_outputs[["SAMPLE_ID", "mode"]] = sample_name.split(".")
    stability_outputs[["SAMPLE_ID", "mode"] + list(stability_outputs.columns[:-2])].to_csv(f"{sample_name}.profile_stability.tsv", sep = "\t", index = False)


    # divide
    mut_probability = mutation_matrix.divide( contextmut_depth )

    # normalize
    mut_probability = mut_probability / mut_probability.sum()

    # if there is any channel with 0 probability we need to add a pseudocount
    if not all(mut_probability[sample_name].values > 0):
        # find the minimum value greater than 0
        min_value_non_zero = mut_probability[mut_probability > 0].min()
        # print(min_value_non_zero)

        # add a dynamic pseudocount of one third of the minimum number greater than 0
        mut_probability = mut_probability + (min_value_non_zero / 3)

    mut_probability = mut_probability / mut_probability.sum()


    # reindex to ensure the right order
    mut_probability = mut_probability.reindex(contexts_formatted)
    mut_probability.index.name = "CONTEXT_MUT"
    mut_probability = mut_probability.reset_index()

    mut_probability.to_csv(f"{sample_name}.profile.tsv",
                            header = True,
                            index = False,
                            sep = "\t")

    if plot:

        # plot the profile as the relative probability of each trinucleotide change to mutate vs the others
        plot_profile(dict(zip(mut_probability["CONTEXT_MUT"], mut_probability[sample_name])),
                        title=f'{sample_name} ({round(total_mutations)} muts)',
                        yaxis_name= "Relative mutation\nprobability",
                        output_f = f'{sample_name}.profile.pdf')

        # plot the profile as a percentage of SBS mutations seen in our sequenced panel
        plot_profile(dict(zip(mutation_matrix_proportions.index, mutation_matrix_proportions[sample_name])),
                        title=f'{sample_name} ({round(total_mutations)} muts)',
                        output_f = f'{sample_name}.profile.percentage.pdf')


    if wgs:
        if not wgs_trinucleotide_counts:
            print("Invalid wgs_trinucleotide_counts, provide a correct file: ", wgs_trinucleotide_counts)
            sys.exit("Invalid wgs_trinucleotide_counts, provide a correct file.")

        mut_probability["CONTEXT"] = mut_probability["CONTEXT_MUT"].apply( lambda x : x[:3])
        ref_trinuc32 = pd.read_csv(wgs_trinucleotide_counts,
                                    sep = "\t", header = 0)

        profile_trinuc_merge = mut_probability.merge(ref_trinuc32, on = "CONTEXT")
        profile_trinuc_merge["MUTS_WGS"] = profile_trinuc_merge[sample_name] * profile_trinuc_merge["COUNT"]
        profile_trinuc_merge["SAMPLE_MUTS_WGS"] = profile_trinuc_merge["MUTS_WGS"] / np.sum(profile_trinuc_merge["MUTS_WGS"]) * total_mutations
        profile_trinuc_clean = profile_trinuc_merge[["CONTEXT_MUT", "SAMPLE_MUTS_WGS"]].set_index("CONTEXT_MUT")
        profile_trinuc_clean.index.name = "CONTEXT_MUT"
        profile_trinuc_clean = profile_trinuc_clean.reindex(contexts_formatted)
        profile_trinuc_clean.columns = [sample_name]

        profile_trinuc_clean.to_csv(f"{sample_name}.matrix.WGS.tsv",
                                    header = True,
                                    index = True,
                                    sep = "\t")

        profile_trinuc_clean_proportion = profile_trinuc_clean.copy()
        profile_trinuc_clean_proportion[sample_name] = profile_trinuc_clean_proportion[sample_name] / profile_trinuc_clean_proportion[sample_name].sum()
        profile_trinuc_clean_proportion.to_csv(f"{sample_name}.proportion_mutations.WGS.tsv",
                                                header = True,
                                                index = True,
                                                sep = "\t")


        # plot the profile as a percentage of SBS mutations seen after sequencing one WGS
        # if mutations were occuring with the same probabilities as they occur in our sequenced panel
        plot_profile(dict(zip(profile_trinuc_clean_proportion.index.values, profile_trinuc_clean_proportion[sample_name].values)),
                        title=f'{sample_name} ({round(total_mutations)} muts)',
                        output_f = f'{sample_name}.profile.percentage_WGS.pdf')

        if sigprofiler:
            profile_trinuc_clean.index = contexts_formatted_sigprofiler
            profile_trinuc_clean.index.name = "CONTEXT_MUT"
            profile_trinuc_clean = profile_trinuc_clean.reset_index().sort_values(by = "CONTEXT_MUT")
            profile_trinuc_clean.to_csv(f"{sample_name}.matrix.WGS.sigprofiler.tsv",
                                            header = True,
                                            index = False,
                                            sep = "\t")





@click.command()
@click.argument('mode', type=click.Choice(['matrix', 'profile']))
@click.option('--sample_name', type=str, help='Name of the sample being processed.')
@click.option('--mut_file', type=click.Path(exists=True), help='Input mutation file')
@click.option('--out_matrix', type=click.Path(), help='Output mutation matrix file')
@click.option('--method', type=click.Choice(['unique', 'multiple']), default='unique')
@click.option('--pseud', type=float, default=0)
@click.option('--per_sample', is_flag=True, help='Create a column for each sample in the input')

@click.option('--mutation_matrix', type=click.Path(exists=True), help='Mutation matrix file (for profile mode)')
@click.option('--trinucleotide_counts', type=click.Path(exists=True), help='Trinucleotide counts file (for profile mode)')
@click.option('--plot', is_flag=True, help='Generate plot and save as PDF')
@click.option('--wgs', is_flag=True, help='Store matrix of mutation counts at WGS level')
@click.option('--wgs_trinucleotide_counts', type=click.Path(exists=True), help='Trinucleotide counts file of the WGS (for profile mode if WGS active)')
@click.option('--smoothed', is_flag=True, help='Apply Bayesian smoothing to the mutation counts using a prior profile')
@click.option('--prior_profile', type=click.Path(exists=True), help='Prior profile file to use for Bayesian smoothing (required if --smoothed is set)')


@click.option('--sigprofiler', is_flag=True, help='Store the index column using the SigProfiler format')

def main(mode, sample_name, mut_file, out_matrix, method, pseud, sigprofiler, per_sample, mutation_matrix,
            trinucleotide_counts, plot, wgs, wgs_trinucleotide_counts, smoothed, prior_profile):

    if mode == 'matrix':
        click.echo(f"Running in matrix mode...")
        click.echo(f"Using the method: {method}")
        click.echo(f"Using the pseudocount: {pseud}")
        compute_mutation_matrix(sample_name, mut_file, out_matrix, method, pseud, sigprofiler, per_sample)

    elif mode == 'profile':
        click.echo(f"Running in profile mode...")
        compute_mutation_profile(sample_name, mutation_matrix, trinucleotide_counts, plot, wgs, wgs_trinucleotide_counts, sigprofiler, smoothed, prior_profile)
        click.echo("Profile computation completed.")

    else:
        click.echo("Unrecognized mode. Choose 'matrix' or 'profile'.")

if __name__ == '__main__':
    main()

