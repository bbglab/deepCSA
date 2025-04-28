#!/usr/bin/env python

# -- Import libraries -- #
import os
import json
import re
import click
import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


from utils_regressions import *

# -- Auxiliary functions -- #


def linreg_analysis(data_df,
                    response_vars, predictors,
                    predictors4plot_dict, responses4plot_dict,
                    # predictors2plot,
                    plot_config_dict,
                    pdf,
                    random_effects = None, make2 = True, regression_type = "univariate", save_tables_file = None,
                    correct_pvals = True, sign_threshold = 0.05, return_results = False, rules = None):
    """
    Given a set of parameters, calculates a linear regression
    and plots the results.

    Parameters
    ----------
    data_df: pandas Dataframe
        Contains the response_vars and predictors values
        of a set of samples
    response_vars: list
        Each response value will be used to calculate an
        independent regression
    predictors: list
        Variables to be used as predictors in the linear
        regression. For univariate, one regression is calculated
        per predictor independently; for multivariate, a single
        regression is calculated for all variables together
    predictors4plot_dict: dictionary
        Contains the equivalence between the predictor names in
        data_df and those used when plotting the results
    predictors2plot: list
    responses4plot_dict: dictionary
        Contains the equivalence between the response names in
        data_df and those used when plotting the results
    plot_config_dict: dictionary
        Contains the plot configuration details, specifically the
        numbers of rows and cols per plot and the height and width
        of the figures. Also, the plot title.
    pdf: PdfPages object
        pdf object where the plots of the regressions
        will be stored
    random_effects: list (default: None)
        Contains the predictor variables that will be used as random
        effects if the model is a mixed-effects multivariate regression.
    make2: Boolean (default: True)
        Whether to generate a second grid in which each subplot contains
        the results of each predictor variable
    regression_type: str (default: univariate)
        Type of linear model to compute. By default, it does a univariate
        linear regression with all the predictors and responses.
    save_tables_dir: str (default: None)
        Path to save the results of the linear models (coefficients,
        confidence intervals, p-values) in tab-separated files.
        By default, does not save these tables.
    correct_pvals: Boolean (default: True)
        Whether to apply multiple testing correction of the p-values.
        Does the correction with the false discovery rate (FDR)
    sign_threshold: float (default: 0.05)
        FDR threshold to consider a value statistically significant.
        By default, the FDR is set to 0.05.
    return_results (TODO: write description)

    Returns
    -------
    None
    """
    print("linreg_analysis")
    # Compute the coefficients, confidence intervals and p-values of the linear regression
    if regression_type == "univariate":
        coeffs_df, lowi_df, highi_df, pvals_df, intercepts_df = univar_linear_regression(response_vars_df = data_df,
                                                                      name_response_vars = response_vars,
                                                                      name_predictor_vars = predictors,
                                                                      correct_pvals = correct_pvals)
    elif regression_type == "multivariate":
        coeffs_df, lowi_df, highi_df, pvals_df, intercepts_df = multivar_linear_regression(response_vars_df = data_df,
                                                                      name_response_vars = response_vars,
                                                                      name_predictor_vars = predictors,
                                                                      correct_pvals = correct_pvals,
                                                                      rules = rules) #TODO: implemente rules in the normal multi

    elif regression_type == "univariateME":
        coeffs_df, lowi_df, highi_df, pvals_df, intercepts_df = univar_mixedeffects_linear_regression(response_vars_df = data_df,
                                                                      name_response_vars = response_vars,
                                                                      name_predictor_vars = predictors,
                                                                      name_random_effects = random_effects,
                                                                      correct_pvals = correct_pvals)

    elif regression_type == "multivariateME":
        coeffs_df, lowi_df, highi_df, pvals_df, intercepts_df = multivar_mixedeffects_linear_regression(response_vars_df = data_df,
                                                                      name_response_vars = response_vars,
                                                                      name_predictor_vars = predictors,
                                                                      name_random_effects = random_effects,
                                                                      correct_pvals = correct_pvals,
                                                                      rules = rules) #TODO: finish documenting this option
    # Save results if applicable
    if save_tables_file != None:
        coeffs_df.to_csv(f"{save_tables_file}.coeffs.tsv", sep = "\t")
        lowi_df.to_csv(f"{save_tables_file}.lowi.tsv", sep = "\t")
        highi_df.to_csv(f"{save_tables_file}.highi.tsv", sep = "\t")
        pvals_df.to_csv(f"{save_tables_file}.pvals.tsv", sep = "\t")
        intercepts_df.to_csv(f"{save_tables_file}.intercepts.tsv", sep = "\t")


    # Plot results
    # if isinstance(predictors2plot, list):
    #     predictors4plot_dict_f = {pred:predictors4plot_dict[pred] for pred in predictors4plot_dict if pred in predictors2plot}
    # else:
        # predictors4plot_dict_f = predictors4plot_dict
    predictors4plot_dict_f = predictors4plot_dict
    ## mode 1: one subplot per response variable
    plot_linreg_coeffs(coeffs_df, lowi_df, highi_df, pvals_df,
                    ncols = plot_config_dict["plot_responses_ncols"],
                    nrows = plot_config_dict["plot_responses_nrows"],
                    title = plot_config_dict["plot_title"],
                    fig_width = plot_config_dict["fig_responses_width"],
                    fig_height = plot_config_dict["fig_responses_height"],
                    main_plotnames_dict = responses4plot_dict,
                    coeff_plotnames_dict = predictors4plot_dict_f,
                    pdf = pdf,
                    sign_threshold = sign_threshold)

    ## mode 2: one subplot per predictor variable
    if make2:
        ### change plot config dictionaries for inverse plotting of response and predictors
        predictors4plot_dict_inv = {predvar: predictors4plot_dict_f[predvar][1] for predvar in predictors4plot_dict_f}
        responses4plot_dict_inv = {respvar: ["#C4BCB7", responses4plot_dict[respvar]] for respvar in responses4plot_dict}
        plot_linreg_coeffs(coeffs_df.T, lowi_df.T, highi_df.T, pvals_df.T,
                        ncols = plot_config_dict["plot_predictors_ncols"],
                        nrows = plot_config_dict["plot_predictors_nrows"],
                        title = plot_config_dict["plot_title"],
                        fig_width = plot_config_dict["fig_predictors_width"],
                        fig_height = plot_config_dict["fig_predictors_height"],
                        main_plotnames_dict = predictors4plot_dict_inv,
                        coeff_plotnames_dict = responses4plot_dict_inv,
                        pdf = pdf,
                        sign_threshold = sign_threshold)

    if return_results:
        return coeffs_df, lowi_df, highi_df, pvals_df
    else:
        return None

def do_regression_analysis(info_row, pdf,
                           response_subplots = True, total_plot = True, response_and_total_subplots = True,
                           make2 = True, save_tables_dir = None, correct_pvals = True, sign_threshold = 0.05):
    """
    Does a univariate coupled with hypothesis-driven multivariate
    linear analysis for the association between two variables or
    more variables, one acting as response and the other/s as predictor/s.
    Initially thought to analyze how well a or several clinical variables
    (predictors) predict a mutagenesis/selection metric (response)
    computed with deepCSA pipeline.
    The specifications of the analysis are in info_row.
    The input metric file consists of a value per sample and response.

    Execute with a lambda function:
    config_df.apply(lambda row: do_regression_analysis(), axis = 1)

    Parameters
    ----------
    info_row: pandas DataFrame row
        Contains the information to perform the analysis.
        Columns: method, metric, metric_file_path,
        responses_subset, samples_subset
        predictors, predictors4plot, predictors_file_path
        random_effects
    pdf: PdfPages object
        pdf object where the plots of the regressions
        will be stored
    response_subplots: Boolean (default: True)
        Whether to create a grid in which each subplot contains
        the results per response and, if make2=True, another one
        with the results per predictor
    total_plot: Boolean (default: True)
        Whether to create a grid in which the only subplot contains
        the results of the total response variable
    response_and_total_subplots: Boolean (default: True)
        Whether to create a grid in which each subplot contains
        the results per response (including the total response variable)
        and, if make2=True, another one with the results per predictor
    make2: Boolean (default: True)
        Whether to generate a second grid in which each subplot contains
        the results of each predictor variable
    save_tables_dir: str (default: None)
        Path to save the results of the linear models (coefficients,
        confidence intervals, p-values) in tab-separated files.
        By default, does not save these tables.
    correct_pvals: Boolean (default: True)
        Whether to apply multiple testing correction of the p-values.
        Does the correction with the false discovery rate (FDR)
    sign_threshold: float (default: 0.05)
        FDR threshold to consider a value statistically significant.
        By default, the FDR is set to 0.05.

    Returns
    -------
    None
    """

    print(f'Analysis starting for {info_row["method"]}: {info_row["metric"]}')
    ## -- DATA AND PROCESSING -- ##
    # Load data
    ## metric dataframe (load as many dataframes if there is a joining rule)
    print("Loading metric data")
    if not isinstance(info_row["joining_rule"], str):
        metric_df = pd.read_csv(info_row["metric_file_path"], sep = "\t", index_col = "gene")
    else:
        files_paths = info_row["metric_file_path"].split(",")
        metrics = [set(submetrics.split("_")) for submetrics in info_row["metric"].split(",")]
        metrics_flat = [val for submetrics in metrics for val in submetrics]
        nonrepeated_metrics = [val for val in metrics_flat if metrics_flat.count(val) < 2]
        metrics_short = ["_".join(list(submetrics.intersection(nonrepeated_metrics))) for submetrics in metrics] #TODO. I think this does not do what I want, revise
        metric_df = pd.concat((pd.read_csv(file_path, sep = "\t", index_col = "gene").T.add_suffix(f"_{metric}").T for file_path, metric in zip(files_paths, metrics_short)),
                            ignore_index = False)

    metric_df = metric_df.sort_index() # sort to display the same gene order for all the analysis

    ## metric totals dataframe
    print("Loading metric data: total")
    if not isinstance(info_row["joining_rule"], str):
        metric_total_path = f"{'.'.join(info_row['metric_file_path'].split('.')[:2])}.total_sample.tsv"
        metric_total_df = pd.read_csv(metric_total_path, sep = "\t", index_col = 0, header = 0, names = ["total"])
    else:
        metric_total_paths = [f"{'.'.join(file_path.split('.')[:2])}.total_sample.tsv" for file_path in files_paths]
        metric_total_df = pd.concat((pd.read_csv(metric_total_path, sep = "\t", header = 0, names = ["total"]).add_suffix(f"_{metric}") for metric_total_path, metric in zip(metric_total_paths, metrics_short)),
                        ignore_index = False, axis = 1)

    ## predictors dataframe: needs to have categorical data numerically codified
    print("Loading predictors")
    predictors_df = pd.read_csv(info_row["predictors_file_path"], sep = "\t", index_col = "SAMPLE_ID")

    # Subset data if applicable (both responses and samples); note: I check for the dtype because checking for null value does not really work
    if isinstance(info_row["responses_subset"], list):
        responses_subset = [resp_df for resp_df in metric_df.index if any(resp_df for resp in info_row["responses_subset"] if resp in resp_df)]
        metric_df = metric_df.loc[responses_subset]
    if isinstance(info_row["samples_subset"], list):
        print("Selecting a subset of samples")
        samples_subset = list(set(info_row["samples_subset"]).intersection(set(metric_df.columns.tolist())))
        metric_df = metric_df[samples_subset]
    ## save used responses for later
    responses = metric_df.index.tolist()

    # Merge metric and predictors dataframes according to metric df index (which was subset previously)
    print("Merging all the data required for computing the models")
    metric_df = metric_df.T # transpose to merge with the other tables
    metric_predictors_df = metric_df.merge(predictors_df, right_index = True, left_index = True, how = "left").merge(
        metric_total_df, right_index = True, left_index = True, how = "left")
    metric_predictors_df = metric_predictors_df.dropna(how = "all", axis = 1) # drop col if all values are NA to avoid errors (TODO: is this needed?)

    ## -- DEFINE GENERAL SETTINGS FOR THE PLOTS -- ##
    print("Defining general plotting settings")
    # Number of columns and rows in the grid
    plot_config = {}
    plot_config["plot_responses_ncols"] = plot_config["plot_predictors_ncols"] = 4 # keep always four columns (for now)
    plot_config["plot_responses_nrows"] = len(responses) // 4 + 1
    plot_config["plot_predictors_nrows"] = len(info_row["predictors"]) // 4 + 1
    plot_config["fig_responses_height"] = 3.5 * plot_config["plot_responses_nrows"]
    plot_config["fig_responses_width"] = 4 * plot_config["plot_responses_ncols"]
    plot_config["fig_predictors_height"] = 7 * plot_config["plot_predictors_nrows"]
    plot_config["fig_predictors_width"] = 4 * plot_config["plot_predictors_ncols"]

    # Informative title for grid
    attributes_list = []
    for attrib in re.split(r'[,_]', info_row["metric"]):
        attributes_list.append(keyword2title[attrib])
    attributes_list = list(set(attributes_list))

    ## -- DO UNIVARIATE REGRESSIONS -- ##
    if isinstance(info_row["random_effects"], str):
        regression_modality = "univariateME"
    else:
        regression_modality = "univariate"
    print(f"Starting univariate analysis. Modality: {regression_modality}")

    # Define path to store the results of the univariate analysis if applicable
    if save_tables_dir != None:
        if not os.path.isdir(save_tables_dir):
            os.makedirs(save_tables_dir)
        save_tables_file = os.path.join(save_tables_dir,
                                        f"{info_row['method']}.{info_row['metric']}.{regression_modality}")
        save_tables_file = save_tables_file.replace(",", "_")
    else:
        save_tables_file = None

    # Define the title of this grid: each feature of the analysis goes in a separated line (make sure they are not cut in the pdf)
    print("Defining information to include in the title of each subplot")
    attributes_list_ = attributes_list.copy()
    attributes_list_.append(regression_modality)
    plot_config["plot_title"] = f"{keyword2title[info_row['method']]} ("
    for attrib in attributes_list_:
        plot_config["plot_title"] += f"{attrib},\n"
    plot_config["plot_title"] = f'{plot_config["plot_title"][:-2]})\n'

    # Create plotting information for each modality that is ON
    print("Defining plotting settings for each activated plotting modality")
    plotting_info_dict = {}
    if response_subplots: # TODO: check this one is working
        responses4plot_dict = {response:response for response in responses}
        plotting_info_dict["response_subplots"] = [responses, responses4plot_dict]
    if total_plot: # TODO: check this one is working
        responses4plot_dict = {}
        total_columns = [col for col in metric_predictors_df.columns if "total" in col]
        for total_col in total_columns:
            responses4plot_dict[total_col] = f'All responses {", ".join(total_col.split("_")[1:])}'
        plotting_info_dict["total_plot"] = [total_columns, responses4plot_dict]
    if response_and_total_subplots:
        responses4plot_dict = {response:response for response in responses}
        total_columns = [col for col in metric_predictors_df.columns if "total" in col]
        for total_col in total_columns:
            responses4plot_dict[total_col] = f'All responses {", ".join(total_col.split("_")[1:])}'
        plotting_info_dict["response_and_total_plot"] = [responses + total_columns, responses4plot_dict]
    if not plotting_info_dict:
        return print("None of the plotting modalities was selected. Aborting analysis!")

    # Make plots
    print("Computing regressions...")
    #TODO: doing it like this, we calculate the models again per each plotting modality, which is nonsense
    for mode in plotting_info_dict:
        print(f"Plotting modality: {mode}")
        coeffs_df, lowi_df, highi_df, pvals_df = linreg_analysis(metric_predictors_df,
                                                                response_vars = plotting_info_dict[mode][0],
                                                                predictors = info_row["predictors"],
                                                                predictors4plot_dict = info_row["predictors4plot"],
                                                                # predictors2plot = info_row["predictors2plot"],
                                                                responses4plot_dict = plotting_info_dict[mode][1],
                                                                plot_config_dict = plot_config,
                                                                pdf = pdf,
                                                                random_effects = info_row["random_effects"],
                                                                make2 = make2,
                                                                regression_type = regression_modality,
                                                                save_tables_file = save_tables_file,
                                                                correct_pvals = correct_pvals,
                                                                sign_threshold = sign_threshold,
                                                                return_results = True,
                                                                )

        ## -- DO MULTIVARIATE REGRESSIONS -- ##
        print("Testing significant effects with multivariate analysis")
        if isinstance(info_row["random_effects"], str):
            regression_modality = "multivariateME"
        else:
            regression_modality = "multivariate"

        # Redefine path to store the results of the multivariate analysis if applicable
        if save_tables_dir != None:
            save_tables_file = os.path.join(save_tables_dir,
                                            f"{info_row['method']}.{info_row['metric']}.{regression_modality}")
            save_tables_file = save_tables_file.replace(",", "_")
        else:
            save_tables_file = None

        # Redefine the title of this grid: each feature of the analysis goes in a separated line (make sure they are not cut in the pdf)
        print("Redefining information to include in the title of each subplot")
        attributes_list_ = attributes_list.copy()
        attributes_list_.append(regression_modality)
        plot_config["plot_title"] = f"{keyword2title[info_row['method']]} ("
        for attrib in attributes_list_:
            plot_config["plot_title"] += f"{attrib},\n"
        plot_config["plot_title"] = f'{plot_config["plot_title"][:-2]})\n'

        print("Defining tests based on user defined rules and univariate analysis results")
        def set_multi_tests(row, predictors, sign_threshold):
            predictors2test = []
            for pred in predictors:
                if row[pred] < sign_threshold:
                    predictors2test.append(pred)

            if predictors2test:
                row["predictors2test"] = predictors2test

            return row

        pvals_df = pvals_df.apply(lambda row: set_multi_tests(row, info_row["predictors"], sign_threshold), axis = 1)
        if "predictors2test" not in pvals_df.columns:
            print("No significant effects found in the univariate analysis. Skipping multivariate")
        else:
            multi_tests = pvals_df.loc[~pvals_df["predictors2test"].isna()]["predictors2test"].to_dict()
            multi_rules = [set(info_row["multivariate_rules"][rule].split(", ")) for rule in info_row["multivariate_rules"]]

            # Function to update the tests
            def update_tests():
                changed = False
                for test in multi_tests:
                    test_predictors = set(multi_tests[test])
                    for rule in multi_rules:
                        if rule.intersection(test_predictors):
                            new_predictors = rule.union(test_predictors)
                            if new_predictors != test_predictors:
                                multi_tests[test] = list(new_predictors)
                                changed = True
                return changed

            # Keep updating the tests until no more changes occur
            while update_tests():
                pass

            # for rule in info_row["multivariate_rules"]: #TODO: rules have an effect in each other, order matters here and should not
            #     rule = set(info_row["multivariate_rules"][rule].split(", "))
            #     for test in multi_tests:
            #         test_predictors = set(multi_tests[test])
            #         if rule.intersection(test_predictors):
            #             multi_tests[test] = list(rule.union(test_predictors))

            ### delete test for which a single variable remains to be tested
            tests2delete = [test for test in multi_tests if len(multi_tests[test]) == 1]
            for test in tests2delete:
                del multi_tests[test]

            predictors_multi = list(set([val for vals in multi_tests.values() for val in vals]))
            # if isinstance(info_row["predictors2plot"], list):
            #     predictors2plot_multi = [pred for pred in info_row["predictors2plot"] if pred in predictors_multi]
            # else:
            #     predictors2plot_multi = predictors_multi
            responses_multi =  multi_tests.keys()
            predictors4plot_multi = {pred: info_row["predictors4plot"][pred] for pred in info_row["predictors4plot"] if pred in predictors_multi}
            responses4plot_multi = {resp: plotting_info_dict[mode][1][resp] for resp in plotting_info_dict[mode][1] if resp in multi_tests}
            print("Calculating multivariate models")
            linreg_analysis(metric_predictors_df,
                            response_vars = responses_multi,
                            predictors = predictors_multi,
                            predictors4plot_dict = predictors4plot_multi,
                            responses4plot_dict = responses4plot_multi,
                            # predictors2plot = predictors2plot_multi,
                            plot_config_dict = plot_config,
                            pdf = pdf,
                            random_effects = info_row["random_effects"],
                            make2 = make2,
                            regression_type = regression_modality,
                            save_tables_file = save_tables_file,
                            correct_pvals = correct_pvals,
                            sign_threshold = sign_threshold,
                            return_results = True,
                            rules = multi_tests
                            )

            print("Analysis finished. Moving to next metric")


    return None


# -- Main function options -- #

@click.command()

@click.option('--pdf_path',
	 		  '-pdf',
              required = True,
              help = "path to the pdf file where the regression plots will be stored")

@click.option('--metric',
              required = True,
              default = "",
              help = "")

@click.option('--data',
              required = True,
              default = "",
              help = "")

@click.option('--metric_params',
              required = False,
              default = "",
              help = "")

@click.option('--responses_subset',
              required = False,
              default = "",
              help = "")

@click.option('--responses_excluded',
              required = False,
              default = "",
              help = "")

@click.option('--samples_subset',
              required = False,
              default = "",
              help = "")

@click.option('--predictors_file',
              required = False,
              default = "",
              help = "")

@click.option('--predictors_plot_config',
              required = False,
              default = "",
              help = "")

@click.option('--random_effects_vars',
              required = False,
              default = "",
              help = "")

@click.option('--multipletesting_join',
              required = False,
              default = "",
              help = "")

@click.option('--multivariate_rules',
              required = False,
              default = "",
              help = "")

@click.option('--response_subplots',
              type = click.BOOL,
              default = False,
              show_default = True,
              required = False,
              help = "whether to generate per response and/or per predictor subplots")

@click.option('--total_plot',
              type = click.BOOL,
              default = False,
              show_default = True,
              required = False,
              help = "whether to generate a single total subplot")

@click.option('--response_and_total_subplots',
              type = click.BOOL,
              default = True,
              show_default = True,
              required = False,
              help = "whether to generate per response and/or per predictor subplots including a total subplot")

@click.option('--make2',
              type = click.BOOL,
              default = True,
              show_default = True,
              required = False,
              help = "whether to generate a second grid in which each subplot contains the results of each predictor variable")

@click.option('--save_tables_dir',
	 		  '-save',
              required = False,
              help = "path to the directory where to store the results of the linear models")

@click.option('--correct_pvals',
              type = click.BOOL,
              default = True,
              show_default = True,
              required = False,
              help = "whether to apply multiple testing correction of the p-values")

@click.option('--sign_threshold',
              required = False,
              show_default = True,
              default = 0.05,
              help = "FDR threshold to consider a value statistically significant")

@click.option('--missing_values',
              required = False,
              show_default = True,
              default = 'ignore',
              help = "How to deal with missing values in the regressions")

# -- Main function  -- #

def main(pdf_path, metric,
        data, metric_params,
        responses_subset, responses_excluded, samples_subset,
        predictors_file, predictors_plot_config, random_effects_vars,
        multipletesting_join, multivariate_rules, response_subplots,
        total_plot, response_and_total_subplots, make2, save_tables_dir,
        correct_pvals, sign_threshold,
        missing_values
        ):
    """

    Generates the input dataframe for do_regression_analysis
    based on the information provided in config.ini. Opens
    a pdf to save the results.

    Parameters
    ----------
    config_file: str
        Path to the config.ini file.
    pdf_path: str
        Path to the pdf where the regression plots
        will be stored
    response_subplots: Boolean (default: True)
        Whether to create a grid in which each subplot contains
        the results per response and, if make2=True, another one
        with the results per predictor
    total_plot: Boolean (default: True)
        Whether to create a grid in which the only subplot contains
        the results of the total response variable
    response_and_total_subplots: Boolean (default: True)
        Whether to create a grid in which each subplot contains
        the results per response (including the total response variable)
        and, if make2=True, another one with the results per predictor
    make2: Boolean (default: True)
        Whether to generate a second grid in which each subplot contains
        the results of each predictor variable
    save_tables_dir: str (default: None)
        Path to save the results of the linear models (coefficients,
        confidence intervals, p-values) in tab-separated files.
        By default, does not save these tables.
    correct_pvals: Boolean (default: True)
        Whether to apply multiple testing correction of the p-values.
        Does the correction with the false discovery rate (FDR)
    sign_threshold: float (default: 0.05)
        FDR threshold to consider a value statistically significant.
        By default, the FDR is set to 0.05.

    Returns
    -------
    """

    print("Starting association's analysis...\n")
    print("### 1. PREPARE ANALYSIS INFORMATION ###\n")

    # Initialize config dataframe
    config_df = pd.DataFrame(columns = ["method", "metric", "metric_file_path"])
    new_row = {}

    # Load configuration data and create a dictionary with the information
    print("\tGenerating config dataframe with analysis information according to the rules provided as command line arguments...")
    config = {
        "metrics.data": {metric: data},
        "metrics.config": {metric: metric_params},
        "variables": {
            "responses_subset": responses_subset,
            "responses_excluded": responses_excluded,
            # "responses_excluded": '{"omega.missense":["CDKN1A","ARID1A","KMT2C"],"omega.truncating":["TP53","RB1"]}',
            "samples_subset": samples_subset,
            "predictors_file": predictors_file,
            "predictors": predictors_plot_config,
            # "predictors": '{"had_prior_chemotherapy":["#FAEC82","Had prior chemotherapy"],"history_drinking":["#C5FA8F","Was/is drinker"],"history_smoking":["#B27DFA","Was/is smoker"],"bmi_scaled":["#FACC79","BMI 0-1 scaled"],"is_male":["#9EC9F9","Is Male"],"age_decades":["#F7645E","Age decade"]}',
            "random_effects": random_effects_vars
        },
        "advance": {
            "multipletesting_join": multipletesting_join,
            # "multipletesting_join": '{"multipletest_rule1":["truncating|nonsynonymoussplice","allprof","uniquemuts","nosignificant"]}',
            "multivariate_rules": multivariate_rules
            # "multivariate_rules": '{"multirule3":"age_decades, bmi_scaled","multirule1":"history_smoking, is_male","multirule2":"history_smoking, age_decades"}'
        }
        }
    print("\tConfig dictionary created\n")

    # Iterate through all metrics to be analyzed and create one row per analysis modality
    print("\tStart populating config dataframe. Input files for each method will be generated applying subsetting rules for samples and responses")

    ## Create directory to store input files
    os.mkdir("inputs")

    ## Load subsetting rules
    print("\tThis subsetting rules will be applied when creating the input files:")
    samples = [var.strip() for var in config["variables"]["samples_subset"].split(",")]
    if samples == ['']:
        print(os.getcwd())
        predictors_df = pd.read_csv(config["variables"]["predictors_file"], sep = "\t", index_col = "SAMPLE_ID")  #TODO: not ideal this df is loaded here and also in the function that executes the analysis
        samples = predictors_df.index.unique()
        print(f'\t\tSubset of samples: NA (all samples in {config["variables"]["predictors_file"]} will be used)')
    else:
        print(f"\t\tSubset of samples: {samples}")

    responses = [var.strip() for var in config["variables"]["responses_subset"].split(",")]
    if responses == ['']:
        print("\t\tSubset of responses: NA (all responses calculated in each method will be used)")
    else:
        print(f"\t\tSubset of responses: {responses}")

    responses_excluded = json.loads(config["variables"]["responses_excluded"])

    ## Go through each method
    for method in config["metrics.config"]:
        method_config = [conf.strip() for conf in config["metrics.config"][method].split(",")]
        method_data = config["metrics.data"][method]
        ### follow only those for which there is data source and params in the config dict
        if method_config == [""] or method_data == "":
            print(f"\tMETHOD -> {method}: data source and/or method parameters not submitted. Analysis will not be done for this method.")
        else:
            print(f"\tMETHOD -> {method}: data source and method parameters were submitted. Input files will be generated and analysis performed for this method.")

            ### create subdirectory to store input files for this method
            method_directory = f"inputs/{method}"
            os.mkdir(method_directory)

            ### create input file diferently depending on the method
            if method == "mutrate":
                process_mutrate(mutrate_file = config["metrics.data"]["mutrate"],
                                mutrate_config = method_config,
                                rows_names = responses,
                                cols_names = samples,
                                save_files_dir = method_directory,
                                missing_values_treatment = missing_values,
                                metric = "mutrate")
            elif method == "oncodrivefml":
                process_oncodrivefml(oncodrivefml_data = config["metrics.data"]["oncodrivefml"],
                                    oncodrivefml_config = method_config,
                                    total_cols_by = "mean",
                                    rows_names = responses,
                                    cols_names = samples,
                                    missing_values_treatment = missing_values,
                                    save_files_dir = method_directory)
            elif method == "omega":
                process_omega(omega_data = config["metrics.data"]["omega"],
                            omega_config = method_config,
                            total_cols_by = "mean",
                            rows_names = responses,
                            cols_names = samples,
                            save_files_dir = method_directory,
                            missing_values_treatment = missing_values,
                            omega_modality = "mle",
                            global_loc = False)

            elif method == "omegagloballoc":
                process_omega(omega_data = config["metrics.data"]["omegagloballoc"],
                            omega_config = method_config,
                            total_cols_by = "mean",
                            rows_names = responses,
                            cols_names = samples,
                            save_files_dir = method_directory,
                            missing_values_treatment = missing_values,
                            omega_modality = "mle",
                            global_loc = True)

            ## for each metric, there are several calculations; take only those specify in the config file
            ### if method_config=[""] because no value was submitted, all the calculation modalities will be loaded
            print(os.listdir(method_directory))
            for file in os.listdir(method_directory):
                print(file)
                if all(
                    any(sub_conf in file for sub_conf in conf.split("-"))
                    for conf in method_config
                    ):
                    if "total_gene" not in file and "total_sample" not in file: # keep a single file per metric
                        print(file)
                        new_row["method"] = file.split(".")[0]
                        new_row["metric"] = file.split(".")[1]
                        new_row["metric_file_path"] = os.path.join(method_directory, file)
                        for k in responses_excluded:
                            k_elems = k.split(".") # should contain only two elements: method.keyword-for-metric
                            if k_elems[0] == new_row["method"]:
                                metric_elems = new_row["metric"].split("_")
                                for elem in metric_elems:
                                    if elem == k_elems[1]:
                                        metric_df2edit = pd.read_csv(new_row["metric_file_path"], sep = "\t")
                                        genes2exclude = responses_excluded[k]
                                        metric_df2edit = metric_df2edit.loc[~metric_df2edit["gene"].isin(genes2exclude)]
                                        os.remove(new_row["metric_file_path"])
                                        metric_df2edit.to_csv(new_row["metric_file_path"], sep = "\t")

                        print(f'Method: {new_row["method"]}, metric: {new_row["metric"]}')
                        print("Loading general config...")
                        ## add general configuration to all rows (needs to be done per row because it is not possible to assign lists to full columns in pandas)
                        new_row["responses_subset"] = [var.strip() for var in config["variables"]["responses_subset"].split(",")]
                        new_row["samples_subset"] = [var.strip() for var in config["variables"]["samples_subset"].split(",")]
                        new_row["predictors_file_path"] = config["variables"]["predictors_file"]
                        new_row["predictors"] = [var for var in json.loads(config["variables"]["predictors"]).keys()]
                        new_row["predictors4plot"] = json.loads(config["variables"]["predictors"])
                        new_row["random_effects"] = config["variables"]["random_effects"] # a single random effect can be provided
                        new_row["multivariate_rules"] = json.loads(config["advance"]["multivariate_rules"])
                        # new_row["predictors2plot"] = [var.strip() for var in config["advance"]["predictors2plot"].split(",")]

                        config_df = pd.concat([config_df, pd.DataFrame([new_row])], ignore_index = True)

    ## join metrics that should be analyzed together for the multiple testing correction
    print("Defining if there are any joining rules")
    joining_rules = json.loads(config["advance"]["multipletesting_join"])
    config_df["joining_rule"] = config_df.index.to_series().apply(lambda x: f'norule_{x}')
    if joining_rules:
        for rule in joining_rules:
            mask = pd.Series([True] * len(config_df))
            for subrule in joining_rules[rule]:
                mask = config_df["metric_file_path"].str.contains(subrule) # format: ["subrule1|subrule2|subrule3", "subrule4"]
            config_df.loc[mask, "joining_rule"] = rule
        config_df = config_df.groupby("joining_rule", dropna = False).agg({'method': lambda x: x.iloc[0] if not x.isnull().all() else np.nan, # choose first instance, as it is always the same
                                                        'metric': ','.join,
                                                        'metric_file_path': ','.join,
                                                        'responses_subset': lambda x: x.iloc[0] if not x.isnull().all() else np.nan,
                                                        'samples_subset': lambda x: x.iloc[0] if not x.isnull().all() else np.nan,
                                                        'predictors_file_path': lambda x: x.iloc[0] if not x.isnull().all() else np.nan,
                                                        'predictors': lambda x: x.iloc[0] if not x.isnull().all() else np.nan,
                                                        'predictors4plot': lambda x: x.iloc[0] if not x.isnull().all() else np.nan,
                                                        'random_effects': lambda x: x.iloc[0] if not x.isnull().all() else np.nan,
                                                        'multivariate_rules': lambda x: x.iloc[0] if not x.isnull().all() else np.nan,
                                                        # 'predictors2plot': lambda x: x.iloc[0] if not x.isnull().all() else np.nan
                                                        }
                                                        ).reset_index()
        config_df.loc[config_df["joining_rule"].str.contains("norule_"), "joining_rule"] = np.nan

    ### change to NA any empty string, list or dictionary
    def replace_with_nan(value):
        if value == "" or value == [""] or value == {"", ""}:
            return np.nan
        return value
    config_df = config_df.applymap(replace_with_nan)

    ### check if columns that cannot have a NA have it
    print("Checking if all mandatory values were submitted")
    mandatory_cols = ["predictors_file_path", "predictors"]
    for col in mandatory_cols:
        if config_df[col].isnull().all():
            return print(f"Analysis aborted. You need to specify a value for {col} to proceed")

    print("DONE!")

    print("2. Performing analysis per config group of settings...")
    print(config_df["metric"].unique())
    # Open pdf and do the analysis
    with PdfPages(pdf_path) as pdf:
        config_df.apply(lambda row: do_regression_analysis(info_row = row,
                                                            pdf = pdf,
                                                            response_subplots = response_subplots,
                                                            total_plot = total_plot,
                                                            response_and_total_subplots = response_and_total_subplots,
                                                            make2 = make2,
                                                            save_tables_dir = save_tables_dir,
                                                            correct_pvals = correct_pvals,
                                                            sign_threshold = sign_threshold),
                                                        axis = 1)

    print("All analysis finished. Bye! :)")
    return None


if __name__ == "__main__":
    main()

