import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd
import matplotlib.pyplot as plt
import os

def create_metric_table(metric_df, metric_var, rows_var, cols_var,
                        rows_names, cols_names,
                        total_cols_by, total_rows_by,
                        metric_var4file, save_files_dir,
                        keep_rows_ordered = True, keep_cols_ordered = True):
    """
    Converts metric_df in a pivoted table using
    metric_var, rows_var and cols_var. Saves this
    dataframe with columns as cols_names and
    index as rows_names, regardless of any of the
    provided names is not in metric_df (thus, that
    column/row is filled with NA). Orders rows and/or
    columns if indicated according to cols_names/rows_names
    order. Also, creates two dataframes with the totals
    of the rows and the columns, according to total_by_cols
    and total_by_rows. Keeps as well all cols_names and rows_names

    Parameters
    ----------
    metric_df: pandas Dataframe
        Initial dataframe with at least 3 columns named
        as provided in metric_var, rows_var, cols_var
    metric_var: str
        Column name for the values to be used to fill
        the pivoted dataframe
    rows_var: str
        Column name for the values to be used as row names
        of the pivoted dataframe
    cols_var: str
        Column name for the values to be used as column names
        of the pivoted dataframe
    rows_names: list
        List of values to be used as the row names of the
        pivoted dataframe. Used both to subset the df if needed
        and to add rows with NA we want to keep for downstream
        analysis
    cols_names: list
        List of values to be used as the column names of the
        pivoted dataframe. Used both to subset the df if needed
        and to add columns with NA we want to keep for downstream
        analysis
    total_cols_by: str
        How to calculate the total for the columns of the dataframe.
        If "sum", "mean" or "median", calculates the total from the
        main pivoted dataframe. If other value, the total is provided
        in metric_df and total_cols_by value should correspond with the
        value used in the specific metric_df to gather all the columns
        into one value.
    total_rows_by: str
        How to calculate the total for the rows of the dataframe.
        If "sum", "mean" or "median", calculates the total from the
        main pivoted dataframe. If other value, the total is provided
        in metric_df and total_rows_by value should correspond with the
        value used in the specific metric_df to gather all the rows
    metric_var4file: str
        Name to be used in the file naming for the metric.
        TODO: maybe provide a dictionary with the equivalences instead
        (this is for the regressions, for which the file name is
        informative)
    save_files_dir: str
        Path where to store the generated files
    keep_rows_ordered: Boolean (default: True)
        Whether to keep rows ordered as provided in row_names
        If False, orders alphabetically
    keep_cols_ordered: Boolean (default: True)
        Whether to keep columns ordered as provided in col_names
        If False, orders alphabetically

    Returns
    -------
    None


    """

    ## -- MAIN DATAFRAME -- ##
    # pivot and reindex so that it contains all col_names and rows_names
    metric_df_p = metric_df.pivot(values = metric_var,
                                  index = rows_var,
                                  columns = cols_var).reindex(index = rows_names,
                                                              columns = cols_names)
    # fill NA with gene's mean
    metric_df_p = metric_df_p.apply(lambda row: row.fillna(row.mean()), axis = 1)

    # remove genes for which all the values are NA
    genes2remove = metric_df_p.loc[metric_df_p.isna().all(axis = 1)].index
    print("These genes will be removed from the analysis because no value was calculated:", genes2remove.tolist())
    metric_df_p = metric_df_p.dropna(how = 'all')

    # reorder alphabetically if provided ordered is not kept
    if not keep_rows_ordered:
        metric_df_p = metric_df_p.sort_index(axis = 0)
    if not keep_cols_ordered:
        metric_df_p = metric_df_p.sort_index(axis = 1)

    ## -- TOTAL COLUMNS DATAFRAME -- ##
    # if sum, mean or median, calculate from the pivoted df
    if total_cols_by == "sum":
        metric_cols_total_df = metric_df_p.sum(axis = 0, skipna = True).to_frame(metric_var) #skipna=True to avoid getting NA as the total when there are NAs in the df
    elif total_cols_by == "mean":
        metric_cols_total_df = metric_df_p.mean(axis = 0, skipna = True).to_frame(metric_var)
    elif total_cols_by == "median":
        metric_cols_total_df = metric_df_p.median(axis = 0, skipna = True).to_frame(metric_var)
    # otherwise, get from metric_df
    else:
        metric_cols_total_df = metric_df.loc[(metric_df[rows_var] == total_cols_by)
                                            & (metric_df[cols_var].isin(cols_names))][[cols_var, metric_var]].set_index(
                                            cols_var).reindex(index = cols_names)
        ## reorder alphabetically if provided ordered is not kept
        if not keep_cols_ordered:
            metric_cols_total_df = metric_cols_total_df.sort_index(axis = 0)

    ## -- TOTAL ROWS DATAFRAME -- ##
    # if sum, mean or median, calculate from the pivoted df
    if total_rows_by == "sum":
        metric_rows_total_df = metric_df_p.sum(axis = 1, skipna = True).to_frame(metric_var)
    elif total_rows_by == "mean":
        metric_rows_total_df = metric_df_p.mean(axis = 1, skipna = True).to_frame(metric_var)
    elif total_rows_by == "median":
        metric_rows_total_df = metric_df_p.median(axis = 1, skipna = True).to_frame(metric_var)
    # otherwise, get from metric_df
    else:
        metric_rows_total_df = metric_df.loc[(metric_df[cols_var] == total_rows_by)
                                            & (metric_df[rows_var].isin(rows_names))][[rows_var, metric_var]].set_index(
                                            rows_var).reindex(index = rows_names)
        ## reorder alphabetically if provided ordered is not kept
        if not keep_rows_ordered:
            metric_rows_total_df = metric_rows_total_df.sort_index(axis = 0)

    ## -- SAVE OR RETURN -- ##
    metric_df_p.to_csv(os.path.join(save_files_dir, f"{metric_var4file.lower()}.tsv"), sep = "\t")
    assert metric_df_p.shape == (len(rows_names), len(cols_names))
    metric_rows_total_df.to_csv(os.path.join(save_files_dir, f"{metric_var4file.lower()}.total_{rows_var.lower()}.tsv"), sep = "\t")
    assert metric_rows_total_df.shape == (len(rows_names), 1)
    metric_cols_total_df.to_csv(os.path.join(save_files_dir, f"{metric_var4file.lower()}.total_{cols_var.lower()}.tsv"), sep = "\t") # having points here may collapse with the file naming in the main function
    assert metric_cols_total_df.shape == (len(cols_names), 1)

    return None

def process_mutrate(mutrate_file, mutrate_config,
                    rows_names, cols_names,
                    save_files_dir, metric = "mutrate"):
    """
    Generates and saves pivoted dataframes of mutation rates,
    with columns as samples and rows as genes. Does the same with
    a two column table for the total of the genes and the total
    of the samples. Creates as many versions as those specified
    in mutrate_config

    Parameters
    ----------
    mutrate_file: str
        Path to the mutation rates file generated by MUTRATE
        process (deepCSA pipeline)
    mutrate_config: list
        List containing the specific versions to be created
        regarding region and mutations used
        Allowed values: all, nonproteinaffecting, proteinaffecting,
        snv, alltypes, complex-deletion-insertion-mnv, indel
    rows_names: list
        List of values to be used as the row names of the
        pivoted dataframe. Used both to subset the df if needed
        and to add rows with NA we want to keep for downstream
        analysis. If rows_names=[''], uses all the genes in mutrate_file_f
    cols_names: list
        List of values to be used as the column names of the
        pivoted dataframe. Used both to subset the df if needed
        and to add columns with NA we want to keep for downstream
        analysis
    save_files_dir: str
        Path where to store the generated files
    metric: str (default: "mutrate")
        Mutation rate metric to use, which can be mutation rate or
        mutated reads rate.
        Allowed values: mutrate, mutreadsrate

    Returns
    -------
    None
    """

    # load master mutrate file
    mutrate_df = pd.read_csv(mutrate_file, sep = "\t")

    # set table settings
    mutrate_config_fixed = []
    for val in mutrate_config:
        vals = val.split("-")
        mutrate_config_fixed.extend(vals)

    print(mutrate_config_fixed)

    # set version settings
    ## main metric
    if metric == "mutrate":
        metric = "MUTRATE_MB"
    elif metric == "mutreadsrate":
        metric = "MUTREADSRATE_MB"

    ## region
    regressions2mutrate_regions = {"all": "all",
                                   "nonproteinaffecting": "non_protein_affecting",
                                   "proteinaffecting": "protein_affecting"}
    regions = []
    for region in regressions2mutrate_regions:
        if region in mutrate_config_fixed:
            regions.append(regressions2mutrate_regions[region])

    ## muts used
    regressions2mutrate_muts = {"snv": "SNV",
                                "alltypes": "all_types",
                                "complex-deletion-insertion-mnv": "COMPLEX-DELETION-INSERTION-MNV",
                                "indel": "DELETION-INSERTION"}
    muttypes = []
    for muttype in regressions2mutrate_muts:
        if muttype in mutrate_config_fixed:
            muttypes.append(regressions2mutrate_muts[muttype])

    # create tables
    for region in regions:
        for muttype in muttypes:

            mutrate_df_f = mutrate_df.loc[(mutrate_df["REGIONS"] == region)
                                            & (mutrate_df["MUTTYPES"] == muttype)]
            mutrate_df_f = mutrate_df_f.rename({"GENE": "gene", "SAMPLE_ID": "sample"}, axis = 1)

            if rows_names == ['']:
                rows_names = [gene for gene in list(mutrate_df_f["gene"].unique()) if gene != "ALL_GENES"]

            if muttype == "DELETION-INSERTION": # to correctly name the file below
                muttype = "indel"

            create_metric_table(metric_df = mutrate_df_f,
                                metric_var = metric,
                                rows_var = "gene", cols_var = "sample",
                                rows_names = rows_names, cols_names = cols_names,
                                total_cols_by = "ALL_GENES",
                                total_rows_by = "all_samples",
                                metric_var4file = f'{metric.split("_")[0]}.{"".join(region.split("_"))}_{muttype.lower().replace("_", "")}',
                                save_files_dir = save_files_dir,
                                keep_rows_ordered = True, keep_cols_ordered = True)

    return None

def process_oncodrivefml(oncodrivefml_data, oncodrivefml_config, total_cols_by,
                         rows_names, cols_names,
                         save_files_dir):
    """
    Generates and saves pivoted dataframes of OncodriveFML metrics
    (z-score and difference between observed and expected mean
    deleteriousness scores), with columns as samples and rows as genes.
    Does the same with a two column table for the total of
    the genes and the total of the samples.
    Creates as many versions as those specified in
    oncodrivefml_config

    Parameters
    ----------
    oncodrivefml_data: str
        Path to the directory where the output files of OncodriveFML
        are stored or comma-separated string with paths to each sample's directory
    oncodrivefml_config: list
        List containing the specific versions to be created
        regarding metric, profile and significance level
        Allowed values: zscore, diffobsvsexp, allprof,
        nonprotaffprof, nosignificant, significant
    total_cols_by: str
        How to calculate the total for the columns of the dataframe.
        Can be "sum", "mean" or "median"
    rows_names: list
        List of values to be used as the row names of the
        pivoted dataframe. Used both to subset the df if needed
        and to add rows with NA we want to keep for downstream
        analysis. If rows_names=[''], uses all the genes in oncodrivefml_df_f
    cols_names: list
        List of values to be used as the column names of the
        pivoted dataframe. Used both to subset the df if needed
        and to add columns with NA we want to keep for downstream
        analysis
    save_files_dir: str
        Path where to store the generated files

    Returns
    -------
    None
    """

    # load OncodriveFML results in a df
    oncodrivefml_df = pd.DataFrame()
    oncodrivefml_data = [item.strip() for item in oncodrivefml_data.split(",")]

    ## a single path was provided with all the directories per sample
    if len(oncodrivefml_data) == 1:
        oncodrivefml_dir = oncodrivefml_data[0]
        oncodrivefml_data = [f"{oncodrivefml_dir}/{file}/{file.split('.')[0]}-oncodrivefml.tsv.gz" for file in os.listdir(oncodrivefml_dir)]
        for file in oncodrivefml_data:
            sample_df = pd.read_csv(file, sep = "\t", header = 0)
            sample_df["sample"] = file.split("/")[-2]
            oncodrivefml_df = pd.concat((oncodrivefml_df, sample_df)).reset_index(drop = True)
    ## a list of directory paths was provided
    else:
        for file_dir in oncodrivefml_data:
            sample_df = pd.read_csv(f"{file_dir}/{file_dir.split('.')[0]}-oncodrivefml.tsv.gz", sep = "\t", header = 0)
            sample_df["sample"] = file_dir
            oncodrivefml_df = pd.concat((oncodrivefml_df, sample_df)).reset_index(drop = True)

    # set version settings
    ## main metric
    metrics = []
    if "zscore" in oncodrivefml_config:
        metrics.append("Z-SCORE")
    if "diffobsvsexp" in oncodrivefml_config:
        ### compute difference between expected and observed mean deleteriousness scores
        oncodrivefml_df["DIFF-OBSvsEXP"] = oncodrivefml_df["AVG_SCORE_OBS"] - oncodrivefml_df["POPULATION_MEAN"]
        metrics.append("DIFF-OBSvsEXP")

    ## profile
    profiles = []
    regressions2fml_profile = {"allprof": ".all",
                               "nonprotaffprof": ".non_prot_aff"}
    for prof in regressions2fml_profile:
        if prof in oncodrivefml_config:
            profiles.append(regressions2fml_profile[prof])

    # create tables
    for metric_var in metrics:
        for profile in profiles:

            oncodrivefml_df_f = oncodrivefml_df.loc[oncodrivefml_df["sample"].str.contains(profile)]
            oncodrivefml_df_f["sample"] = oncodrivefml_df_f.apply(
                lambda row: row["sample"].split(".")[0], axis = 1)    # after use, remove profile info to every sample id
            oncodrivefml_df_f = oncodrivefml_df_f.rename({"GENE_ID": "gene"}, axis = 1)
            if rows_names == ['']:
                rows_names = [gene for gene in list(oncodrivefml_df_f["gene"].unique()) if gene != "ALL_GENES"]

            ## table w/ all oncodrivefml values, regardless of significance
            if "nosignificant" in oncodrivefml_config:
                create_metric_table(metric_df = oncodrivefml_df_f,
                                    metric_var = metric_var,
                                    rows_var = "gene", cols_var = "sample",
                                    rows_names = rows_names, cols_names = cols_names,
                                    total_cols_by = total_cols_by,
                                    total_rows_by = "all_samples",
                                    metric_var4file = f'oncodrivefml.{"".join(metric_var.split("-"))}_{"".join(profile[1:].split("_"))}prof_nosignificant',
                                    save_files_dir = save_files_dir,
                                    keep_rows_ordered = True, keep_cols_ordered = True)

            ## table w/ only significant values
            if "significant" in oncodrivefml_config:
                ### aproach 1: filter out non-significant values, generates NA
                # oncodrivefml_df_f.loc[(oncodrivefml_df_f["Q_VALUE"] > 0.05)
                                    # | ((oncodrivefml_df_f["Q_VALUE"].isna()) & (oncodrivefml_df_f["P_VALUE"] > 0.05)), metric_var] = np.nan                                          | ((oncodrivefml_df_f["Q_VALUE"].isna()) & (oncodrivefml_df_f["P_VALUE"] < 0.05))]
                ### aproach 2: fill non-significant values with zero
                oncodrivefml_df_f.loc[(oncodrivefml_df_f["Q_VALUE"] > 0.05)
                                    | ((oncodrivefml_df_f["Q_VALUE"].isna()) & (oncodrivefml_df_f["P_VALUE"] > 0.05)), metric_var] = 0
                create_metric_table(metric_df = oncodrivefml_df_f,
                                    metric_var = metric_var,
                                    rows_var = "gene", cols_var = "sample",
                                    rows_names = rows_names, cols_names = cols_names,
                                    total_cols_by = total_cols_by,
                                    total_rows_by = "all_samples",
                                    metric_var4file = f'oncodrivefml.{"".join(metric_var.split("-"))}_{"".join(profile[1:].split("_"))}prof_significant',
                                    save_files_dir = save_files_dir,
                                    keep_rows_ordered = True, keep_cols_ordered = True)

    return None

def process_omega(omega_data, omega_config,
                      total_cols_by,
                      rows_names, cols_names,
                      save_files_dir,
                      omega_modality = "mle",
                      global_loc = True):
    """
    Generates and saves pivoted dataframes of omega (MLE or bayes),
    with columns as samples and rows as genes.
    Does the same with a two column table for the total of
    the genes and the total of the samples.
    Creates as many versions as those specified in
    omega_config

    Parameters
    ----------
    omega_dir: str
        Path to the directory where the output files of omega
        are stored or comma-separated string with file names
    total_cols_by: str
        How to calculate the total for the columns of the dataframe.
        Can be "sum", "mean" or "median"
    rows_names: list
        List of values to be used as the row names of the
        pivoted dataframe. Used both to subset the df if needed
        and to add rows with NA we want to keep for downstream
        analysis. If rows_names=[''], uses all the genes in omega_df_f
    cols_names: list
        List of values to be used as the column names of the
        pivoted dataframe. Used both to subset the df if needed
        and to add columns with NA we want to keep for downstream
        analysis
    save_files_dir: str
        Path where to store the generated files
    omega_modality: str (default: "mle")
        Omega modality used to compute omegas. Can be "mle" or "bayes".

    Returns
    -------
    None
    """

    # load omega results in a df
    omega_df = pd.DataFrame()
    omega_data = [item.strip() for item in omega_data.split(",")]
    if len(omega_data) == 1:
        omega_dir = omega_data[0]
        omega_data = [f"{omega_dir}/{file}" for file in os.listdir(omega_dir) if omega_modality in file]
        for file in omega_data:
            try :
                sample_df = pd.read_csv(file, sep = "\t", header = 0)
            except pd.errors.EmptyDataError:
                sample_df = pd.DataFrame()
            sample_df["sample"] = file.split("/")[-1]
            omega_df = pd.concat((omega_df, sample_df)).reset_index(drop = True)
    else:
        for file in omega_data:
            try :
                sample_df = pd.read_csv(file, sep = "\t", header = 0)
            except pd.errors.EmptyDataError:
                sample_df = pd.DataFrame()
            sample_df["sample"] = file
            omega_df = pd.concat((omega_df, sample_df)).reset_index(drop = True)

    # set table settings
    omega_config_fixed = []
    for val in omega_config:
        vals = val.split("-")
        omega_config_fixed.extend(vals)

    print(omega_config_fixed)
    metric_var = "dnds"

    ## mutation impacts
    regressions2omega_impacts = {"essentialsplice": "essential_splice",
                                "essentialspliceplus": "essential_splice_plus",
                                "missense": "missense",
                                "nonsense": "nonsense",
                                "nonsynonymoussplice": "nonsynonymous_splice",
                                "truncating": "truncating",
                                "truncatingplus": "truncating_plus"}
    impacts = []
    for impact in regressions2omega_impacts:
        if impact in omega_config_fixed:
            impacts.append(regressions2omega_impacts[impact])
    print(impacts)

    ## profile
    regressions2omega_profile = {"allprof": "NoValue",
                                "nonprotaffprof": "non_prot_aff"}
    omega2regressions_profile = {"NoValue": "allprof",
                                "non_prot_aff": "nonprotaffprof"}
    profiles = []
    for prof in regressions2omega_profile:
        if prof in omega_config_fixed:
            profiles.append(regressions2omega_profile[prof])

    ## unique or multi muts
    regressions2omega_muts = {"uniquemuts": "NoValue",
                            "multimuts": "multi"}
    omega2regressions_muts = {"NoValue": "uniquemuts",
                            "multi": "multimuts"}
    uniqueormulti_muts = []
    for muts in regressions2omega_muts:
        if muts in omega_config_fixed:
            uniqueormulti_muts.append(regressions2omega_muts[muts])

    samples = list(set([".".join(file.split(".")[:2]) for file in omega_df["sample"]]))
    for profile in profiles:
        for uniqormulti in uniqueormulti_muts:
            if global_loc:
                analysis_info = [f'{sampl}.{profile}.{uniqormulti}.global_loc.tsv'.replace('.NoValue', '') for sampl in samples]
            else:
                analysis_info = [f'{sampl}.{profile}.{uniqormulti}.tsv'.replace('.NoValue', '') for sampl in samples] # this way we do the subset properly (checked)

            for impact in impacts:
                omega_df_f = omega_df.loc[(omega_df["sample"].isin(analysis_info)) &
                                            (omega_df["impact"] == impact)]
                print(omega_df_f)
                omega_df_f["sample"] = omega_df_f.apply(lambda row: row["sample"].split(".")[1], axis = 1)   # after use, remove analysis info to every sample id
                omega_df_f = omega_df_f.rename({"GENE_ID": "gene"}, axis = 1)
                if rows_names == ['']:
                    rows_names = [gene for gene in list(omega_df_f["gene"].unique()) if gene != "ALL_GENES"]

                ## table w/ all omega values, regardless of significance
                if "nosignificant" in omega_config_fixed:
                    metric_var4file = f"omega.{metric_var}_{omega_modality}_{impact.replace('_', '')}_{omega2regressions_profile[profile]}_{omega2regressions_muts[uniqormulti]}_nosignificant"
                    create_metric_table(metric_df = omega_df_f,
                                        metric_var = metric_var,
                                        rows_var = "gene", cols_var = "sample",
                                        rows_names = rows_names, cols_names = cols_names,
                                        total_cols_by = total_cols_by,
                                        total_rows_by = "all_samples",
                                        metric_var4file = metric_var4file,
                                        save_files_dir = save_files_dir,
                                        keep_rows_ordered = True, keep_cols_ordered = True)

                ## table w/ only significant omega values
                if "significant" in omega_config_fixed:
                    metric_var4file = f"omega.{metric_var}_{omega_modality}_{impact.replace('_', '')}_{omega2regressions_profile[profile]}_{omega2regressions_muts[uniqormulti]}_significant"
                    ### aproach 1: filter out non-significant values, generates NA
                    # omega_df_f = omega_df_f.loc[(omega_df_f["pvalue"] < 0.05)]
                    ### aproach 2: fill non-significant values with zero
                    omega_df_f.loc[(omega_df_f["pvalue"] > 0.05), metric_var] = 0
                    create_metric_table(metric_df = omega_df_f,
                                        metric_var = metric_var,
                                        rows_var = "gene", cols_var = "sample",
                                        rows_names = rows_names, cols_names = cols_names,
                                        total_cols_by = total_cols_by,
                                        total_rows_by = "all_samples",
                                        metric_var4file = metric_var4file,
                                        save_files_dir = save_files_dir,
                                        keep_rows_ordered = True, keep_cols_ordered = True)


    return None

def univar_linear_regression(response_vars_df, name_response_vars, name_predictor_vars, correct_pvals = True):
    """
    Calculates a univariate linear regression for several response and predictor variables
    and returns the coefficients, confidence intervals and corrected p values
    in a dataframe.

    Parameters
    ----------
    response_vars_df: pandas DataFrame
        Dataframe containing the data points of the response variables and
        predictor variable per sample. Both predictor and response variables must be
        in the dataframe columns, and samples in the rows.
    name_response_vars: list
        Response variables
    name_predictor_vars: list
        Predictor variables

    Returns
    -------
    coeffs_df: pandas DataFrame
        Dataframe containing each predictor-response coefficient
    lowi_df: pandas DataFrame
        Dataframe containing each predictor-response low confidence interval
        of the coefficient
    highi_df: pandas DataFrame
        Dataframe containing each predictor-response high confidence interval
        of the coefficient
    corr_pvals_df: pandas DataFrame
        Dataframe containing each predictor-response corrected p value of the coefficient
        (Benjamini/Hochberg method)
    """

    # initialize param model dfs to store results
    coeffs_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    lowi_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    highi_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    pvals_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    intercepts_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)

    # list of variables that need intercept passing through zero
    vars_zero_intercep = ["age_decades"]

    # compute linear regressions for each response variable
    for response_var in name_response_vars:

        # access each covariate
        for pred_var in name_predictor_vars:

            if pred_var in vars_zero_intercep:
                interc = "- 1"
            else:
                interc = "+ 1"

            print(f"MODEL formula (univariable): {response_var} ~ {pred_var} {interc}")
            mod = smf.ols(formula = f'{response_var} ~ {pred_var}', data = response_vars_df, missing = "drop")
            res = mod.fit()

            # extract predictor variable coefficient, intervals, and p-value
            coeffs_df.loc[response_var, pred_var] = res.params[pred_var]
            lowi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][0]
            highi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][1]
            pvals_df.loc[response_var, pred_var] = res.pvalues[pred_var]
            # new! extract also the intercept
            if pred_var in vars_zero_intercep:
                intercept = 0
            else:
                intercept = res.params["Intercept"]
            intercepts_df.loc[response_var, pred_var] = intercept

    # correct p values for false discovery rate
    if correct_pvals:
        _, corr_pvals = fdrcorrection(pvals_df.values.flatten(), alpha = 0.05, method = 'indep', is_sorted = False)
        corr_pvals_df = pd.DataFrame(corr_pvals.reshape(pvals_df.shape), index = pvals_df.index, columns = pvals_df.columns)
    else:
        corr_pvals_df = pvals_df

    return coeffs_df, lowi_df, highi_df, corr_pvals_df, intercepts_df

def univar_mixedeffects_linear_regression(response_vars_df, name_response_vars, name_predictor_vars, name_random_effects, correct_pvals = True):
    """
    Calculates a  univariate mixed-effects linear regression for several response and predictor variables
    and returns the coefficients, confidence intervals and p values in a dataframe.

    Parameters
    ----------
    response_vars_df: pandas DataFrame
        Dataframe containing the data points of the response variables and
        predictor variable per sample. Both predictor and response variables must be
        in the dataframe columns, and samples in the rows.
    name_response_vars: list
        Response variables
    name_predictor_vars: list
        Predictor variables
    name_random_effects: str
        Random effects we want take into account in the model

    Returns
    -------
    coeffs_df: pandas DataFrame
        Dataframe containing each predictor-response coefficient
    lowi_df: pandas DataFrame
        Dataframe containing each predictor-response low confidence interval
        of the coefficient
    highi_df: pandas DataFrame
        Dataframe containing each predictor-response high confidence interval
        of the coefficient
    corr_pvals_df: pandas DataFrame
        Dataframe containing each predictor-response corrected p value of the coefficient
        (Benjamini/Hochberg method)
    """

    # initialize param model dfs to store results
    coeffs_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    lowi_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    highi_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    pvals_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    intercepts_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)

    # list of variables that need intercept passing through zero
    vars_zero_intercep = ["age_decades"]

    # compute linear regressions for each response variable
    for response_var in name_response_vars:

            # access each covariate
            for pred_var in name_predictor_vars:

                if pred_var in vars_zero_intercep:
                    interc = "- 1"
                else:
                    interc = "+ 1"

                print(f"MODEL formula (univariable): {response_var} ~ {pred_var} {interc}")
                mod = smf.mixedlm(formula = f'{response_var} ~ {pred_var} {interc}', data = response_vars_df,
                                groups = response_vars_df[name_random_effects], missing = "drop") #raises error if NA otherwise
                res = mod.fit()

                # extract predictor variable coefficient, intervals, and p-value
                coeffs_df.loc[response_var, pred_var] = res.params[pred_var]
                lowi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][0]
                highi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][1]
                pvals_df.loc[response_var, pred_var] = res.pvalues[pred_var]
                # new! extract also the intercept
                if pred_var in vars_zero_intercep:
                    intercept = 0
                else:
                    intercept = res.params["Intercept"]
                intercepts_df.loc[response_var, pred_var] = intercept


    # correct p values for false discovery rate
    if correct_pvals:
        _, corr_pvals = fdrcorrection(pvals_df.values.flatten(), alpha = 0.05, method = 'indep', is_sorted = False)
        corr_pvals_df = pd.DataFrame(corr_pvals.reshape(pvals_df.shape), index = pvals_df.index, columns = pvals_df.columns)
    else:
        corr_pvals_df = pvals_df

    return coeffs_df, lowi_df, highi_df, corr_pvals_df, intercepts_df

def multivar_linear_regression(response_vars_df, name_response_vars, name_predictor_vars, correct_pvals = True, rules = None):
    """
    Calculates a multivariate linear regression for several response variables
    and returns the coefficients, confidence intervals and corrected p values in a dataframe.

    Parameters
    ----------
    response_vars_df: pandas DataFrame
        Dataframe containing the data points of the response variables and
        predictor variable per sample. Both predictor and response variables must be
        in the dataframe columns, and samples in the rows.
    name_response_vars: list
        Response variables
    name_predictor_vars: list
        Predictor variables

    Returns
    -------
    coeffs_df: pandas DataFrame
        Dataframe containing each predictor-response coefficient
    lowi_df: pandas DataFrame
        Dataframe containing each predictor-response low confidence interval
        of the coefficient
    highi_df: pandas DataFrame
        Dataframe containing each predictor-response high confidence interval
        of the coefficient
    corr_pvals_df: pandas DataFrame
        Dataframe containing each predictor-response corrected p value of the coefficient
        (Benjamini/Hochberg method)
    """

    # initialize param model dfs to store results
    coeffs_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    lowi_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    highi_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    pvals_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    intercepts_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)

    # list of variables that need intercept passing through zero
    vars_zero_intercep = ["age_decades"]

    # compute a multivariate linear regressions for each response variable
    for response_var in name_response_vars:

        # generate formula
        formula = f'{response_var} ~ '
        if rules == None:
            for pred_var in name_predictor_vars:
                formula = f'{formula} + {pred_var}'

        else:
            name_predictor_vars = rules[response_var]
            for pred_var in name_predictor_vars:
                formula = f'{formula} + {pred_var}'

        formula = formula.replace("+", "", 1)

        if any(pred_var in vars_zero_intercep for pred_var in name_predictor_vars):
            formula = f"{formula} -1"

        # compute model
        print(f"MODEL formula (multivariable): {formula}")
        mod = smf.ols(formula = formula, data = response_vars_df, missing = "drop")
        res = mod.fit()

        ## extract covar coefficient, intervals, and p-value for each predictor
        for pred_var in name_predictor_vars:
            coeffs_df.loc[response_var, pred_var] = res.params[pred_var]
            lowi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][0]
            highi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][1]
            pvals_df.loc[response_var, pred_var] = res.pvalues[pred_var]
            # new! extract also the intercept
            if any(pred_var in vars_zero_intercep for pred_var in name_predictor_vars):
                    intercept = 0
            else:
                intercept = res.params["Intercept"]
            intercepts_df.loc[response_var, pred_var] = intercept

    # correct p values for false discovery rate
    if correct_pvals:
        _, corr_pvals = fdrcorrection(pvals_df.values.flatten(), alpha = 0.05, method = 'indep', is_sorted = False)
        corr_pvals_df = pd.DataFrame(corr_pvals.reshape(pvals_df.shape), index = pvals_df.index, columns = pvals_df.columns)
    else:
        corr_pvals_df = pvals_df

    return coeffs_df, lowi_df, highi_df, corr_pvals_df, intercepts_df

def multivar_mixedeffects_linear_regression(response_vars_df, name_response_vars,
                                            name_predictor_vars, name_random_effects,
                                            correct_pvals = True, rules = None):
    """
    Calculates a mixed-effects linear regression for several response variables
    and returns the coefficients, confidence intervals and corrected p values in a dataframe.

    Parameters
    ----------
    response_vars_df: pandas DataFrame
        Dataframe containing the data points of the response variables and
        predictor variable per sample. Both predictor and response variables must be
        in the dataframe columns, and samples in the rows.
    name_response_vars: list
        Response variables
    name_predictor_vars: list
        Predictor variables
    name_random_effects: str
        Random effects we want take into account in the model


    Returns
    -------
    coeffs_df: pandas DataFrame
        Dataframe containing each predictor-response coefficient
    lowi_df: pandas DataFrame
        Dataframe containing each predictor-response low confidence interval
        of the coefficient
    highi_df: pandas DataFrame
        Dataframe containing each predictor-response high confidence interval
        of the coefficient
    corr_pvals_df: pandas DataFrame
        Dataframe containing each predictor-response corrected p value of the coefficient
        (Benjamini/Hochberg method)
    """

    # initialize param model dfs to store results
    coeffs_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    lowi_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    highi_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    pvals_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)
    intercepts_df = pd.DataFrame(index = name_response_vars, columns = name_predictor_vars)

    # list of variables that need intercept passing through zero
    vars_zero_intercep = ["age_decades"]

    # compute a multivariate linear regressions for each response variable
    for response_var in name_response_vars:

        # generate formula
        formula = f'{response_var} ~ '
        if rules == None:
            for pred_var in name_predictor_vars:
                formula = f'{formula} + {pred_var}'

        else:
            name_predictor_vars = rules[response_var]
            for pred_var in name_predictor_vars:
                formula = f'{formula} + {pred_var}'

        formula = formula.replace("+", "", 1)

        if any(pred_var in vars_zero_intercep for pred_var in name_predictor_vars):
            formula = f"{formula} -1"

        # compute model
        print(f"MODEL formula (multivariable): {formula}")
        mod = smf.mixedlm(formula = formula, data = response_vars_df, groups = response_vars_df[name_random_effects], missing = "drop")
        res = mod.fit()

        ## extract covar coefficient, intervals, and p-value for each predictor
        for pred_var in name_predictor_vars:
            coeffs_df.loc[response_var, pred_var] = res.params[pred_var]
            lowi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][0]
            highi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][1]
            pvals_df.loc[response_var, pred_var] = res.pvalues[pred_var]
            # new! extract also the intercept
            if any(pred_var in vars_zero_intercep for pred_var in name_predictor_vars):
                    intercept = 0
            else:
                intercept = res.params["Intercept"]
            intercepts_df.loc[response_var, pred_var] = intercept

    # correct p values for false discovery rate
    if correct_pvals:
        _, corr_pvals = fdrcorrection(pvals_df.values.flatten(), alpha = 0.05, method = 'indep', is_sorted = False)
        corr_pvals_df = pd.DataFrame(corr_pvals.reshape(pvals_df.shape), index = pvals_df.index, columns = pvals_df.columns)
    else:
        corr_pvals_df = pvals_df

    return coeffs_df, lowi_df, highi_df, corr_pvals_df, intercepts_df

def plot_linreg_coeffs(coeffs_df, lowi_df, highi_df, pvals_df,
                    ncols, nrows, title, fig_width, fig_height,
                    main_plotnames_dict, coeff_plotnames_dict,
                    pdf,
                    dot_size = 150, sharey = True,
                    sign_threshold = 0.05):
    """
    Plots as many subplots as keys in main_plotparams_dict
    with the results of linear regression coefficients,
    including confidence intervals and p value

    Parameters
    ----------
    coeffs_df: pandas DataFrame
        Dataframe containing each predictor-response coefficient
    lowi_df: pandas DataFrame
        Dataframe containing each predictor-response lower boundary
        of the confidence interval of the coefficient
    highi_df: pandas DataFrame
        Dataframe containing each predictor-response upper boundary
        of the confidence interval of the coefficient
    pvals_df: pandas DataFrame
        Dataframe containing each predictor-response p value
    ncols: int
        Number of columns in the plot grid
    nrows: int
        Number of rows in the plot grid
    title: str
        Title of the plot grid
    fig_width: int
        Width of the plot grid
    fig_height: int
        Height of the plot grid
    main_plotnames_dict: dictionary
        Dictionary with the equivalence between the
        name used in the results tables and the one that
        will be used as subplot title. Mandatory: {key: subplot_title}
    coeff_plotnames_dict: dictionary
        Dictionary with the equivalence between the
        name used in the results tables and the one that
        will be used as name for the coefficient row
        in the subplot. Mandatory: {key: [color, ylabel]}
    pdf: PdfPages object
        pdf object where the plots of the regressions
        will be stored
    dot_size: int (Default: 150)
        Size of each coefficient dot.
    sharey: Boolean (Default: True)
        Whether to share the y axis between all the subplots in the
        same plot grid row
    sign_threshold: float (default: 0.05)
        Significance threshold to highlight the coefficient with
        a black edgecolor.

    Returns
    -------
    None
    """

    # General configuration
    # config_params(12)
    fig, axs = plt.subplots(nrows, ncols, figsize = (fig_width, fig_height), sharey = sharey)
    plt.suptitle(title, fontsize = 14)
    # Make a subplot per each variable in main_plotnames_dict
    for main_var, ax in zip(main_plotnames_dict.keys(), axs.flat):

        ## plot a coefficient line for each variable in coeff_plotnames_dict
        for i, coeff_var in enumerate(coeff_plotnames_dict.keys()):

            ### check whether the coeff is significant to highlight the point in black
            if pvals_df.loc[main_var][coeff_var] < sign_threshold:
                edgecolors = "black"
            else:
                edgecolors = "face"

            ### plot dot + confidence interval
            ax.scatter(coeffs_df.loc[main_var][coeff_var], i, color = coeff_plotnames_dict[coeff_var][0],
                            s = dot_size, edgecolors = edgecolors, linewidths = 2)
            ax.hlines(i, lowi_df.loc[main_var][coeff_var], highi_df.loc[main_var][coeff_var],
                        color = coeff_plotnames_dict[coeff_var][0])

        ## add labels to the plotted coefficients
        y_ticks = [i for i in range(len(coeff_plotnames_dict.keys()))]
        y_labels = [coeff_plotnames_dict[coeff_var][1] for coeff_var in coeff_plotnames_dict]
        ax.set_yticks(y_ticks, y_labels)

        ## set zero effect in x axis
        ax.vlines(0, -1, len(coeff_plotnames_dict.keys())+0.3, ls = '--', color = 'grey')

        ## subplot config
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel('Effect size')
        ax.set_ylim(-1, len(coeff_plotnames_dict.keys())+0.5)
        ax.set_title(main_plotnames_dict[main_var])

    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.show()
    plt.close()


keyword2title = {
    # methods
    'mutrate': "Mutation rate",
    'mutreadsrate': "Mutated reads rate",
    'oncodrivefml': "OncodriveFML",
    'omega': "Omega",
    'omegagloballoc': "Omega global",
    'proportionmutepithelium': "Proportion mutated epithelium",

    # metrics
    'diffobsvsexp': "observed-expected mean deleteriousness score",
    'meandnds': "mean dNdS",
    'dnds': 'dNdS',
    'zscore': "z-score",

    # regions
    'all': 'all mutations',
    'allprof': 'all mutations profile',
    'nonproteinaffecting': 'non-protein-affecting mutations',
    'nonprotaffprof': 'non-protein-affecting profile',
    'proteinaffecting': 'protein-affecting mutations',
    # 'excess': 'excess',

    # others
    'multimuts': "mutated reads used",
    'significant': "significant values only",
    'nosignificant': "no significance threshold",
    'uniquemuts': "unique mutations used",
    'bayes': "bayes",
    "mle": "MLE",

    # omega: mutations used
    'essentialsplice': 'essential splice',
    'essentialspliceplus': 'essential splice extended',
    'missense': 'missense',
    'nonsense': 'nonsense',
    'nonsynonymoussplice': 'nonsynonymous and splice',
    'truncating': 'truncating',
    'truncatingplus': 'truncating and splice', #TODO: truncating in config.ini will also take truncatingplus

    # mutrate: mutations used
    'snv': 'SNVs',
    'alltypes': 'all mutation types',
    'complex-deletion-insertion-mnv': 'indels, MNVs and complex mutations',

    # mutated epithelium
    'lowerbound': 'lower bound CI'

                }

