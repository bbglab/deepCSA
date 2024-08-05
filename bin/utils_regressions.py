import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd
import matplotlib.pyplot as plt


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

    # compute linear regressions for each response variable
    for response_var in name_response_vars:

        # access each covariate
        for pred_var in name_predictor_vars:

            mod = smf.ols(formula = f'{response_var} ~ {pred_var}', data = response_vars_df, missing = "drop")
            res = mod.fit()

            # extract predictor variable coefficient, intervals, and p-value
            coeffs_df.loc[response_var, pred_var] = res.params[pred_var]
            lowi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][0]
            highi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][1]
            pvals_df.loc[response_var, pred_var] = res.pvalues[pred_var]

    # correct p values for false discovery rate
    if correct_pvals:
        _, corr_pvals = fdrcorrection(pvals_df.values.flatten(), alpha = 0.05, method = 'indep', is_sorted = False)
        corr_pvals_df = pd.DataFrame(corr_pvals.reshape(pvals_df.shape), index = pvals_df.index, columns = pvals_df.columns)
    else:
        corr_pvals_df = pvals_df

    return coeffs_df, lowi_df, highi_df, corr_pvals_df

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

    # compute linear regressions for each response variable
    for response_var in name_response_vars:

            # access each covariate
            for pred_var in name_predictor_vars:

                mod = smf.mixedlm(formula = f'{response_var} ~ {pred_var}', data = response_vars_df,
                                  groups = response_vars_df[name_random_effects], missing = "drop") #raises error if NA otherwise
                res = mod.fit()

                # extract predictor variable coefficient, intervals, and p-value
                coeffs_df.loc[response_var, pred_var] = res.params[pred_var]
                lowi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][0]
                highi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][1]
                pvals_df.loc[response_var, pred_var] = res.pvalues[pred_var]

    # correct p values for false discovery rate
    if correct_pvals:
        _, corr_pvals = fdrcorrection(pvals_df.values.flatten(), alpha = 0.05, method = 'indep', is_sorted = False)
        corr_pvals_df = pd.DataFrame(corr_pvals.reshape(pvals_df.shape), index = pvals_df.index, columns = pvals_df.columns)
    else:
        corr_pvals_df = pvals_df

    return coeffs_df, lowi_df, highi_df, corr_pvals_df

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

    # compute a multivariate linear regressions for each response variable
    for response_var in name_response_vars:

        # generate formula
        formula = f'{response_var} ~ '
        for pred_var in name_predictor_vars:
            formula = f'{formula} + {pred_var}'
        formula = formula.replace("+", "", 1)

        # compute model
        mod = smf.ols(formula = formula, data = response_vars_df, missing = "drop")
        res = mod.fit()

        ## extract covar coefficient, intervals, and p-value for each predictor
        for pred_var in name_predictor_vars:
            coeffs_df.loc[response_var, pred_var] = res.params[pred_var]
            lowi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][0]
            highi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][1]
            pvals_df.loc[response_var, pred_var] = res.pvalues[pred_var]

    # correct p values for false discovery rate
    if correct_pvals:
        _, corr_pvals = fdrcorrection(pvals_df.values.flatten(), alpha = 0.05, method = 'indep', is_sorted = False)
        corr_pvals_df = pd.DataFrame(corr_pvals.reshape(pvals_df.shape), index = pvals_df.index, columns = pvals_df.columns)
    else:
        corr_pvals_df = pvals_df

    return coeffs_df, lowi_df, highi_df, corr_pvals_df

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

        # compute model
        mod = smf.mixedlm(formula = formula, data = response_vars_df, groups = response_vars_df[name_random_effects], missing = "drop")
        res = mod.fit()

        ## extract covar coefficient, intervals, and p-value for each predictor
        for pred_var in name_predictor_vars:
            coeffs_df.loc[response_var, pred_var] = res.params[pred_var]
            lowi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][0]
            highi_df.loc[response_var, pred_var] = res.conf_int().loc[pred_var][1]
            pvals_df.loc[response_var, pred_var] = res.pvalues[pred_var]

    # correct p values for false discovery rate
    if correct_pvals:
        _, corr_pvals = fdrcorrection(pvals_df.values.flatten(), alpha = 0.05, method = 'indep', is_sorted = False)
        corr_pvals_df = pd.DataFrame(corr_pvals.reshape(pvals_df.shape), index = pvals_df.index, columns = pvals_df.columns)
    else:
        corr_pvals_df = pvals_df

    return coeffs_df, lowi_df, highi_df, corr_pvals_df

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

