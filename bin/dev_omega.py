import os
from functools import reduce, partial

import numpy as np
import pandas as pd

import numpy as np
from scipy.optimize import root_scalar, minimize, minimize_scalar, brentq
from scipy.stats import chi2, norm

import statsmodels.api as sm
import statsmodels.formula.api as smf

SEED = 31

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['TF_DETERMINISTIC_OPS'] = 'true'
os.environ['TF_CUDNN_DETERMINISTIC']= '1'


def get_reparameterized_negative_binomial(mean, overdispersion):
    """
    This parameterization is consistent with:
    variance = mean + overdispersion * mean^2
    reference: https://github.com/tensorflow/probability/issues/372
    """
    import tensorflow_probability.python.distributions as tfd

    total_count = 1 / overdispersion
    p = mean / (mean + total_count)
    return tfd.NegativeBinomial(total_count, probs=p)


def get_synthetic_data(l, true_omega, true_overdispersion):

    neg_binom = get_reparameterized_negative_binomial(l * true_omega, true_overdispersion)
    return neg_binom.sample()


def get_bounds(profile_objective, omega_hat, lower_bound=0.0, upper_bound=5e2):
    """
    Computes Wilks' Profile Likelihood CIs for a single parameter.

    goal is to find the roots of the function profile_objective(w)
    profile_objective(w) is something of the form:
    ll(w) - (ll(w_hat) - 1/2*chi2(1).ppf(0.95))
    w belongs to the CI iff profile_objective(w) >= 0

    Parameters:
        profile_objective : log-likelihood profile objective function f(omega) -> float
        omega_hat      : float, MLE of dN/dS
        lower_bound    : float, physical lower bound (default 0.0 for dN/dS)
    """

    # -------------------------------------------------------------
    # 1. Lower Bound 
    # -------------------------------------------------------------
    
    # geometrically shrink lower bound a profile_objective(a) < 0

    a = omega_hat
    while profile_objective(a) >= 0:
        a /= 2
        if a < 1e-2:
            ci_lower = lower_bound
            break

    else:

        res_lower = root_scalar(profile_objective, bracket=[a, omega_hat], method='brentq')
        ci_lower = res_lower.root

        
    # -------------------------------------------------------------
    # 2. Upper Bound 
    # -------------------------------------------------------------
    
    # geometrically expand upper bound b until profile objective(b) < 0
    
    b = omega_hat + 0.1
    while profile_objective(b) >= 0:
        b *= 2
        if b > upper_bound:
            return ci_lower, upper_bound
            
    res_upper = root_scalar(profile_objective, bracket=[omega_hat, b], method='brentq')
    ci_upper = res_upper.root

    return ci_lower, ci_upper


def fit_glm_negative_binomial(formula, data, offset_col=None):

    offset = data[offset_col] if offset_col else None
    
    # 1. Define negative log-likelihood as a function of log(alpha)
    def neg_log_likelihood(log_alpha):
        alpha = np.exp(log_alpha)  # Ensures alpha > 0
        family = sm.families.NegativeBinomial(alpha=alpha)
        model = smf.glm(formula, data=data, family=family, offset=np.log(offset))
        res = model.fit()
        return -res.llf

    # 2. Find the alpha that maximizes log-likelihood
    opt = minimize_scalar(neg_log_likelihood, bounds=(-10, 5), method='bounded')
    best_alpha = float(np.exp(opt.x))

    # 3. Fit final GLM with the optimal alpha
    final_family = sm.families.NegativeBinomial(alpha=best_alpha)
    final_model = smf.glm(formula, data=data, family=final_family, offset=np.log(offset))
    final_results = final_model.fit(maxiter=1000)
    
    return final_results, best_alpha


class Background:

    def __init__(self, data):

        # data is a pandas dataframe with columns: 'offset_syn', 'n_syn', 'covariates'
        self.data = data
        self.model = None
        self.results = None

    def fit(self, verb=True, covariates=False, samples=False):

        # negative binomial regression with log link function

        if covariates:

            joinplus = lambda x, y: x + ' + ' + y
            pcs = reduce(joinplus, [f'PC{i}' for i in range(1, 21)])
            formula = f"n_syn ~ {pcs}"
            
            if samples:
                formula = f"n_syn ~ {pcs} + sample_id"
        
        else:

            formula = "n_syn ~ 1"
            if samples:
                formula = "n_syn ~ 1 + sample_id"
    
        self.results, overdispersion = fit_glm_negative_binomial(formula, self.data, 'offset_syn')
        
        if verb:
            print(self.results.summary2())
        
        gene_coefficients = self.results.params.filter(like='gene_id')
        sample_coefficients = self.results.params.filter(like='sample_id')

        # mean predictions and confidence intervals for the fitted model

        pred_obj = self.results.get_prediction(
            self.data, offset=np.log(self.data['offset_syn']))
        pred_df = pred_obj.summary_frame(alpha=0.05)
        augmented_data = pd.concat([self.data.reset_index(drop=True), pred_df], axis=1)

        # gamma prior parameters: shape (alpha) and rate (beta)

        augmented_data['alpha'] = 1 / overdispersion
        augmented_data['beta'] = 1 / (augmented_data['mean'] * overdispersion)

        return augmented_data, overdispersion, gene_coefficients, sample_coefficients


class Omega:
    
    def __init__(self, n_syn, n_mis, n_trunc, offset_syn, offset_mis, offset_trunc, alpha, beta, mu_syn):

        self.tf, self.tfd = self._lazy_import_tf()
        self.theta = 1 / alpha
        self.alpha = alpha
        self.beta = beta
        self.n_syn = n_syn
        self.n_mis = n_mis
        self.n_trunc = n_trunc
        self.offset_syn = offset_syn
        self.offset_mis = offset_mis
        self.offset_trunc = offset_trunc
        self.mu_syn = mu_syn
        self.omega_hat = None
        
    def _lazy_import_tf(self):

        """Lazy import TensorFlow only when needed"""
        import tensorflow as tf
        import tensorflow_probability as tfp

        tfd = tfp.distributions
        
        tf.random.set_seed(SEED)
        tfp.util.SeedStream(SEED, salt="random_beta")

        return tf, tfd

    def optimal_mu(self, n_obs, n_exp_rel):

        # Compute MAP estimate of mutation density (à la dNdScv)
        
        return (n_obs + self.alpha - 1) / (n_exp_rel + self.beta)


    def fit_dndscv(self, confidence_level=0.95):

        n = self.n_syn + self.n_mis + self.n_trunc
        gamma_model = self.tfd.Gamma(concentration=self.alpha, rate=self.beta)
        
        # ll_joint: syn is neutral, mis and trunc are possibly under selection

        neutral_obs = self.n_syn
        map_mu_joint = self.optimal_mu(neutral_obs, 1)
        t_hat_joint = map_mu_joint / self.offset_syn
        se_log = 1 / np.sqrt(self.n_syn + self.alpha - 1)  # SE(log t_hat) ~ (1/t_hat) * SE(t_hat), with SE(t_hat) via observed Fisher information
        t_hat_low, t_hat_high = self.ci_log_wald(t_hat_joint, se_log)

        omega_mis = self.n_mis / (self.offset_mis * t_hat_joint)
        omega_trunc = self.n_trunc / (self.offset_trunc * t_hat_joint)

        def ll_joint_func(w1, w2):
            mu = self.offset_syn * t_hat_joint + self.offset_mis * t_hat_joint * w1 + self.offset_trunc * t_hat_joint * w2
            poisson_model = self.tfd.Poisson(mu)
            return poisson_model.log_prob(n) + gamma_model.log_prob(map_mu_joint)

        ll_joint = ll_joint_func(omega_mis, omega_trunc)

        # ll_0: all mutations are meant to be neutral

        neutral_obs = self.n_syn + self.n_mis + self.n_trunc
        k = (self.offset_syn + self.offset_mis + self.offset_trunc) / self.offset_syn
        map_mu = self.optimal_mu(neutral_obs, k)
        poisson_model = self.tfd.Poisson(map_mu)
        ll_0 = poisson_model.log_prob(n) + gamma_model.log_prob(map_mu)

        # ll_mis: syn and mis are neutral, trunc is possibly under selection

        neutral_obs = self.n_syn + self.n_mis
        k = (self.offset_syn + self.offset_mis) / self.offset_syn
        map_mu_mis = self.optimal_mu(neutral_obs, k)
        t_hat_mis = map_mu_mis / self.offset_syn
        omega_trunc = self.n_trunc / (self.offset_trunc * t_hat_mis)
        
        def ll_mis_func(w1, w2):
            mu = self.offset_syn * k * t_hat_mis * w1 + self.offset_trunc * t_hat_mis * w2
            poisson_model = self.tfd.Poisson(mu)
            return poisson_model.log_prob(n) + gamma_model.log_prob(map_mu)
        
        ll_mis = ll_mis_func(1, omega_trunc)

        # ll_trunc: syn and trunc are neutral, mis is possibly under selection

        neutral_obs = self.n_syn + self.n_trunc
        k = (self.offset_syn + self.offset_trunc) / self.offset_syn
        map_mu_trunc = self.optimal_mu(neutral_obs, k)
        t_hat_trunc = map_mu_trunc / self.offset_syn
        omega_mis = self.n_mis / (self.offset_mis * t_hat_trunc)
        
        def ll_trunc_func(w1, w2):
            mu = self.offset_mis * t_hat_trunc * w1 + self.offset_syn * k * t_hat_trunc * w2
            poisson_model = self.tfd.Poisson(mu)
            return poisson_model.log_prob(n) + gamma_model.log_prob(map_mu_trunc)

        ll_trunc = ll_trunc_func(omega_mis, 1)

        # statistical testing

        chi = 2 * (ll_joint - ll_trunc)
        pval_trunc = (self.tfd.Chi2(1.).survival_function(self.tf.cast(chi, self.tf.float32))).numpy()

        chi = 2 * (ll_joint - ll_mis)
        pval_mis = (self.tfd.Chi2(1.).survival_function(self.tf.cast(chi, self.tf.float32))).numpy()

        # Wald confidence intervals
        if omega_mis == 0.:
            lower_mis, upper_mis = 0., 0.
        else:
            se_log = np.sqrt(self.variance_log_omega(omega_mis, self.offset_mis / self.offset_syn))
            lower_mis, upper_mis = self.ci_log_wald(omega_mis, se_log, confidence_level=confidence_level)
        
        if omega_trunc == 0.:
            lower_trunc, upper_trunc = 0., 0.
        else:
            se_log = np.sqrt(self.variance_log_omega(omega_trunc, self.offset_trunc / self.offset_syn))
            lower_trunc, upper_trunc = self.ci_log_wald(omega_trunc, se_log, confidence_level=confidence_level)

        return omega_mis, lower_mis, upper_mis, pval_mis, np.nan, np.nan, omega_trunc, lower_trunc, upper_trunc, pval_trunc, np.nan, np.nan, t_hat_joint, t_hat_low, t_hat_high


    def log_like_map(self, omega, n_neutral, n_nonneutral, offset_neutral, offset_nonneutral):

        """
        Negative binomial model:
        - Overdispersion as in the NB prior across genes
        - Mean governed by MAP mean estimate -- not a marginal model
        """

        k = offset_neutral / self.offset_syn
        mu_map = self.optimal_mu(n_neutral, k)
        t_hat = max(1e-8, mu_map / self.offset_syn)
        mu = t_hat * offset_neutral + t_hat * offset_nonneutral * omega
        theta = self.theta
        if theta > 0:
            f = 1 / theta
            p = mu / (mu + f)
            model = self.tfd.NegativeBinomial(f, probs=p)
        else:
            model = self.tfd.Poisson(mu)
        return model.log_prob(n_neutral + n_nonneutral).numpy().item()
    

    def log_like_posterior_marginal(self, omega, n_neutral, n_nonneutral, offset_neutral, offset_nonneutral):

        """
        Negative binomial model:
        - Marginal over mutdensities drawn from gamma updated with gene specific syn mutdensity 
        - NB model with mean of the form (alpha/beta)*(ks + kn*omega) 

        """

        ks = 1
        kn = offset_nonneutral / offset_neutral
        alpha_update = n_neutral + self.alpha - 1 
        beta_update = 1 + self.beta
        k_eff = kn * omega + ks
        f = alpha_update
        p = k_eff / (beta_update + k_eff)
        model = self.tfd.NegativeBinomial(f, probs=p)
        return model.log_prob(n_neutral + n_nonneutral).numpy().item()

    
    def log_like_loc(self, omega, n_neutral, n_nonneutral, offset_neutral, offset_nonneutral):

        """
        - omega loc model
        - assumes n_syn is positive -- in case of globalloc, n_syn must be already adjusted
        """

        nonneutral_rate = n_neutral * (offset_nonneutral / offset_neutral) * omega
        mu = n_neutral + nonneutral_rate
        model = self.tfd.Poisson(rate=mu)
        return model.log_prob(n_neutral + n_nonneutral)


    def log_delta_method(self, mu, sigma2, nu, tau2):

        """
        Computes the variance of X/Y using the 1st order delta-method
        X ~ (mu, sigma2)
        Y ~ (nu, tau2)
        """

        b = sigma2 / mu ** 2 
        c = tau2 / nu ** 2

        return b + c

    
    def law_total_variance(self, omega, k):

        """Variance of counts distribution subject to omega and k"""

        return k * omega * ((k * omega * self.alpha / self.beta ** 2) + (self.alpha / self.beta))
        
    
    def variance_log_omega(self, omega, k):

        # omega_hat = nn / k * lambda_hat, where lambda_hat = (ns + alpha - 1) / (1 + beta)
        
        # x == nn
        mu = k * omega * self.alpha / self.beta
        sigma2 = self.law_total_variance(omega, k)

        # y == lambda_hat
        C = 1 / (self.beta + 1)
        D = (self.alpha - 1) / (self.beta + 1)
        nu = C * (self.alpha / self.beta) + D
        tau2 = C ** 2 * (self.alpha / self.beta ** 2) * (self.beta + 1)

        # delta method
        variance_log_delta = self.log_delta_method(mu, sigma2, nu, tau2)

        return variance_log_delta


    def ci_log_wald(self, estimate, se_log, confidence_level=0.95):

        a = 1 - confidence_level
        z = norm.ppf(1 - a/2)
        ci_lower = estimate * np.exp(-z * se_log)
        ci_upper = estimate * np.exp(z * se_log)
        
        return ci_lower, ci_upper


    def compute_lrt_pvalue(self, log_like_func, omega_hat, dof=1.0):
        """
        Computes robust LRT p-value.
        """

        ll_alt = self.tf.cast(log_like_func(omega_hat), self.tf.float64)
        ll_null = self.tf.cast(log_like_func(self.tf.constant(1.0, dtype=self.tf.float64)), self.tf.float64)

        if self.tf.math.is_nan(ll_alt) or self.tf.math.is_nan(ll_null):
            return 1.0  # Fallback to non-significant if likelihood is undefined

        # Compute test statistic and clamp at 0.0 to fix numerical precision artifacts
        lambda_stat = self.tf.maximum(0.0, 2.0 * (ll_alt - ll_null))

        # Compute p-value using Chi2 survival function
        chi2_dist = self.tfd.Chi2(df=self.tf.constant(dof, dtype=self.tf.float32))
        pval = chi2_dist.survival_function(self.tf.cast(lambda_stat, self.tf.float32))

        return pval.numpy().item()


    def fit_posterior_marginal(self, confidence_level=0.95):

        # synonymous mutation density

        mu_map = self.optimal_mu(self.n_syn, 1)
        t_hat = mu_map / self.offset_syn
        se_log = 1 / np.sqrt(self.n_syn + self.alpha - 1)  # SE(log t_hat) ~ (1/t_hat) * SE(t_hat), with SE(t_hat) via observed Fisher information
        t_hat_low, t_hat_high = self.ci_log_wald(t_hat, se_log)
        
        # missense
        # --------
        
        n = self.n_syn + self.n_mis
        ks = 1
        kn = self.offset_mis / self.offset_syn
        alpha_update = self.n_syn + self.alpha - 1 
        beta_update = self.beta + 1

        # MLE omega hat

        omega_mis = max(0, (n * (beta_update / alpha_update) - ks) / kn)

        def log_like_mis(omega):
            return self.log_like_posterior_marginal(omega, self.n_syn, self.n_mis, self.offset_syn, self.offset_mis)

        # pvalue

        pvalue_mis = self.compute_lrt_pvalue(log_like_mis, omega_mis, dof=1.0)

        # CI

        max_ll = log_like_mis(omega_mis)
        target_lrt = chi2.ppf(confidence_level, df=1) / 2.0
        target_ll = max_ll - target_lrt

        def profile_objective_mis(w):
            return log_like_mis(w) - target_ll

        lower_mis, upper_mis = get_bounds(profile_objective_mis, omega_mis)

        # p-value negative selection test
        pneg_mis = (omega_mis >= 1) * 1 + (omega_mis < 1) * pvalue_mis * 0.5
        
        # conditional cumulative probability P(n <= n_obs | omega=1) of observed under neutrality
        pcum_mis = 0
        for n in range(self.n_mis + 1):
            pcum_mis += np.exp(self.log_like_posterior_marginal(1., self.n_syn, n, self.offset_syn, self.offset_mis))


        # truncating
        # ----------

        n = self.n_syn + self.n_trunc
        ks = 1
        kn = self.offset_trunc / self.offset_syn
        alpha_update = self.n_syn + self.alpha - 1 
        beta_update = self.beta + 1

        # MLE omega hat

        omega_trunc = max(0, (n * (beta_update / alpha_update) - ks) / kn)

        def log_like_trunc(omega):
            return self.log_like_posterior_marginal(omega, self.n_syn, self.n_trunc, self.offset_syn, self.offset_trunc)

        # pvalue

        pvalue_trunc = self.compute_lrt_pvalue(log_like_trunc, omega_trunc, dof=1.0)

        # CI

        max_ll = log_like_trunc(omega_trunc)
        target_lrt = chi2.ppf(confidence_level, df=1) / 2.0
        target_ll = max_ll - target_lrt

        def profile_objective_trunc(w):
            return log_like_trunc(w) - target_ll

        lower_trunc, upper_trunc = get_bounds(profile_objective_trunc, omega_trunc)

        # p-value negative selection test
        pneg_trunc = (omega_trunc >= 1) * 1 + (omega_trunc < 1) * pvalue_trunc * 0.5
        
        # conditional cumulative probability P(n <= n_obs | omega=1) of observed under neutrality
        pcum_trunc = 0
        for n in range(self.n_mis + 1):
            pcum_trunc += np.exp(self.log_like_posterior_marginal(1., self.n_syn, n, self.offset_syn, self.offset_trunc))

        return omega_mis, lower_mis, upper_mis, pvalue_mis, pneg_mis, pcum_mis, omega_trunc, lower_trunc, upper_trunc, pvalue_trunc, pneg_trunc, pcum_trunc, t_hat, t_hat_low, t_hat_high


    
    def fit_map(self, confidence_level=0.95):
        
        # mean of MAP model if of the form t_hat * (ks + kn * omega)

        mu_map = self.optimal_mu(self.n_syn, 1)
        t_hat = mu_map / self.offset_syn
        se_log = 1 / np.sqrt(self.n_syn + self.alpha - 1)  # SE(log t_hat) ~ (1/t_hat) * SE(t_hat), with SE(t_hat) via observed Fisher information
        t_hat_low, t_hat_high = self.ci_log_wald(t_hat, se_log)
        
        # Missense
        # --------

        ks = self.offset_syn
        kn = self.offset_mis
        n = self.n_syn + self.n_mis

        def log_like_mis(omega):
            return self.log_like_map(omega, self.n_syn, self.n_mis, self.offset_syn, self.offset_mis)

        omega_mis = max(0, ((n / t_hat) - ks) / kn)
        
        # alternative versions
        # omega_mis = self.n_mis / (self.offset_mis * t_hat)
        # omega_mis_mle = float(minimize_scalar(lambda w: -float(log_like_mis(w)), bounds=(0.0, 100.0), method='bounded').x)

        # confidence interval

        max_ll = log_like_mis(omega_mis)
        target_lrt = chi2.ppf(confidence_level, df=1) / 2.0
        target_ll = max_ll - target_lrt

        def profile_objective_mis(w):
            return log_like_mis(w) - target_ll

        lower_mis, upper_mis = get_bounds(profile_objective_mis, omega_mis)

        # alternative version
        # se_log = np.sqrt(self.variance_log_omega(omega_mis, self.offset_mis / self.offset_syn))
        # lower_mis, upper_mis = self.ci_log_wald(omega_mis, se_log, confidence_level=confidence_level)

        # p-value
        pvalue_mis = self.compute_lrt_pvalue(log_like_mis, omega_mis, dof=1.0)
        
        # p-value negative selection test
        pneg_mis = (omega_mis >= 1) * 1 + (omega_mis < 1) * pvalue_mis * 0.5
        
        # conditional cumulative probability P(n <= n_obs | omega=1) of observed under neutrality
        pcum_mis = 0
        for n in range(self.n_mis + 1):
            pcum_mis += np.exp(self.log_like_map(1., self.n_syn, n, self.offset_syn, self.offset_mis))
        
        # Truncating
        # ----------

        ks = self.offset_syn
        kn = self.offset_trunc
        n = self.n_syn + self.n_trunc

        def log_like_trunc(omega):
            return self.log_like_map(omega, self.n_syn, self.n_trunc, self.offset_syn, self.offset_trunc)

        omega_trunc = max(0, ((n / t_hat) - ks) / kn)

        # alternative versions:
        # omega_trunc = self.n_trunc / (self.offset_trunc * t_hat)
        # omega_trunc_mle = float(minimize_scalar(lambda w: -float(log_like_trunc(w)), bounds=(0.0, 100.0), method='bounded').x)

        # confidence interval

        max_ll = log_like_trunc(omega_trunc)
        target_lrt = chi2.ppf(confidence_level, df=1) / 2.0
        target_ll = max_ll - target_lrt

        def profile_objective_trunc(w):
            return log_like_trunc(w) - target_ll

        lower_trunc, upper_trunc = get_bounds(profile_objective_trunc, omega_trunc)
        
        # alternative version
        # lower_trunc, upper_trunc = self.ci_wald(omega_trunc, self.offset_trunc / self.offset_syn, confidence_level=confidence_level)

        # p-value
        pvalue_trunc = self.compute_lrt_pvalue(log_like_trunc, omega_trunc, dof=1.0)

        # p-value negative selection test
        pneg_trunc = (omega_trunc >= 1) * 1 + (omega_trunc < 1) * pvalue_trunc * 0.5

        # conditional cumulative probability P(n <= n_obs | omega=1) of observed under neutrality
        pcum_trunc = 0
        for n in range(self.n_trunc + 1):
            pcum_trunc += np.exp(self.log_like_map(1., self.n_syn, n, self.offset_syn, self.offset_trunc))

        return omega_mis, lower_mis, upper_mis, pvalue_mis, pneg_mis, pcum_mis, omega_trunc, lower_trunc, upper_trunc, pvalue_trunc, pneg_trunc, pcum_trunc, t_hat, t_hat_low, t_hat_high


    def fit_loc(self, confidence_level=0.95):

        # synonymous mutation density

        t_hat = self.n_syn / self.offset_syn

        # missense
        # --------

        k = self.offset_mis / self.offset_syn
        n = self.n_syn + self.n_mis
        
        # MLE omega hat

        if self.n_syn > 0:
            omega_mis = max(0, ((n / self.n_syn) - 1) / k)
        else:
            omega_mis = 0
            # raise ValueError("syn count provided must be > 0")

        def log_like_mis(omega):
            return self.log_like_loc(omega, self.n_syn, self.n_mis, self.offset_syn, self.offset_mis)

        # pvalue

        pvalue_mis = self.compute_lrt_pvalue(log_like_mis, omega_mis, dof=1.0)

        # CI

        max_ll = log_like_mis(omega_mis)
        target_lrt = chi2.ppf(confidence_level, df=1) / 2.0
        target_ll = max_ll - target_lrt

        def profile_objective_mis(w):
            return log_like_mis(w) - target_ll

        lower_mis, upper_mis = get_bounds(profile_objective_mis, omega_mis)

        # p-value negative selection test
        pneg_mis = (omega_mis >= 1) * 1 + (omega_mis < 1) * pvalue_mis * 0.5
        
        # conditional cumulative probability P(n <= n_obs | omega=1) of observed under neutrality
        pcum_mis = 0
        for n in range(self.n_mis + 1):
            pcum_mis += np.exp(self.log_like_loc(1., self.n_syn, n, self.offset_syn, self.offset_mis))


        # truncating
        # ----------

        k = self.offset_trunc / self.offset_syn
        n = self.n_syn + self.n_trunc
        
        # MLE omega hat

        if self.n_syn > 0:
            omega_trunc = max(0, ((n / self.n_syn) - 1) / k)
        else:
            omega_trunc = 0
            # raise ValueError("syn count provided must be > 0")

        def log_like_trunc(omega):
            return self.log_like_loc(omega, self.n_syn, self.n_trunc, self.offset_syn, self.offset_trunc)

        # pvalue

        pvalue_trunc = self.compute_lrt_pvalue(log_like_trunc, omega_trunc, dof=1.0)

        # CI

        max_ll = log_like_trunc(omega_trunc)
        target_lrt = chi2.ppf(confidence_level, df=1) / 2.0
        target_ll = max_ll - target_lrt

        def profile_objective_trunc(w):
            return log_like_trunc(w) - target_ll

        lower_trunc, upper_trunc = get_bounds(profile_objective_trunc, omega_trunc)

        # p-value negative selection test
        pneg_trunc = (omega_trunc >= 1) * 1 + (omega_trunc < 1) * pvalue_trunc * 0.5
        
        # conditional cumulative probability P(n <= n_obs | omega=1) of observed under neutrality
        pcum_trunc = 0
        for n in range(self.n_mis + 1):
            pcum_trunc += np.exp(self.log_like_loc(1., self.n_syn, n, self.offset_syn, self.offset_trunc))

        return omega_mis, lower_mis, upper_mis, pvalue_mis, pneg_mis, pcum_mis, omega_trunc, lower_trunc, upper_trunc, pvalue_trunc, pneg_trunc, pcum_trunc, t_hat, np.nan, np.nan



    def get_fitting_methods(self):
        """Finds all callable methods starting with 'fit_'."""
        return [
            getattr(self, attr)
            for attr in dir(self)
            if attr.startswith("fit_") and callable(getattr(self, attr))
        ]