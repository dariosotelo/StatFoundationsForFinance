import numpy as np
import pandas as pd
from scipy.special import kv, gammaln
from scipy.optimize import minimize
from numpy.linalg import inv, det, slogdet
import matplotlib.pyplot as plt
import random

# Function to fit a mixture of two bivariate Laplace distributions
def fit_mixture_laplace(data):
    """
    Fit a mixture of two bivariate Laplace distributions to the data.
    Returns the negative log-likelihood at the optimum, the number of parameters (k),
    and the optimized parameters.
    """
    n_samples = data.shape[0]
    
    def sigmoid(gamma):
        return 1 / (1 + np.exp(-gamma))
    
    def construct_covariance_matrix(params):
        a, b, c = params
        # Ensure positive definiteness
        a = np.exp(a)
        c = np.exp(c)
        Sigma = np.array([[a, b], [b, c]])
        if det(Sigma) <= 0:
            Sigma += np.eye(2) * 1e-6
        return Sigma
    
    def bivariate_laplace_pdf(x, mu, Sigma, lam=1.0):
        d = len(x)
        delta = x - mu
        Q = delta.T @ inv(Sigma) @ delta
        sqrt_Q = np.sqrt(Q)
        det_Sigma = det(Sigma)
        if det_Sigma <= 0 or sqrt_Q == 0:
            return 1e-10  # Small value to avoid log(0)
        prefactor = (lam ** d) / ((2 * np.pi) ** (d / 2) * det_Sigma ** 0.5)
        bessel_term = kv((d - 2) / 2, lam * sqrt_Q)
        if bessel_term == 0 or np.isnan(bessel_term):
            return 1e-10
        pdf = prefactor * (bessel_term / (sqrt_Q ** ((d - 2) / 2)))
        return pdf
    
    def mixture_pdf(x, params):
        pi = params[0]
        mu1 = np.array(params[1:3])
        mu2 = np.array(params[3:5])
        Sigma1_params = params[5:8]
        Sigma2_params = params[8:11]
        Sigma1 = construct_covariance_matrix(Sigma1_params)
        Sigma2 = construct_covariance_matrix(Sigma2_params)
        lam1 = params[11]
        lam2 = params[12]
        f1 = bivariate_laplace_pdf(x, mu1, Sigma1, lam1)
        f2 = bivariate_laplace_pdf(x, mu2, Sigma2, lam2)
        pdf = pi * f1 + (1 - pi) * f2
        return pdf
    
    def negative_log_likelihood(params):
        gamma = params[0]
        pi = sigmoid(gamma)
        mu1 = params[1:3]
        mu2 = params[3:5]
        Sigma1_params = params[5:8]
        Sigma2_params = params[8:11]
        lam1 = np.exp(params[11])
        lam2 = np.exp(params[12])
        nll = 0
        for x in data:
            pdf_value = mixture_pdf(x, [pi, *mu1, *mu2, *Sigma1_params, *Sigma2_params, lam1, lam2])
            if pdf_value <= 0 or np.isnan(pdf_value):
                nll += 1e6
            else:
                nll -= np.log(pdf_value)
        return nll
    
    # Initial parameter guesses
    gamma_init = 0
    mu1_init = np.mean(data, axis=0) - 1
    mu2_init = np.mean(data, axis=0) + 1
    Sigma1_params_init = [0, 0, 0]
    Sigma2_params_init = [0, 0, 0]
    lam1_init = 0
    lam2_init = 0
    initial_params = np.hstack((gamma_init, mu1_init, mu2_init, Sigma1_params_init, Sigma2_params_init, lam1_init, lam2_init))
    
    # Number of parameters: gamma, 2*mu, 2*3 Sigma params, 2*lambda
    k = 1 + 2 + 2 + 3 + 3 + 1 + 1  # Total 13 parameters
    
    result = minimize(negative_log_likelihood, initial_params, method='BFGS', options={'disp': False})
    nll_opt = result.fun
    optimized_params = result.x
    return nll_opt, k, optimized_params

# Function to fit a multivariate t-distribution
def fit_multivariate_t(data):
    """
    Fit a multivariate t-distribution to the data.
    Returns the negative log-likelihood at the optimum, the number of parameters (k),
    and the optimized parameters.
    """
    n_samples, d = data.shape
    
    def mvt_log_likelihood(params):
        nu = params[0]  # Degrees of freedom
        mu = params[1:1+d]
        Sigma_params = params[1+d:]
        # Construct Sigma
        L = np.zeros((d, d))
        idx = np.tril_indices(d)
        L[idx] = Sigma_params
        Sigma = L @ L.T  # Sigma = L * L'
        sign, logdet = slogdet(Sigma)
        if sign <= 0 or nu <= 0:
            return 1e6
        delta = data - mu
        Q = np.sum(delta @ inv(Sigma) * delta, axis=1)
        ll = gammaln((nu + d) / 2) - gammaln(nu / 2) - (d / 2) * np.log(nu * np.pi) - 0.5 * logdet
        ll -= ((nu + d) / 2) * np.log(1 + Q / nu)
        nll = -np.sum(ll)
        if np.isnan(nll):
            return 1e6
        return nll
    
    # Initial parameter guesses
    nu_init = 5.0
    mu_init = np.mean(data, axis=0)
    # Cholesky decomposition of covariance matrix
    cov_init = np.cov(data, rowvar=False)
    L_init = np.linalg.cholesky(cov_init)
    Sigma_params_init = L_init[np.tril_indices(d)]
    initial_params = np.hstack((nu_init, mu_init, Sigma_params_init))
    
    # Number of parameters: nu, mu (d), Sigma (d*(d+1)/2)
    k = 1 + d + d * (d + 1) // 2  # Total parameters
    
    result = minimize(mvt_log_likelihood, initial_params, method='BFGS', options={'disp': False})
    nll_opt = result.fun
    optimized_params = result.x
    return nll_opt, k, optimized_params

# Main simulation function
def simulate_stock_returns_analysis(stock_returns, reps=50):
    """
    Performs the simulation analysis on stock returns data.
    Args:
        stock_returns: DataFrame or ndarray of stock returns (rows: observations, columns: stocks)
        reps: Number of repetitions
    Returns:
        results_df: DataFrame containing AIC and BIC values for each model across repetitions
    """
    n_stocks = stock_returns.shape[1]
    n_samples = stock_returns.shape[0]
    results = []

    for rep in range(reps):
        # Randomly select 2 stocks
        stock_indices = random.sample(range(n_stocks), 2)
        data = stock_returns[:, stock_indices]

        # Fit mixture of Laplace distributions
        nll_lap, k_lap, _ = fit_mixture_laplace(data)
        # Compute AIC and BIC
        aic_lap = 2 * k_lap + 2 * nll_lap
        bic_lap = k_lap * np.log(n_samples) + 2 * nll_lap

        # Fit multivariate t-distribution
        nll_t, k_t, _ = fit_multivariate_t(data)
        # Compute AIC and BIC
        aic_t = 2 * k_t + 2 * nll_t
        bic_t = k_t * np.log(n_samples) + 2 * nll_t

        # Store results
        results.append({
            'Rep': rep + 1,
            'AIC_Laplace': aic_lap,
            'BIC_Laplace': bic_lap,
            'AIC_t': aic_t,
            'BIC_t': bic_t
        })

        print(f"Iteration {rep + 1}/{reps} completed.")

    results_df = pd.DataFrame(results)
    return results_df

# Plotting function
def plot_aic_bic(results_df):
    """
    Generates plots comparing the AIC and BIC values for both models.
    """
    reps = results_df['Rep']
    plt.figure(figsize=(12, 6))
    plt.plot(reps, results_df['AIC_Laplace'], label='AIC Laplace', marker='o')
    plt.plot(reps, results_df['AIC_t'], label='AIC t-Distribution', marker='x')
    plt.title('Comparison of AIC Values')
    plt.xlabel('Repetition')
    plt.ylabel('AIC')
    plt.legend()
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(12, 6))
    plt.plot(reps, results_df['BIC_Laplace'], label='BIC Laplace', marker='o')
    plt.plot(reps, results_df['BIC_t'], label='BIC t-Distribution', marker='x')
    plt.title('Comparison of BIC Values')
    plt.xlabel('Repetition')
    plt.ylabel('BIC')
    plt.legend()
    plt.grid(True)
    plt.show()


# Load your actual stock return data here
# For demonstration, we'll generate synthetic data
stock_returns = pd.read_csv('DJIA30stockreturns.csv').values

# Perform the simulation analysis
results_df = simulate_stock_returns_analysis(stock_returns, reps=50)

# Output the results table
print("AIC and BIC values across repetitions:")
print(results_df)

# Plot the AIC and BIC comparisons
plot_aic_bic(results_df)
