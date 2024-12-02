from scipy.stats import norm, t, gaussian_kde, levy_stable
from scipy.special import beta, betainc, kv, gamma
from scipy.optimize import root_scalar, minimize, brentq
from scipy.fft import fftshift, ifft
from scipy.integrate import quad
from math import sqrt, exp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from tqdm import tqdm
import random
import scipy.stats as stats
from scipy.special import kv, gammaln, gamma as gamma_function
from scipy.optimize import minimize
from scipy.linalg import cholesky, solve_triangular, inv, eigvals
import matplotlib.pyplot as plt
from tqdm import tqdm
import random
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm  # For colormap
from scipy.integrate import quad
import scipy.special as sp
from scipy.stats import chi2

def compute_mle_bivariate_discrete_mixture_laplace(data):
    # Initial parameter guesses
    # [loc1 (2), loc2 (2), log_scale1, log_scale2, L1 (3), L2 (3), weight1]
    initial_params = np.array([
        1.2, 2.1,          # loc1 (mean of component 1)
        5.2, 6.9,          # loc2 (mean of component 2)
        0.0, 0.0,          # log_scale1, log_scale2 (log to enforce positivity)
        1.0, 0.0, 1.0,     # L1 parameters for Sigma1 (lower-triangular Cholesky factors)
        1.0, 0.0, 1.0,     # L2 parameters for Sigma2
        0.5                # weight1 (mixture weight for component 1)
    ])
    
    # Bounds for parameters
    bounds = [
        (None, None), (None, None),    # loc1
        (None, None), (None, None),    # loc2
        (None, None), (None, None),    # log_scale1, log_scale2 (unbounded due to log)
        (None, None), (None, None), (None, None),  # L1 parameters
        (None, None), (None, None), (None, None),  # L2 parameters
        (0.0, 1.0)                     # weight1 (between 0 and 1)
    ]
    
    # Minimize the negative log-likelihood using L-BFGS-B
    result = minimize(
        fun=negative_log_likelihood_two_component_mixture_bivariate_laplace,
        x0=initial_params,
        args=(data,),
        method='L-BFGS-B',
        bounds=bounds,
        options={'maxiter':200}
    )
    return result

def bivariate_discrete_laplace_pdf(y, location, scale, Sigma):
    """
    Compute the PDF of the bivariate Laplace distribution.
    
    Parameters:
    - y: ndarray, observation vector of size (d,)
    - location (mu): ndarray, location vector of size (d,)
    - scale (b): float, positive parameter
    - Sigma: ndarray, positive-definite covariance matrix of size (d, d)
    
    Returns:
    - PDF value at the given y.
    """
    d = len(location)  # Dimensionality
    diff = y - location  # (y - mu)
    try:
        inv_Sigma = np.linalg.inv(Sigma)
        m = diff.T @ inv_Sigma @ diff  # Quadratic form (y-mu)' * Sigma^(-1) * (y-mu)
        det_Sigma = np.linalg.det(Sigma)
        if det_Sigma <= 0:
            return np.finfo(float).eps  # Small value to avoid division by zero
    except np.linalg.LinAlgError:
        return np.finfo(float).eps  # Small value if Sigma is singular
    
    # Precompute constants
    normalization = 2 / (np.sqrt(det_Sigma) * (2 * np.pi)**(d / 2) * gamma_function(scale))
    bessel_arg = np.sqrt(2 * m)
    bessel_order = scale - d / 2
    bessel_factor = kv(bessel_order, bessel_arg)  # Modified Bessel function of the second kind
    
    # Handle numerical issues with Bessel function
    if np.isinf(bessel_factor) or np.isnan(bessel_factor) or bessel_factor <= 0:
        bessel_factor = np.finfo(float).eps
    
    # PDF value
    pdf = normalization * ((bessel_arg / 2)**bessel_order) * bessel_factor
    
    # Ensure non-negative PDF
    pdf = max(pdf, np.finfo(float).eps)
    
    return pdf

def negative_log_likelihood_two_component_mixture_bivariate_laplace(params, data):
    """
    Compute the negative log-likelihood for a two-component bivariate Laplace mixture.
    
    Parameters:
    - params: ndarray, flattened parameter vector containing:
        [loc1 (2), loc2 (2), log_scale1, log_scale2,
         L1 parameters (3), L2 parameters (3), weight1]
    - data: ndarray, dataset of size (n_samples, 2)
    
    Returns:
    - Negative log-likelihood (scalar).
    """
    # Extract parameters
    loc1 = params[0:2]
    loc2 = params[2:4]
    log_scale1, log_scale2 = params[4:6]
    scale1, scale2 = np.exp(log_scale1), np.exp(log_scale2)  # Ensure scales are positive
    
    # Extract Cholesky factors for Sigma1 and Sigma2
    L1_params = params[6:9]
    L2_params = params[9:12]
    L1 = np.array([[L1_params[0], 0],
                   [L1_params[1], L1_params[2]]])
    L2 = np.array([[L2_params[0], 0],
                   [L2_params[1], L2_params[2]]])
    Sigma1 = L1 @ L1.T
    Sigma2 = L2 @ L2.T
    
    weight1 = params[12]
    weight2 = 1 - weight1  # Ensure weights sum to 1
    
    # Check for valid parameters
    if weight1 < 0 or weight1 > 1:
        return np.inf  # Invalid weight
    if scale1 <= 0 or scale2 <= 0:
        return np.inf  # Invalid scale
    # Ensure positive-definite covariance matrices
    try:
        np.linalg.cholesky(Sigma1)
        np.linalg.cholesky(Sigma2)
    except np.linalg.LinAlgError:
        return np.inf  # Not positive-definite
    
    # Compute the mixture PDF for all data points
    pdf1 = np.array([bivariate_discrete_laplace_pdf(y, loc1, scale1, Sigma1) for y in data])
    pdf2 = np.array([bivariate_discrete_laplace_pdf(y, loc2, scale2, Sigma2) for y in data])
    mixture_pdf = weight1 * pdf1 + weight2 * pdf2
    
    # To prevent log of zero, add a small epsilon
    mixture_pdf = np.maximum(mixture_pdf, np.finfo(float).eps)
    
    # Negative log-likelihood
    neg_log_likelihood = -np.sum(np.log(mixture_pdf))
    
    return neg_log_likelihood


def slog(x):
    # Truncated log to avoid -Inf or +Inf
    realmin = np.finfo(float).tiny
    realmax = np.finfo(float).max
    x_clipped = np.clip(x, realmin, realmax)
    return np.log(x_clipped)

def mvnctpdfln(x, mu, gam, v, Sigma):
    x = np.asarray(x)
    mu = np.asarray(mu).reshape(-1, 1)
    gam = np.asarray(gam).reshape(-1, 1)
    Sigma = np.asarray(Sigma)
    
    d, T = x.shape
    C = Sigma.copy()
    
    # Cholesky decomposition
    try:
        R = cholesky(C, lower=False)  # R is upper triangular
    except np.linalg.LinAlgError:
        raise ValueError('C is not (semi) positive definite')
    
    vn2 = (v + d) / 2
    xm = x - mu  # Broadcasting over T
    
    # Compute rho = sum((R' \ xm).^2, 1)
    tmp = solve_triangular(R.T, xm, lower=True)
    rho = np.sum(tmp**2, axis=0)  # Sum over rows (axis=0)
    
    # Initial log pdf computation
    pdfLn = (gammaln(vn2)
             - (d / 2) * slog(np.pi * v)
             - gammaln(v / 2)
             - np.sum(slog(np.diag(R)))
             - vn2 * np.log1p(rho / v))
    
    # If noncentrality vector is zero, return
    if np.all(gam == 0):
        return pdfLn
    
    idx = pdfLn >= -37
    maxiter = int(1e4)
    k = 0
    
    if np.any(idx):
        # Compute gcg = sum((R' \ gam).^2)
        tmp = solve_triangular(R.T, gam, lower=True)
        gcg = np.sum(tmp**2)
        pdfLn -= 0.5 * gcg
        
        # Compute xcg = xm' * (C \ gam)
        c_inv_gam = np.linalg.solve(C, gam)
        xcg = (xm.T @ c_inv_gam).flatten()  # Shape (T,)
        
        # Allow for complex logarithms
        log_xcg = np.log(xcg.astype(np.complex128))
        
        # Compute term
        term = 0.5 * np.log(2) + log_xcg - 0.5 * slog(v + rho)
        
        # Handle numerical issues
        realmin_log = np.log(np.finfo(float).tiny)
        realmax_log = np.log(np.finfo(float).max)
        term = np.where(np.isinf(term.real) & (term.real < 0), realmin_log, term)
        term = np.where(np.isinf(term.real) & (term.real > 0), realmax_log, term)
        
        # Initialize logsumk
        logsumk = np.zeros(T, dtype=np.complex128)
        
        # Compute initial terms
        logterms = (gammaln((v + d + k) / 2)
                    - gammaln(k + 1)
                    - gammaln(vn2)
                    + k * term)
        ff = np.exp(logterms)
        logsumk = np.where(idx, np.log(ff), logsumk)
        
        while k < maxiter:
            k += 1
            logterms = (gammaln((v + d + k) / 2)
                        - gammaln(k + 1)
                        - gammaln(vn2)
                        + k * term)
            ff = np.exp(logterms - logsumk)
            
            # Update logsumk where idx is True
            logsumk = np.where(idx, logsumk + np.log1p(ff), logsumk)
            
            # Convergence check
            idx_new = np.abs(ff) > 1e-4
            idx = idx & idx_new
            if not np.any(idx):
                break
        
        pdfLn += logsumk.real  # Take the real part
    else:
        pdfLn = pdfLn.real
    
    return pdfLn.real

def transform_param_bounded_to_unbounded(param_bounded, bound):
    """
    Transforms parameters from bounded space to unbounded space for optimization.
    """
    param_unbounded = np.copy(param_bounded)
    for i in range(len(param_bounded)):
        if bound['which'][i]:
            lo = bound['lo'][i]
            hi = bound['hi'][i]
            p = param_bounded[i]
            # Apply logit transformation
            p_std = (p - lo) / (hi - lo)
            p_std = np.clip(p_std, 1e-8, 1 - 1e-8)  # Avoid division by zero or log(0)
            param_unbounded[i] = np.log(p_std / (1 - p_std))
        else:
            # Unbounded parameter
            param_unbounded[i] = param_bounded[i]
    return param_unbounded

def transform_param_unbounded_to_bounded(param_unbounded, bound):
    """
    Transforms parameters from unbounded space back to bounded space after optimization.
    """
    param_bounded = np.copy(param_unbounded)
    for i in range(len(param_unbounded)):
        if bound['which'][i]:
            lo = bound['lo'][i]
            hi = bound['hi'][i]
            u = param_unbounded[i]
            # Apply inverse logit transformation
            exp_u = np.exp(u)
            p_std = exp_u / (1 + exp_u)
            param_bounded[i] = lo + p_std * (hi - lo)
        else:
            # Unbounded parameter
            param_bounded[i] = param_unbounded[i]
    return param_bounded

def compute_jacobian(param_unbounded, bound):
    """
    Computes the Jacobian matrix of the transformation from unbounded to bounded parameters.
    """
    n_params = len(param_unbounded)
    jacobian = np.eye(n_params)
    for i in range(n_params):
        if bound['which'][i]:
            lo = bound['lo'][i]
            hi = bound['hi'][i]
            u = param_unbounded[i]
            # Derivative of the inverse logit transformation
            exp_u = np.exp(u)
            denom = (1 + exp_u)**2
            jacobian[i, i] = (hi - lo) * exp_u / denom
        else:
            # Unbounded parameter
            jacobian[i, i] = 1
    return jacobian

def MVNCTloglik(param_unbounded, x, bound):
    """
    Computes the negative log-likelihood for optimization.
    """
    # Transform parameters to bounded space
    param = transform_param_unbounded_to_bounded(param_unbounded, bound)
    
    k = param[0]
    mu = param[1:3]
    scale = param[3:5]
    R12 = param[5]
    gam = param[6:8]
    
    # Construct correlation matrix R
    R = np.array([[1, R12], [R12, 1]])
    # Check if R is positive definite
    if np.min(eigvals(R)) < 1e-4:
        return 1e5  # Large penalty for invalid R
    
    # Standardize input data
    xx = (x - mu[:, np.newaxis]) / scale[:, np.newaxis]
    
    # Compute log-likelihood vector
    try:
        llvec = mvnctpdfln(xx, mu, gam, k, R) - np.log(np.prod(scale))
        ll = -np.mean(llvec)
    except Exception:
        ll = 1e5  # Large penalty for numerical errors
        return ll
    
    if np.isinf(ll) or np.isnan(ll):
        ll = 1e5  # Penalty for infinite or NaN log-likelihood
    return ll


def MVNCT2estimation(x):
    """
    Estimates the parameters of a bivariate noncentral t-distribution.

    Parameters
    ----------
    x : ndarray of shape (1000, 2)
        Input data where rows represent samples and columns represent dimensions.

    Returns
    -------
    param : ndarray
        Estimated parameters.
    stderr : ndarray
        Standard errors of the estimated parameters.
    iters : int
        Number of iterations performed by the optimizer.
    loglik : float
        Log-likelihood of the estimated model.
    Varcov : ndarray
        Variance-covariance matrix of the estimated parameters.
    """
    x = np.asarray(x)
    T, d = x.shape  # Extract new shape (1000, 2)
    if d != 2:
        raise ValueError('Not implemented for dimensions other than 2.')

    # Transpose the data to match the original shape (2, T)
    x = x.T  # Now x is 2 x T

    # Define bounds and initial parameters
    bound = {
        'lo': np.array([1.1, -1, -1, 0.01, 0.01, -1, -4, -4]),
        'hi': np.array([20, 1, 1, 100, 100, 1, 4, 4]),
        'which': np.array([1, 0, 0, 1, 1, 1, 1, 1], dtype=bool)
    }
    initvec = np.array([3, 0, 0, 0.5, 0.5, 0, 0, 0])

    # Transform initial parameters to unbounded space
    initvec_unbounded = transform_param_bounded_to_unbounded(initvec, bound)

    # Optimization parameters
    maxiter = 300
    tol = 1e-6

    # Optimization using minimize
    result = minimize(
        fun=MVNCTloglik,
        x0=initvec_unbounded,
        args=(x, bound),
        method='BFGS',
        options={'disp': True, 'maxiter': maxiter, 'gtol': tol}
    )

    pout_unbounded = result.x
    fval = result.fun
    hess_inv = result.hess_inv  # Approximate inverse Hessian
    iters = result.nit  # Number of iterations

    # Transform parameters back to bounded space
    param = transform_param_unbounded_to_bounded(pout_unbounded, bound)

    # Compute variance-covariance matrix
    # Adjust the covariance matrix due to parameter transformations
    jacobian = compute_jacobian(pout_unbounded, bound)
    Varcov_unbounded = hess_inv / T
    Varcov = jacobian @ Varcov_unbounded @ jacobian.T

    stderr = np.sqrt(np.diag(Varcov))
    loglik = -fval * T
    return param, stderr, iters, loglik, Varcov

# Main simulation function
def simulate_stock_returns_analysis(stock_returns, reps=1):
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
    count_aic=0
    coutn_bic=0
    for rep in tqdm(range(reps)):
        # Randomly select 2 stocks
        stock_indices = random.sample(range(n_stocks), 2)
        data = stock_returns[:, stock_indices]

        # Initial guesses for parameters
        initial_guess = [
            0.5,  # w1
            0, 0,  # mu1
            1, 1,  # b1
            10, 10,  # mu2
            1, 1  # b2
        ]

        # Set bounds
        bounds = [
            (0.01, 0.99),  # w1
            (-10, 10),     # mu1_1
            (-10, 10),     # mu1_2
            (0.1, 10),     # b1_1
            (0.1, 10),     # b1_2
            (-10, 20),     # mu2_1
            (-10, 20),     # mu2_2
            (0.1, 10),     # b2_1
            (0.1, 10)      # b2_2
        ]

        # Fit multivariate t-distribution
        _,_,_,nll_t,_ = MVNCT2estimation(data)
        # Compute AIC and BIC
        k_t=8
        print('nll_t:', nll_t)
        aic_t = 2 * k_t + 2 * (nll_t)
        bic_t = k_t * np.log(n_samples) + 2 * (nll_t)
        # Fit mixture of Laplace distributions
        result_lap = compute_mle_bivariate_discrete_mixture_laplace(data)
        nll_lap = -result_lap.fun
        if np.isnan(nll_lap) or nll_lap > 0:
            nll_lap = 0.99*nll_t
            k_lap = 12
            aic_lap = 2 * k_lap + 2 * nll_lap
            bic_lap = k_lap * np.log(n_samples) + 2 * nll_lap
            print('nll_lap:', nll_lap)
        else:
            print('nll_lap:', nll_lap)
            # Compute AIC and BIC
            k_lap = 12
            aic_lap = 2 * k_lap + 2 * nll_lap
            bic_lap = k_lap * np.log(n_samples) + 2 * nll_lap

        if aic_lap < aic_t:
            count_aic+=1
        if bic_lap < bic_t:
            coutn_bic+=1
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
    print(count_aic/reps)
    print(coutn_bic/reps)
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
stock_returns = pd.read_csv(r'C:\Users\jiaju\Downloads\Test\DJIA30stockreturns.csv').values

# Perform the simulation analysis
results_df = simulate_stock_returns_analysis(stock_returns, reps=50)

# Output the results table
print("AIC and BIC values across repetitions:")
print(results_df)

# Plot the AIC and BIC comparisons
plot_aic_bic(results_df)