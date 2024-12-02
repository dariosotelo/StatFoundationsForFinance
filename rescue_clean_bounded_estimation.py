#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 11:52:00 2024

@author: darios
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 14:55:03 2024

@author: 
"""

# Libraries
import random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm  # For colormap
from scipy.integrate import quad
import scipy.special as sp
from scipy.optimize import minimize
from scipy.special import gammaln, gamma, kv
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import minimize
from scipy.linalg import inv, eigvals
from scipy.stats import chi2
from tqdm import tqdm
import pandas as pd

#%% II.1

def bivariate_laplace(mu, b, n):
    x1 = np.random.laplace(loc=mu[0], scale=b[0], size=n)
    x2 = np.random.laplace(loc=mu[1], scale=b[1], size=n)
    return np.column_stack((x1, x2))

def simulate_mixture_bivariate_laplace(pi, mu1, b1, Sigma1, mu2, b2, Sigma2, n):
    n1 = int(pi * n)
    n2 = n - n1
    
    samples1 = bivariate_laplace(mu1, b1, n1)
    samples1 = samples1 @ np.linalg.cholesky(Sigma1).T
    
    samples2 = bivariate_laplace(mu2, b2, n2)
    samples2 = samples2 @ np.linalg.cholesky(Sigma2).T

    return np.vstack((samples1, samples2))

# Parameters
pi = 0.7
mu1 = [20, 40]
b1 = [10, 10]
Sigma1 = np.array([[1, 0.5], [0.5, 1]])

mu2 = [0, 0]
b2 = [5, 5]
Sigma2 = np.array([[4, 2], [2, 4]])

n = 1000

# Simulate the mixture
samples = simulate_mixture_bivariate_laplace(pi, mu1, b1, Sigma1, mu2, b2, Sigma2, n)

# 2D Histogram (binning the data)
bins = 20  # Number of bins in each dimension
hist, xedges, yedges = np.histogram2d(samples[:, 0], samples[:, 1], bins=bins)

# Coordinates of bin centers
xpos, ypos = np.meshgrid(xedges[:-1] + np.diff(xedges) / 2, yedges[:-1] + np.diff(yedges) / 2, indexing="ij")
xpos = xpos.ravel()
ypos = ypos.ravel()
zpos = np.zeros_like(xpos)

# Heights of the bars (histogram counts)
heights = hist.ravel()

# Dimensions of the bars
dx = dy = np.diff(xedges)[0]  # Width of each bin
dz = heights  # Heights of the bars

# Normalize heights for colormap
norm = plt.Normalize(vmin=heights.min(), vmax=heights.max())
colors = cm.viridis(norm(heights))  # Use the "viridis" colormap

# Create the 3D bar plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, alpha=0.7)

# Set labels and title
ax.set_title("3D Bar Plot of 2-Component Bivariate Mixture of Laplace")
ax.set_xlabel("X1")
ax.set_ylabel("X2")
ax.set_zlabel("Frequency")

# Add a colorbar for density
sm = cm.ScalarMappable(cmap=cm.viridis, norm=norm)
sm.set_array([])  # Dummy array for the colorbar
plt.colorbar(sm, ax=ax, shrink=0.5, aspect=10, label="Density")

plt.show()

#%% II.2


# Part a

# Our definition of the bessel_k
def bessel_k_int(nu, x):
    integrand = lambda u: 0.5*u**(nu-1)*np.exp(-x/2*(1/u+u))
    result, _ = quad(integrand, 0, np.inf)
    return result

# There is no constraint in nu, feel free to change it as you wish
nu = 2
x = 2.0

# Python's definition of the bessel_k
sp_result = sp.kv(nu, x)
our_result = bessel_k_int(nu, x)

# This function only graphs the previous code over a range of values.
def graph_difference(nu, x_values):
    y_values = [abs(bessel_k_int(nu, x_val) - sp.kv(nu, x_val)) for x_val in x_values]
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, y_values, marker='o', linestyle='-', color='blue', label='Difference')
    plt.title('Difference between bessel_k_int and scipy bessel_k')
    plt.xlabel('x values')
    plt.ylabel('Absolute Difference')
    plt.legend()
    plt.grid(True)
    plt.show()


# You can change the range, after 20 the difference is minimal
x_values = np.linspace(-15, 15, 200)
    
graph_difference(nu, x_values)

#%%
# Part b
# MLE of the k=2 component, d=2 (bivariate) discrete mixture of Laplace.
def compute_mle_bivariate_discrete_mixture_laplace(data):
    #initial_params = np.random.rand(13)  # Locations, scales, Sigma terms, weight
    initial_params = np.array([
        1.2, 2.1,          # loc1 (mean of component 1)
        5.2, 6.9,          # loc2 (mean of component 2)
        1.0, 0.9,      # scale1, scale2
        1.0, 0.02, 1.0, # Sigma1 (3 values: [var1, cov, var2])
        1.4, -0.2, 0.6#,# Sigma2 (3 values: [var1, cov, var2])
        #0.5            # weight1 (mixture weight for component 1)
    ])
    bounds_x0 = [(0.9,2.1),(1.5,2.5),(4.5,5.5),(5.5,7.5),(0.8,1.2),(0.8,1.2),(0.8,1.2),(0.019,0.6),(0.8,1.3),(1.3,1.9),(-0.9,-0.1),(0.5,1)]
    #            loc1       loc2        loc1     loc2      scale1  scale2    var1      cov           var2
    # Minimize the negative log-likelihood using BFGS
    result = minimize(
        fun = negative_log_likelihood_two_component_mixture_bivariate_laplace,
        x0=initial_params,
        args=(data,),
        method='L-BFGS-B',
        bounds = bounds_x0
    )
    return result


def bivariate_discrete_laplace_pdf(location, scale, Sigma, y):
    """
    Compute the PDF of the bivariate Laplace distribution as given in the formula.
    
    I got this out of Paolella's book

    Parameters:
    - y: ndarray, observation vector of size (d,)
    - location (mu): ndarray, location vector of size (d,)
    - Sigma: ndarray, positive-definite covariance matrix of size (d, d)
    - scale (b): float, parameter of the gamma distribution (must be > 0)

    Returns:
    - PDF value at the given y.
    """
    d = len(location)  # Dimensionality
    diff = y - location  # (y - mu)
    m = np.dot(diff.T, np.linalg.inv(Sigma)).dot(diff)  # Quadratic form (y-mu)' * Sigma^(-1) * (y-mu)

    # Determinant of Sigma
    det_Sigma = np.linalg.det(Sigma)

    # Precompute constants
    normalization = 2 / (np.sqrt(det_Sigma) * (2 * np.pi)**(d / 2) * gamma(scale))
    bessel_arg = np.sqrt(2 * m)
    bessel_factor = kv(scale - d / 2, bessel_arg)  # Modified Bessel function of the second kind

    # PDF value
    pdf = normalization * ((bessel_arg/2)**(scale/2 - d/4)) * bessel_factor

    return pdf



def negative_log_likelihood_two_component_mixture_bivariate_laplace(params, data):
    """
    Compute the negative log-likelihood for a two-component bivariate Laplace mixture.

    Parameters:
    - params: ndarray, flattened parameter vector containing:
        [location_1 (2 values), location_2 (2 values),
         scale_1, scale_2,
         Sigma_1 (4 values for 2x2 matrix), Sigma_2 (4 values for 2x2 matrix),
         weight_1 (mixture weight for first component)]
    - data: ndarray, dataset of size (n_samples, 2)

    Returns:
    - Negative log-likelihood (scalar).
    """
    # Extract parameters
    loc1 = params[:2]  # Location for component 1
    loc2 = params[2:4]  # Location for component 2
    scale1, scale2 = params[4:6]  # Scales for components 1 and 2
    
    Sigma1 = np.array([[params[6], params[7]], [params[7], params[8]]])  # Sigma_1
    Sigma2 = np.array([[params[9], params[10]], [params[10], params[11]]])  # Sigma_2
    
    #weight1 = params[12]  # Mixture weight for component 1
    weight1 = 0.6
    weight2 = 1 - weight1  # Mixture weight for component 2 (ensures weights sum to 1)
    
    # Compute log-likelihood for each data point
    log_likelihoods = []
    for y in data:
        # Mixture model: weighted sum of the two PDFs
        pdf1 = weight1 * bivariate_discrete_laplace_pdf(loc1, scale1, Sigma1, y)
        pdf2 = weight2 * bivariate_discrete_laplace_pdf(loc2, scale2, Sigma2, y)
        mixture_pdf = pdf1 + pdf2
        
        # Logarithm of the PDF for numerical stability
        if mixture_pdf > 0:
            log_likelihoods.append(np.log(mixture_pdf))
        else:
            log_likelihoods.append(-np.inf)  # Penalize invalid PDFs
    
    # Negative log-likelihood
    return -np.sum(log_likelihoods)

def generate_bivariate_discrete_laplace_with_sigma_param(n, params):
    """
    Generate synthetic data from a k=2 bivariate discrete mixture of Laplace with correlation (via Sigma matrix).
    Args:
        n: Number of samples.
        params: Array of parameters, containing:
            - loc1: Location vector for component 1 (2 values)
            - loc2: Location vector for component 2 (2 values)
            - scale1, scale2: Scales for components 1 and 2 (2 values)
            - Sigma1: Covariance matrix for component 1 (3 independent values)
            - Sigma2: Covariance matrix for component 2 (3 independent values)
            - weight1: Mixture weight for component 1 (1 value)
    Returns:
        Synthetic bivariate data (n x 2 array).
    """
    # Extract parameters
    loc1 = params[:2]  # Location for component 1
    loc2 = params[2:4]  # Location for component 2
    scale1, scale2 = params[4:6]  # Scales for components 1 and 2
    
    Sigma1 = np.array([[params[6], params[7]], [params[7], params[8]]])  # Sigma_1
    Sigma2 = np.array([[params[9], params[10]], [params[10], params[11]]])  # Sigma_2
    
    #weight1 = params[12]  # Mixture weight for component 1
    weight1 = 0.6
    weight2 = 1 - weight1  # Mixture weight for component 2 (ensures weights sum to 1)
    
    # Ensure covariance matrices are valid
    if not (np.all(np.linalg.eigvals(Sigma1) > 0) and np.all(np.linalg.eigvals(Sigma2) > 0)):
        raise ValueError("Covariance matrices must be positive definite.")
    
    # Initialize storage for the data
    data = np.zeros((n, 2))

    for i in range(n):
        # Randomly choose component based on weights
        if i < n*weight1:
            # Component 1
            z = np.random.laplace(0, scale1, size=2)  # Independent Laplace samples
            data[i] = loc1 + np.linalg.cholesky(Sigma1) @ z
        else:
            # Component 2
            z = np.random.laplace(0, scale2, size=2)  # Independent Laplace samples
            data[i] = loc2 + np.linalg.cholesky(Sigma2) @ z
    
    return data

# Example parameters
params = np.array([
    1, 2,          # loc1 (mean of component 1)
    5, 6,          # loc2 (mean of component 2)
    1.0, 1.2,      # scale1, scale2
    1.0, 0.5, 1.0, # Sigma1 (3 values: [var1, cov, var2])
    1.5, -0.3, 0.8 #,# Sigma2 (3 values: [var1, cov, var2])
    #0.6            # weight1 (mixture weight for component 1)
])

# Generate synthetic data
T = 10000
data = generate_bivariate_discrete_laplace_with_sigma_param(T, params)
#data = np.random.rand(T, 2)
results_1 = compute_mle_bivariate_discrete_mixture_laplace(data)
mle_params_lap=results_1.x

print(mle_params_lap)

loc1 = mle_params_lap[:2]  # Location vector for component 1
loc2 = mle_params_lap[2:4]  # Location vector for component 2
scale1, scale2 = mle_params_lap[4:6]  # Scale parameters
Sigma1 = np.array([[mle_params_lap[6], mle_params_lap[7]], [mle_params_lap[7], mle_params_lap[8]]])  # Covariance matrix for component 1
Sigma2 = np.array([[mle_params_lap[9], mle_params_lap[10]], [mle_params_lap[10], mle_params_lap[11]]])  # Covariance matrix for component 2
#weight1 = mle_params[12]  # Mixture weight for component 1
weight1 = 0.6
weight2 = 1 - weight1  # Mixture weight for component 2

# Print results
print("Maximum Likelihood Estimates (MLE):\n")
print(f"Location for Component 1: {loc1}")
print(f"Location for Component 2: {loc2}")
print(f"Scale for Component 1: {scale1:.4f}")
print(f"Scale for Component 2: {scale2:.4f}")
print(f"Covariance Matrix for Component 1:\n{Sigma1}")
print(f"Covariance Matrix for Component 2:\n{Sigma2}")
print(f"Weight for Component 1: {weight1:.4f}")
print(f"Weight for Component 2: {weight2:.4f}")


def compute_aic_bic_laplace(log_likelihood, num_params, num_data):
    aic = -2 * log_likelihood + 2 * num_params
    bic = -2 * log_likelihood + num_params * np.log(num_data)
    return aic, bic

# Compute AIC and BIC
n_params = len(mle_params_lap)
log_likelihood=results_1.fun
print(log_likelihood)
aic, bic = compute_aic_bic_laplace(log_likelihood, n_params, T)

print(f"AIC Lap: {aic:.4f}")
print(f"BIC Lap: {bic:.4f}")

#%% II.3

# Program 1

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

def slog(x):
    # Truncated log to avoid -Inf or +Inf
    realmin = np.finfo(float).tiny
    realmax = np.finfo(float).max
    x_clipped = np.clip(x, realmin, realmax)
    return np.log(x_clipped)

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
    x = np.asarray(x)
    T, d = x.shape  # Extract new shape (1000, 2)
    if d != 2:
        x = x.T
        T, d = x.shape
        #raise ValueError('Not implemented for dimensions other than 2.')

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

T = 10000
np.random.seed(42)
data_new = np.random.rand(10000, 2)
mle_param_t, stderr, iters, loglik_stock, Varcov = MVNCT2estimation(data_new)

results_2 = compute_mle_bivariate_discrete_mixture_laplace(data_new)
mle_params_lap=results_2.x

n_params_lap = len(mle_params_lap)
ll=results_2.fun
print(ll)
aic_lap, bic_lap = compute_aic_bic_laplace(ll, n_params_lap, T)

print(f"AIC Lap: {aic_lap:.4f}")
print(f"BIC Lap: {bic_lap:.4f}")

# Compute AIC and BIC
n_params_t = len(mle_param_t)
print(loglik_stock)
aic_t, bic_t = compute_aic_bic_laplace(-loglik_stock, n_params_t, T)

print(f"AIC t: {aic_t:.4f}")
print(f"BIC t: {bic_t:.4f}")

#%% helpful junk 2

samples = simulate_mixture_bivariate_laplace(pi, mu1, b1, Sigma1, mu2, b2, Sigma2, n)



#%% helpful junk
def bivariate_discrete_laplace_logpdf_previous_version(mean, scale, x):
    """
    Log-probability density of the bivariate discrete Laplace distribution.
    Args:
        mean (tuple): (mu1, mu2), mean of the bivariate Laplace.
        scale (tuple): (b1, b2), scale parameters.
        x (array-like): Observed data points (n x 2 array).
    Returns:
        Log-probability density for each point in x.
    """
    mu1, mu2 = mean
    b1, b2 = scale
    x1, x2 = x[:, 0], x[:, 1]
    return -np.abs(x1 - mu1) / b1 - np.abs(x2 - mu2) / b2 - np.log(2 * b1 * b2)

# Negative log-likelihood for the k=2 mixture model
def negative_log_likelihood_bvlp(params, x):
    """
    Negative log-likelihood for the k=2 component, d=2 discrete mixture of Laplace.
    Args:
        params: [w1, mu1_1, mu1_2, b1_1, b1_2, mu2_1, mu2_2, b2_1, b2_2]
                where w1 is the weight of the first component, and the rest are
                the parameters of the two components (means and scales).
        x: Observed bivariate data points (n x 2 array).
    Returns:
        Negative log-likelihood of the mixture model.
    """
    w1 = params[0]
    mu1 = (params[1], params[2])
    b1 = (params[3], params[4])
    mu2 = (params[5], params[6])
    b2 = (params[7], params[8])
    w2 = 1 - w1  # Enforce weight constraint

    # Compute log-probabilities for each component
    logpdf1 = bivariate_discrete_laplace_logpdf(mu1, b1, x)
    logpdf2 = bivariate_discrete_laplace_logpdf(mu2, b2, x)

    # Mixture log-likelihood
    log_likelihood = np.log(w1 * np.exp(logpdf1) + w2 * np.exp(logpdf2))
    return -np.sum(log_likelihood)  # Negative log-likelihood







# Printing shii


np.random.seed(42)
n_samples = 200
data1 = np.random.randint(-5, 5, size=(n_samples // 2, 2))
data2 = np.random.randint(5, 15, size=(n_samples // 2, 2))
data = np.vstack([data1, data2])

# Example dataset

# Call the MLE computation function
mle_params = compute_mle_bivariate_discrete_mixture_laplace(data)

# Parse and print the results
def print_mle_results(mle_params):
    # Extract parameters
    loc1 = mle_params[:2]  # Location vector for component 1
    loc2 = mle_params[2:4]  # Location vector for component 2
    scale1, scale2 = mle_params[4:6]  # Scale parameters
    Sigma1 = np.array([[mle_params[6], mle_params[7]], [mle_params[7], mle_params[8]]])  # Covariance matrix for component 1
    Sigma2 = np.array([[mle_params[9], mle_params[10]], [mle_params[10], mle_params[11]]])  # Covariance matrix for component 2
    weight1 = mle_params[12]  # Mixture weight for component 1
    weight2 = 1 - weight1  # Mixture weight for component 2

    # Print results
    print("Maximum Likelihood Estimates (MLE):\n")
    print(f"Location for Component 1: {loc1}")
    print(f"Location for Component 2: {loc2}")
    print(f"Scale for Component 1: {scale1:.4f}")
    print(f"Scale for Component 2: {scale2:.4f}")
    print(f"Covariance Matrix for Component 1:\n{Sigma1}")
    print(f"Covariance Matrix for Component 2:\n{Sigma2}")
    print(f"Weight for Component 1: {weight1:.4f}")
    print(f"Weight for Component 2: {weight2:.4f}")

# Print the results
print_mle_results(mle_params)





import numpy as np
import random

def generate_bivariate_discrete_laplace_with_sigma(n, param, loc1, loc2, scale1, scale2, Sigma111, Sigma112, Sigma122, Sigma211, Sigma212, Sigma222, w1):
    """
    Generate synthetic data from a k=2 bivariate discrete mixture of Laplace with correlation (via Sigma matrix).
    Args:
        n: Number of samples.
        w1: Weight of the first component.
        mu1: Mean vector of the first component.
        b1: Scale parameter of the first component.
        sigma1: Covariance matrix of the first component.
        mu2: Mean vector of the second component.
        b2: Scale parameter of the second component.
        sigma2: Covariance matrix of the second component.
    Returns:
        Synthetic bivariate data (n x 2 array).
    """    
    
    # Initialize storage for the data
    data = np.zeros((n, 2))
    
    sigma1 = [[Sigma111, Sigma112],
              [Sigma112, Sigma122]]
    sigma2 = [[Sigma211, Sigma212],
              [Sigma212, Sigma222]]
    for i in range(n):
        # Randomly choose component based on weights
        if np.random.rand() < w1:
            # Component 1
            z = np.random.laplace(0, scale1, size=2)  # Independent Laplace samples
            data[i] = loc1 + np.linalg.cholesky(sigma1) @ z
        else:
            # Component 2
            z = np.random.laplace(0, scale2, size=2)  # Independent Laplace samples
            data[i] = loc2 + np.linalg.cholesky(sigma2) @ z
    return data

# True parameters
w1_true = 0.6
mu1_true = np.array([1, 2])
b1_true = 1.0  # Scale for component 1
sigma1_true = np.array([[1.0, 0.5], [0.5, 1.0]])  # Covariance matrix for component 1

mu2_true = np.array([5, 6])
b2_true = 1.2  # Scale for component 2
sigma2_true = np.array([[1.5, -0.3], [-0.3, 0.8]])  # Covariance matrix for component 2

# Generate synthetic data
n_samples = 1000
data = generate_bivariate_discrete_laplace_with_sigma(
    n=n_samples, w1=w1_true, mu1=mu1_true, b1=b1_true, sigma1=sigma1_true,
    mu2=mu2_true, b2=b2_true, sigma2=sigma2_true
)

# Print first 5 samples for verification
print("Generated data (first 5 samples):")
print(data[:5])