#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 18:58:34 2024

"""

#%% Libraries
from scipy.special import betainc
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

#%% Preparation functions

#PDF
def f_GAt(z, d, nu, theta, K):
    if any(param <= 0 for param in [d, nu, theta]):
        return "d, nu, or theta must be positive."
    if z < 0:
        return K * (1 + (-z * theta) ** d / nu) ** -(nu + 1/d)
    else:
        return K * (1 + (z / theta) ** d / nu) ** -(nu + 1/d)

#Auxiliary function for the ES function
def calculate_L(c, nu, theta, d):
    return nu / (nu + (-c * theta) ** d)

#ES
def S_r(c, r, d, nu, theta):
    if c >= 0:
        return "c must be negative."
    L = calculate_L(c, nu, theta, d)
    B_L_num = betainc(nu - r / d, (r + 1) / d, L)  
    B_L_den = betainc(nu, 1 / d, L)                
    return ((-1) ** r) * (nu ** (r / d)) * ((1 + theta ** 2) / (theta ** r + theta ** (r + 2))) * (B_L_num / B_L_den)

#The ES would be calculated using VaR as c and r=1

#CDF
def GAt(z, d, nu, theta, K=1):
    pdf = f_GAt(z, d, nu, theta, K)
    #Underscore because quad returns the error as the second parameter
    cdf, _ = quad(lambda t: f_GAt(t, d, nu, theta, K), -np.inf, z)
    return pdf, cdf


#%% I.1


# Paolella's code suggestion (it is a kinda revised copy-paste from ChatGPT)

def GAtsim(sim, d, v, theta):
    """
    Simulates data from the GAt distribution.
    
    Parameters:
    sim (int): Number of simulations.
    d (float): Parameter d of the GAt distribution.
    v (float): Parameter v of the GAt distribution.
    theta (float): Parameter theta of the GAt distribution.
    
    Returns:
    np.ndarray: Simulated data from the GAt distribution.
    """
    x = np.zeros(sim)
    lo = 1e-6
    hi = 1 - lo
    
    for i in range(sim):
        u = np.random.rand()
        u = max(u, lo)
        u = min(u, hi)
        x[i] = GAtquantile(u, d, v, theta)
        
    return x

def GAtquantile(p, d, v, theta):
    """
    Computes the quantile of the GAt distribution for a given probability p.
    
    Parameters:
    p (float): Probability (0 < p < 1).
    d (float): Parameter d of the GAt distribution.
    v (float): Parameter v of the GAt distribution.
    theta (float): Parameter theta of the GAt distribution.
    
    Returns:
    float: Quantile corresponding to the probability p.
    """
    # Find lower bound for the quantile
    lobound = 0
    while ff(lobound, p, d, v, theta) >= 0:
        lobound -= 4
    
    # Find upper bound for the quantile
    hibound = 0
    while ff(hibound, p, d, v, theta) <= 0:
        hibound += 4
    
    # Use fsolve to find the root of the function, i.e., the quantile
    tol = 1e-5
    q = fsolve(lambda x: ff(x, p, d, v, theta), [lobound, hibound], xtol=tol)[0]
    
    return q





def ff(x, p, d, v, theta, K=1):
    """
    Helper function to compute the difference between the CDF at x and the probability p.
    
    Parameters:
    x (float or array-like): Point at which to evaluate the CDF.
    p (float): Target probability.
    d (float): Shape parameter.
    v (float): Shape parameter.
    theta (float): Skew parameter.
    
    Returns:
    float: Difference between CDF(x) and p.
    """
    # If x is an array, extract the scalar value (fsolve may pass a one-element array)
    if isinstance(x, np.ndarray):
        x = x[0]
        
    # Calculate the PDF and CDF at x
    _, cdf = GAt(x, d, v, theta, K)
    
    # Return the difference between CDF and p
    return cdf - p


#Something is not working but it is related to the array the GAtsim function returns which is then used in ff
#I was thinking on doing a for loop but i havent revised Paolellas code in detail.


# Set parameters
sim = 5000
d = 2.0
v = 1.5
theta = 0.9

# Simulate data
simulated_data = GAtsim(sim, d, v, theta)


#%% Second request to ChatGPT, this one worked.


import numpy as np
from scipy.optimize import fsolve

def GAtsim(sim, d, v, theta):
    """
    Simulates data from the GAt distribution.
    
    Parameters:
    sim (int): Number of simulations.
    d (float): Parameter d of the GAt distribution.
    v (float): Parameter v of the GAt distribution.
    theta (float): Parameter theta of the GAt distribution.
    
    Returns:
    np.ndarray: Simulated data from the GAt distribution.
    """
    x = np.zeros(sim)
    lo = 1e-6
    hi = 1 - lo
    
    for i in range(sim):
        u = np.random.rand()
        u = max(u, lo)
        u = min(u, hi)
        x[i] = GAtquantile(u, d, v, theta)
        
    return x

def GAtquantile(p, d, v, theta):
    """
    Computes the quantile of the GAt distribution for a given probability p.
    
    Parameters:
    p (float): Probability (0 < p < 1).
    d (float): Parameter d of the GAt distribution.
    v (float): Parameter v of the GAt distribution.
    theta (float): Parameter theta of the GAt distribution.
    
    Returns:
    float: Quantile corresponding to the probability p.
    """
    # Find lower bound for the quantile
    lobound = 0
    while ff(lobound, p, d, v, theta) >= 0:
        lobound -= 4
    
    # Find upper bound for the quantile
    hibound = 0
    while ff(hibound, p, d, v, theta) <= 0:
        hibound += 4
    
    # Use fsolve to find the root of the function, i.e., the quantile
    tol = 1e-5
    q = fsolve(lambda x: ff(x, p, d, v, theta), x0=(lobound + hibound) / 2, xtol=tol)[0]
    
    return q

def ff(x, p, d, v, theta):
    """
    Helper function to compute the difference between the CDF at x and the probability p.
    
    Parameters:
    x (float): Point at which to evaluate the CDF.
    p (float): Target probability.
    d (float): Parameter d of the GAt distribution.
    v (float): Parameter v of the GAt distribution.
    theta (float): Parameter theta of the GAt distribution.
    
    Returns:
    float: Difference between CDF(x) and p.
    """
    _, cdf = GAt(x, d, v, theta)
    return cdf - p

def GAt(z, d, v, theta, K=1):
    """
    Computes the PDF and approximates the CDF of the GAt distribution.
    
    Parameters:
    z (float): Point at which to evaluate the PDF and CDF.
    d (float): Shape parameter.
    v (float): Shape parameter.
    theta (float): Skew parameter.
    K (float): Normalizing constant, if known.
    
    Returns:
    tuple: (pdf, cdf) values at point z.
    """
    # Calculate PDF using the GAt distribution formula
    pdf = f_GAt(z, d, v, theta, K)
    
    # Approximate CDF by integrating the PDF from -∞ to z
    from scipy.integrate import quad
    cdf, _ = quad(lambda t: f_GAt(t, d, v, theta, K), -np.inf, z)
    
    return pdf, cdf

def f_GAt(z, d, nu, theta, K=1):
    """
    PDF of the GAt distribution.
    
    Parameters:
    z (float): Point at which to evaluate the PDF.
    d (float): Shape parameter.
    nu (float): Shape parameter.
    theta (float): Skew parameter.
    K (float): Normalizing constant, if known.
    
    Returns:
    float: PDF value at z.
    """
    if any(param <= 0 for param in [d, nu, theta]):
        raise ValueError("d, nu, or theta must be positive.")
    
    if z < 0:
        return K * (1 + (-z * theta) ** d / nu) ** -(nu + 1 / d)
    else:
        return K * (1 + (z / theta) ** d / nu) ** -(nu + 1 / d)


sim = 5000
d = 2.0
v = 1.5
theta = 0.9

simulated_data = GAtsim(sim, d, v, theta)


# %%
from scipy.stats import norm
from math import sqrt, exp
import scipy
import math
import numpy as np
from pylab import plot, show, grid, xlabel, ylabel
from statistics import mean
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import random
from scipy.integrate import quad
from scipy.stats import t
from scipy.optimize import root_scalar
from scipy.special import beta, betainc
from scipy.stats import gaussian_kde
import scipy.stats as stats
from scipy.optimize import minimize
from scipy.special import kv, gamma  
from scipy.fft import fftshift, ifft
import seaborn as sns
import pandas as pd
np.random.seed()
# Question 1 The following is the computation for GAt_pdf and GAt_cdf
def GAt_pdf(z, d, v, theta, loc=0, scale=1):
    # Apply location and scale transformation
    z = (z - loc) / scale
    
    # Calculate the normalization constant C_d,v,theta
    C_inv = (theta**-1 + theta) * d**-1 * v**(1/d) * beta(1/d, v)
    C = 1 / C_inv
    
    # Determine the PDF based on the sign of z
    if isinstance(z, (float, int)):  # Single value case
        if z < 0:
            pdf = C * (1 + ((-z * theta)**d) / v)**(-(v + 1/d))
        else:
            pdf = C * (1 + ((z / theta)**d) / v)**(-(v + 1/d))
    else:  # Handle array-like input
        pdf = np.where(
            z < 0,
            C * (1 + ((-z * theta)**d) / v)**(-(v + 1/d)),
            C * (1 + ((z / theta)**d) / v)**(-(v + 1/d))
        )

    # Adjust for scale
    pdf /= scale

    return pdf

# # Plot the PDF of the generalized asymmetric t-distribution
# z_values = np.linspace(-10, 10, 500)
# pdf_values = GAt_pdf(z_values, d=2, v=1.5, theta=1)

# plt.figure()
# plt.plot(z_values, pdf_values, label='Generalized Asymmetric t-Distribution')

# # Plot the PDF of the t-distribution with 3 degrees of freedom
# t_pdf_values = t.pdf(z_values*sqrt(2), df=3)*sqrt(2)
# plt.plot(z_values, t_pdf_values, label='t-Distribution (df=3)', linestyle='--')

# plt.title('PDF Comparison')
# plt.xlabel('z')
# plt.ylabel('Probability Density')
# plt.legend()
# plt.grid(True)
# plt.show()

def GAt_cdf(z, d, v, theta, loc=0, scale=1):
    z = (z - loc) / scale  # Apply location and scale transformation

    def B_L(a, b):
        return betainc(a, b, L)

    def B_U(a, b):
        return betainc(a, b, U)

    # Compute L and U as piecewise values depend on z
    if isinstance(z, (float, int)):
        if z <= 0:
            L = v / (v + (-z * theta)**d)
            return B_L(v, 1/d) / (1 + theta**2)
        else:
            U = (z / theta)**d / (v + (z / theta)**d)
            return B_U(1/d, v) / (1 + theta**-2) + (1 + theta**2)**-1
    else:
        # Array operations
        cdf = np.zeros_like(z)
        L = v / (v + (-z[z < 0] * theta)**d)
        cdf[z <= 0] = B_L(v, 1/d) / (1 + theta**2)
        U = (z[z > 0] / theta)**d / (v + (z / theta)**d)
        cdf[z > 0] = B_U(1/d, v) / (1 + theta**2)
    return cdf
    
# cdf_values = []
# # Plot the CDF of the generalized asymmetric t-distribution
# for i in z_values:
#     cdf_values.append(GAt_cdf(i, d=1.71, v=1.96, theta=0.879))

# plt.figure()
# plt.plot(z_values, cdf_values, label='Generalized Asymmetric t-Distribution CDF')

# # Plot the CDF of the t-distribution with 3 degrees of freedom
# t_cdf_values = t.cdf(z_values*sqrt(2), df=3)
# plt.plot(z_values, t_cdf_values, label='t-Distribution (df=3) CDF', linestyle='--')

# plt.title('CDF Comparison')
# plt.xlabel('z')
# plt.ylabel('Cumulative Probability')
# plt.legend()
# plt.grid(True)
# plt.show()

def ff(x, p, d, v, theta, loc, scale):
    cdf = GAt_cdf(x, d, v, theta, loc, scale)
    return cdf - p

def GAtquantile(p, d, v, theta, loc=0, scale=1):
    
    # Find initial bounds for root finding
    lobound = loc
    while ff(lobound, p, d, v, theta, loc, scale) >= 0:
        lobound -= 4 * scale
    
    hibound = loc
    while ff(hibound, p, d, v, theta, loc, scale) <= 0:
        hibound += 4 * scale

    # Solve for quantile using root_scalar (numerical root finding)
    result = root_scalar(ff, args=(p, d, v, theta, loc, scale), bracket=[lobound, hibound], method='bisect', xtol=1e-5)
    
    if result.converged:
        return result.root
    else:
        raise ValueError("Root-finding did not converge")
    

p = 0.01  # Target probability
d = 1.71  # Shape parameter
v = 1.96   # Degrees of freedom
theta = 0.879  # Skewness parameter
loc = 0.175  # Location parameter
scale = 1.21  # Scale parameter

# Compute quantile for p
quantile = GAtquantile(p, d, v, theta, loc, scale)
print(f"The quantile for p={p} is q={quantile}")

def GAtsim(sim, d, v, theta, loc=0, scale=1):
    x = np.zeros(sim)
    lo = 1e-6
    hi = 1 - lo

    for i in tqdm(range(sim)):
        u = np.random.uniform(lo, hi)
        x[i] = GAtquantile(u, d, v, theta, loc, scale)
    
    return x

# Test parameters
sim = 500
d = 1.71  # Shape parameter
v = 1.96   # Degrees of freedom
theta = 0.879  # Skewness parameter
loc = 0.175  # Location parameter
scale = 1.21  # Scale parameter

# Generate random samples from the GAt distribution
samples = GAtsim(sim, d, v, theta, loc, scale)

def conditional_expectation_Zr(r=1, c=quantile, d=1.71, v=1.96, theta=0.879, loc=0.175, scale=1.21):
    """
    The conditional expectation E[Z^r | Z < c] with location and scale parameters.
    """
    # Apply location and scale transformation
    c = (c - loc) / scale
    
    # Compute L based on the formula
    L = v / (v + (-c * theta)**d)
    
    # Beta function values
    numerator = betainc(v - r/d, (r + 1)/d, L) * beta(v - r/d, (r + 1)/d)
    denominator = betainc(v, 1/d, L) * beta(v, 1/d)
    
    # Compute the expectation
    term1 = (-1)**r * v**(r/d)
    term2 = (1 + theta**2) / (theta**r + theta**(r + 2))
    expectation = term1 * term2 * numerator / denominator
    
    # Adjust for scale
    expectation *= scale**r
    
    return expectation

def nonparametric_bootstrap(data, ESlevel=0.01, B=500, n=500):
    ES_samples = []
    VaR_samples = []
    for _ in range(B):
        bootstrap_sample = np.random.choice(data, n, replace=True)
        ES = expected_shortfall(bootstrap_sample, ESlevel)
        VaR = np.percentile(bootstrap_sample, ESlevel * 100)
        ES_samples.append(ES)
        VaR_samples.append(VaR)

    return ES_samples, VaR_samples

def mle_student_t(data):
    def neg_log_likelihood(params):
        df, loc, scale = params
        return -np.sum(t.logpdf(data, df, loc, scale))
    initial_params = [10, np.mean(data), np.std(data)]
    bounds = [(2, None), (None, None), (1e-6, None)]
    result = minimize(neg_log_likelihood, initial_params, bounds=bounds)
    return result.x[0], result.x[1], result.x[2]

mle_t_df, mle_t_loc, mle_t_scale = mle_student_t(samples)

def expected_shortfall(data, ESlevel=0.01):
    sorted_data = np.sort(data)
    index = int(ESlevel * len(data))
    ES = np.mean(sorted_data[:index])
    return ES

def parametric_bootstrap_t(data, ESlevel=0.01, B=500, n=500):
    # Generate new set of data samples using the MLE parameters
    df_mle, loc_mle, scale_mle = mle_student_t(data)
    ES_samples = []
    VaR_samples = []
    for _ in range(B):
        bootstrap_sample = t.rvs(df_mle, loc_mle, scale_mle, size=n)
        ES = expected_shortfall(bootstrap_sample, ESlevel)
        VaR = np.percentile(bootstrap_sample, ESlevel * 100)
        ES_samples.append(ES)
        VaR_samples.append(VaR)
    return ES_samples, VaR_samples

def mle_gaussian(data):
    mean_estimate = np.mean(data)
    std_estimate = np.std(data, ddof=1)  # Use ddof=1 for unbiased estimate
    return mean_estimate, std_estimate

def parametric_bootstrap_gaussian(data, ESlevel=0.01, B=500, n=500):
    # Generate new set of data samples using the MLE parameters
    mean_mle, std_mle = mle_gaussian(data)
    ES_samples = []
    VaR_samples = []
    for _ in range(B):
        bootstrap_sample = np.random.normal(mean_mle, std_mle, size=n)
        ES = expected_shortfall(bootstrap_sample, ESlevel)
        VaR = np.percentile(bootstrap_sample, ESlevel * 100)
        ES_samples.append(ES)
        VaR_samples.append(VaR)
    return ES_samples, VaR_samples

def mle_noncentral_t(data):
    # Negative log-likelihood function
    def negative_log_likelihood(params):
        nu, delta, loc, scale = params
        # Enforce constraints: nu > 0, scale > 0
        if nu <= 0 or scale <= 0:
            return np.inf
        # Compute the negative log-likelihood
        nll = -np.sum(stats.nct.logpdf(data, nu, delta, loc=loc, scale=scale))
        return nll

    # Initial guesses for parameters
    nu_init = 5.0
    delta_init = 0.0
    loc_init = 0
    scale_init = 1
    params_init = [nu_init, delta_init, loc_init, scale_init]

    # Bounds: nu > 0, scale > 0; delta and loc are unbounded
    bounds = [
        (1e-6, None),  # nu > 0
        (None, None),  # delta unbounded
        (None, None),  # loc unbounded
        (1e-6, None)   # scale > 0
    ]

    # Optimize the negative log-likelihood
    result = minimize(negative_log_likelihood, params_init, bounds=bounds, method='L-BFGS-B')

    nu_mle, delta_mle, loc_mle, scale_mle = result.x
    return nu_mle, delta_mle, loc_mle, scale_mle

def parametric_bootstrap_nonc(data, ESlevel=0.01, B=500, n=500):
    # Generate new set of data samples using the MLE parameters
    nu_mle, delta_mle, loc_mle, scale_mle = mle_noncentral_t(data)
    ES_samples = []
    VaR_samples = []
    for _ in range(B):
        bootstrap_sample = stats.nct.rvs(nu_mle, delta_mle, loc=loc_mle, scale=scale_mle, size=n)
        ES = expected_shortfall(bootstrap_sample, ESlevel)
        VaR = np.percentile(bootstrap_sample, ESlevel * 100)
        ES_samples.append(ES)
        VaR_samples.append(VaR)
    return ES_samples, VaR_samples

# mle_nct_df, mle_nct_nc, mle_nct_loc, mle_nct_scale = mle_noncentral_t(samples)
# test = stats.nct.rvs(mle_nct_df, mle_nct_nc, loc=mle_nct_loc, scale=mle_nct_scale, size=10000)

# # Filter values between -5 and 5
# test = test[(test >= -5) & (test <= 5)]

# # Plot the histogram of the noncentral t-distribution samples
# plt.figure(figsize=(10, 6))
# plt.hist(test, bins=50, density=True, alpha=0.6, color='g', label='Noncentral t-distribution Samples')

# # Plot the theoretical PDF of the noncentral t-distribution
# x = np.linspace(-5, 5, 1000)
# pdf = stats.nct.pdf(x, mle_nct_df, mle_nct_nc, loc=mle_nct_loc, scale=mle_nct_scale)
# plt.plot(x, pdf, 'r-', label='Theoretical Noncentral t-distribution PDF')

# # Plot the theoretical PDF of the GAt distribution
# GAt_pdf_values = GAt_pdf(x, d=1.71, v=1.96, theta=0.879)
# plt.plot(x, GAt_pdf_values, 'b--', label='Theoretical GAt PDF')

# plt.title('Histogram and PDF of Noncentral t-distribution Samples')
# plt.xlabel('Value')
# plt.ylabel('Density')
# plt.legend()
# plt.grid(True)
# plt.show()
#HOW TF DO I DO MLE FOR GAt???????
# def mle_GAt(data):
#     def neg_log_likelihood(params):
#         d, v, theta = params
#         return -np.sum(np.log(GAt_pdf(data, d, v, theta)))
    
#     initial_params = [2, 1, 0]
#     bounds = [(1e-6, None), (1e-6, None), (None, None)]
#     result = minimize(neg_log_likelihood, initial_params, bounds=bounds)
    
#     if result.success:
#         return result.x
#     else:
#         raise ValueError("MLE optimization did not converge")
    
# print(mle_GAt(samples))

def estimate_mixed_normals(data, max_iter=1000, tol=1e-6):
    data = np.asarray(data)
    n = data.size

    # Initialize parameters
    # Randomly assign data to clusters to initialize responsibilities
    responsibilities = np.random.rand(n, 2)
    responsibilities /= responsibilities.sum(axis=1, keepdims=True)

    # Initial parameter estimates
    pi1 = responsibilities[:, 0].mean()
    pi2 = 1 - pi1

    mu1 = (responsibilities[:, 0] @ data) / responsibilities[:, 0].sum()
    mu2 = (responsibilities[:, 1] @ data) / responsibilities[:, 1].sum()

    sigma1 = np.sqrt(((responsibilities[:, 0] * (data - mu1)**2).sum()) / responsibilities[:, 0].sum())
    sigma2 = np.sqrt(((responsibilities[:, 1] * (data - mu2)**2).sum()) / responsibilities[:, 1].sum())

    log_likelihood_old = None

    for iteration in range(max_iter):
        # E-step: compute responsibilities
        # Compute the probability density for each component
        pdf1 = pi1 * (1 / (np.sqrt(2 * np.pi) * sigma1)) * np.exp(-0.5 * ((data - mu1) / sigma1)**2)
        pdf2 = pi2 * (1 / (np.sqrt(2 * np.pi) * sigma2)) * np.exp(-0.5 * ((data - mu2) / sigma2)**2)

        total_pdf = pdf1 + pdf2
        responsibilities = np.vstack((pdf1, pdf2)).T / total_pdf[:, np.newaxis]

        # M-step: update parameters
        Nk = responsibilities.sum(axis=0)
        pi1 = Nk[0] / n
        pi2 = Nk[1] / n

        mu1 = (responsibilities[:, 0] @ data) / Nk[0]
        mu2 = (responsibilities[:, 1] @ data) / Nk[1]

        sigma1 = np.sqrt((responsibilities[:, 0] * (data - mu1)**2).sum() / Nk[0])
        sigma2 = np.sqrt((responsibilities[:, 1] * (data - mu2)**2).sum() / Nk[1])

        # Check for convergence
        log_likelihood = np.sum(np.log(total_pdf))
        if log_likelihood_old is not None and abs(log_likelihood - log_likelihood_old) < tol:
            break
        log_likelihood_old = log_likelihood

    return pi1, mu1, sigma1, mu2, sigma2

def generate_mixed_normal_samples(pi1, mu1, sigma1, mu2, sigma2, size=1000):

    # Generate component labels based on mixing proportion
    component_labels = np.random.choice([0, 1], size=size, p=[pi1, 1 - pi1])
    
    # Generate samples for each component
    samples = np.where(
        component_labels == 0,
        np.random.normal(mu1, sigma1, size),
        np.random.normal(mu2, sigma2, size)
    )
    
    return samples

def parametric_bootstrap_mixed(data, ESlevel=0.01, B=500, n=500):
    # Generate new set of data samples using the MLE parameters
    pi1, mu1, sigma1, mu2, sigma2 = estimate_mixed_normals(data)
    ES_samples = []
    VaR_samples = []
    for _ in range(B):
        bootstrap_sample = generate_mixed_normal_samples(pi1, mu1, sigma1, mu2, sigma2, size=n)
        ES = expected_shortfall(bootstrap_sample, ESlevel)
        VaR = np.percentile(bootstrap_sample, ESlevel * 100)
        ES_samples.append(ES)
        VaR_samples.append(VaR)
    return ES_samples, VaR_samples

# Generate parametric bootstrap samples
nonparametric_bootstrap_samples, var_nonpara = nonparametric_bootstrap(samples)
bootstrap_samples_t, var_t = parametric_bootstrap_t(samples)
bootstrap_samples_nonc, var_nonc = parametric_bootstrap_nonc(samples)
bootstrap_samples_gaussian, var_gaussian = parametric_bootstrap_gaussian(samples)
bootstrap_samples_mixed, var_mixed = parametric_bootstrap_mixed(samples)

# Calculate expected shortfall for the original data
original_ES = conditional_expectation_Zr(1, quantile, d=1.71, v=1.96, theta=0.879, loc=0.175, scale=1.21)

# Create a DataFrame for plotting
data = {
    'Expected Shortfall': (nonparametric_bootstrap_samples + bootstrap_samples_t + 
                           bootstrap_samples_nonc + bootstrap_samples_gaussian + 
                           bootstrap_samples_mixed),
    'Type': (['Nonparametric Bootstrap'] * len(nonparametric_bootstrap_samples) +
             ['Bootstrap t'] * len(bootstrap_samples_t) +
             ['Bootstrap Noncentral t'] * len(bootstrap_samples_nonc) +
             ['Bootstrap Gaussian'] * len(bootstrap_samples_gaussian) +
             ['Bootstrap Mixed'] * len(bootstrap_samples_mixed))
}
df = pd.DataFrame(data)

# Filter the DataFrame to show expected shortfall from 0 to -15
df = df[(df['Expected Shortfall'] >= -15) & (df['Expected Shortfall'] <= 0)]

# Plot the boxplot with outliers represented by a + sign
plt.figure(figsize=(10, 6))
sns.boxplot(x='Type', y='Expected Shortfall', data=df, flierprops=dict(marker='+', color='red', markersize=8))
plt.axhline(y=original_ES, color='r', linestyle='--', label='Original ES')
plt.title('Boxplot of Expected Shortfall')
plt.xlabel('Type')
plt.ylabel('Expected Shortfall')
plt.legend()
plt.grid(True)
plt.show()

# Create a DataFrame for plotting VaR
data_var = {
    'VaR': (var_nonpara + var_t + var_nonc + var_gaussian + var_mixed),
    'Type': (['Nonparametric Bootstrap'] * len(var_nonpara) +
             ['Bootstrap t'] * len(var_t) +
             ['Bootstrap Noncentral t'] * len(var_nonc) +
             ['Bootstrap Gaussian'] * len(var_gaussian) +
             ['Bootstrap Mixed'] * len(var_mixed))
}
df_var = pd.DataFrame(data_var)

# Filter the DataFrame to show VaR from 0 to -10
df_var = df_var[(df_var['VaR'] >= -10) & (df_var['VaR'] <= 0)]

# Plot the boxplot for VaR with outliers represented by a + sign
plt.figure(figsize=(10, 6))
sns.boxplot(x='Type', y='VaR', data=df_var, flierprops=dict(marker='+', color='red', markersize=8))
plt.axhline(y=quantile, color='r', linestyle='--', label='True VaR')
plt.title('Boxplot of VaR')
plt.xlabel('Type')
plt.ylabel('VaR')
plt.legend()
plt.grid(True)
plt.show()

samples_t = np.random.standard_t(4, 500)

def student_t_ES(df, loc=0, scale=1, ESlevel=0.01):
    return -scipy.stats.t.pdf(scipy.stats.t.ppf(ESlevel, df), df, loc, scale)/scipy.stats.t.cdf(scipy.stats.t.ppf(ESlevel, df), df, loc, scale)*(df+scipy.stats.t.ppf(ESlevel, df)**2)/(df-1)

def VaRt(df, alpha=0.01):
    return scipy.stats.t.ppf(alpha, df)

nonparametric_bootstrap_samples,var_nonpara = nonparametric_bootstrap(samples_t)
bootstrap_samples_t, var_t = parametric_bootstrap_t(samples_t)
bootstrap_samples_nonc, var_nonc = parametric_bootstrap_nonc(samples_t)
bootstrap_samples_gaussian, var_gaussian = parametric_bootstrap_gaussian(samples_t)
bootstrap_samples_mixed, var_mixed = parametric_bootstrap_mixed(samples_t)
# Calculate expected shortfall for the original data
original_ES = student_t_ES(4)


# Create a DataFrame for plotting
data = {
    'Expected Shortfall': (nonparametric_bootstrap_samples + bootstrap_samples_t + 
                           bootstrap_samples_nonc + bootstrap_samples_gaussian + 
                           bootstrap_samples_mixed),
    'Type': (['Nonparametric Bootstrap'] * len(nonparametric_bootstrap_samples) +
             ['Bootstrap t'] * len(bootstrap_samples_t) +
             ['Bootstrap Noncentral t'] * len(bootstrap_samples_nonc) +
             ['Bootstrap Gaussian'] * len(bootstrap_samples_gaussian) +
             ['Bootstrap Mixed'] * len(bootstrap_samples_mixed))
}
df = pd.DataFrame(data)

# Filter the DataFrame to show expected shortfall from 0 to -15
df = df[(df['Expected Shortfall'] >= -15) & (df['Expected Shortfall'] <= 0)]

# Plot the boxplot with outliers represented by a + sign
plt.figure(figsize=(10, 6))
sns.boxplot(x='Type', y='Expected Shortfall', data=df, flierprops=dict(marker='+', color='red', markersize=8))
plt.axhline(y=original_ES, color='r', linestyle='--', label='Original ES')
plt.title('Boxplot of Expected Shortfall')
plt.xlabel('Type')
plt.ylabel('Expected Shortfall')
plt.legend()
plt.grid(True)
plt.show()

data_var = {
    'VaR': (var_nonpara + var_t + var_nonc + var_gaussian + var_mixed),
    'Type': (['Nonparametric Bootstrap'] * len(var_nonpara) +
             ['Bootstrap t'] * len(var_t) +
             ['Bootstrap Noncentral t'] * len(var_nonc) +
             ['Bootstrap Gaussian'] * len(var_gaussian) +
             ['Bootstrap Mixed'] * len(var_mixed))
}
df_var = pd.DataFrame(data_var)

# Filter the DataFrame to show VaR from 0 to -10
df_var = df_var[(df_var['VaR'] >= -10) & (df_var['VaR'] <= 0)]

# Plot the boxplot for VaR with outliers represented by a + sign
plt.figure(figsize=(10, 6))
sns.boxplot(x='Type', y='VaR', data=df_var, flierprops=dict(marker='+', color='red', markersize=8))
plt.axhline(y=VaRt(4), color='r', linestyle='--', label='True VaR')
plt.title('Boxplot of VaR')
plt.xlabel('Type')
plt.ylabel('VaR')
plt.legend()
plt.grid(True)
plt.show()
