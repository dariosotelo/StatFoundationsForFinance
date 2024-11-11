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

# Question 1 The following is the computation for GAt_pdf and GAt_cdf
def GAt_pdf(z, d, v, theta):
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

    return pdf

# Plot the PDF of the generalized asymmetric t-distribution
z_values = np.linspace(-10, 10, 500)
pdf_values = GAt_pdf(z_values, d=2, v=1.5, theta=1)

plt.figure()
plt.plot(z_values, pdf_values, label='Generalized Asymmetric t-Distribution')

# Plot the PDF of the t-distribution with 3 degrees of freedom
t_pdf_values = t.pdf(z_values*sqrt(2), df=3)*sqrt(2)
plt.plot(z_values, t_pdf_values, label='t-Distribution (df=3)', linestyle='--')

plt.title('PDF Comparison')
plt.xlabel('z')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True)
plt.show()

def GAt_cdf(z, d, v, theta):
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
    
cdf_values = []
# Plot the CDF of the generalized asymmetric t-distribution
for i in z_values:
    cdf_values.append(GAt_cdf(i, d=2, v=1.5, theta=1))

plt.figure()
plt.plot(z_values, cdf_values, label='Generalized Asymmetric t-Distribution CDF')

# Plot the CDF of the t-distribution with 3 degrees of freedom
t_cdf_values = t.cdf(z_values*sqrt(2), df=3)
plt.plot(z_values, t_cdf_values, label='t-Distribution (df=3) CDF', linestyle='--')

plt.title('CDF Comparison')
plt.xlabel('z')
plt.ylabel('Cumulative Probability')
plt.legend()
plt.grid(True)
plt.show()

def ff(x, p, d, v, theta):
    cdf = GAt_cdf(x, d, v, theta)
    return cdf - p

def GAtquantile(p, d, v, theta):
    
    # Find initial bounds for root finding
    lobound = 0
    while ff(lobound, p, d, v, theta) >= 0:
        lobound -= 4
    
    hibound = 0
    while ff(hibound, p, d, v, theta) <= 0:
        hibound += 4

    # Solve for quantile using root_scalar (numerical root finding)
    result = root_scalar(ff, args=(p, d, v, theta), bracket=[lobound, hibound], method='bisect', xtol=1e-5)
    
    if result.converged:
        return result.root
    else:
        raise ValueError("Root-finding did not converge")
    

p = 0.05  # Target probability
d = 2.0   # Shape parameter
v = 1.5   # Degrees of freedom
theta = 1  # Skewness parameter

# Compute quantile for p
quantile = GAtquantile(p, d, v, theta)
print(f"The quantile for p={p} is q={quantile}")

def GAtsim(sim, d, v, theta):
    x = np.zeros(sim)
    lo = 1e-6
    hi = 1 - lo

    for i in tqdm(range(sim)):
        u = np.random.uniform(lo, hi)
        x[i] = GAtquantile(u, d, v, theta)
    
    return x

# Test parameters
sim = 5000
d = 2.0
v = 1.5
theta = 1

# Generate random samples from the GAt distribution
samples = GAtsim(sim, d, v, theta)

# Plot kernel density estimate (KDE) of simulated samples
xx = np.linspace(-5, 5, 1000)
pdf_values = [GAt_pdf(x, d, v, theta) for x in xx]  # Compute theoretical PDF
kde = gaussian_kde(samples)

# Plot
plt.figure(figsize=(10, 6))
plt.plot(xx, kde(xx), 'r-', label='Simulated KDE')
plt.plot(xx, pdf_values, 'b--', label='Theoretical GAt PDF')
plt.title('Generalized Asymmetric t-distribution Simulation')
plt.xlabel('x')
plt.ylabel('Density')
plt.legend()
plt.grid(True)
plt.show()

def conditional_expectation_Zr(r, c, d, v, theta):
    """
    The conditional expectation E[Z^r | Z < c].
    """
    # Compute L based on the formula
    L = v / (v + (-c * theta)**d)
    
    # Beta function values
    numerator = betainc(v - r/d, (r + 1)/d, L) * beta(v - r/d, (r + 1)/d)
    denominator = betainc(v, 1/d, L) * beta(v, 1/d)
    
    # Compute the expectation
    term1 = (-1)**r * v**(r/d)
    term2 = (1 + theta**2) / (theta**r + theta**(r + 2))
    expectation = term1 * term2 * numerator / denominator
    
    return expectation

print(conditional_expectation_Zr(1, quantile, d=2, v=1.5, theta=1))

def mle_student_t(data):
    def neg_log_likelihood(params):
        df, loc, scale = params
        return -np.sum(t.logpdf(data, df, loc, scale))
    initial_params = [10, np.mean(data), np.std(data)]
    bounds = [(2, None), (None, None), (1e-6, None)]
    result = minimize(neg_log_likelihood, initial_params, bounds=bounds)
    return result.x[0], result.x[1], result.x[2]

mle_t_df, mle_t_loc, mle_t_scale = mle_student_t(samples)

#T-ES
def expected_shortfall_t(alpha, df, loc, scale):
    VaR = t.ppf(alpha, df, loc, scale)
    integrand = lambda x: x * t.pdf(x, df, loc, scale)
    ES, _ = quad(integrand, -np.inf, VaR)
    return ES / alpha

alpha = 0.05
ES_t = expected_shortfall_t(alpha, mle_t_df, mle_t_loc, mle_t_scale)
print(f"95% Expected Shortfall (t-distribution): {ES_t}")


def mle_gaussian(data):
    mean_estimate = np.mean(data)
    std_estimate = np.std(data, ddof=1)  # Use ddof=1 for unbiased estimate
    return mean_estimate, std_estimate

mle_gaussian_mean, mle_gaussian_std = mle_gaussian(samples)

def mle_noncentral_t(data):
    def neg_log_likelihood(params):
        df, nc, loc, scale = params
        return -np.sum(stats.nct.logpdf(data, df, nc, loc=loc, scale=scale))
    
    initial_params = [10, 0, np.mean(data), np.std(data)]
    bounds = [(2, None), (None, None), (None, None), (1e-6, None)]
    result = minimize(neg_log_likelihood, initial_params, bounds=bounds)
    
    if result.success:
        return result.x
    else:
        raise ValueError("MLE optimization did not converge")

mle_nct_df, mle_nct_nc, mle_nct_loc, mle_nct_scale = mle_noncentral_t(samples)
print(mle_nct_df, mle_nct_nc, mle_nct_loc, mle_nct_scale)

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
# %%
