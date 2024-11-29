#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 14:55:03 2024

@author: 
"""

# Libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm  # For colormap
from scipy.integrate import quad
import scipy.special as sp
from scipy.optimize import minimize
from scipy.special import gammaln, gamma, kv
from scipy.linalg import cholesky



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
mu1 = [0, 0]
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
x_values = np.linspace(0, 20, 200)
    
graph_difference(nu, x_values)


# Part b
# MLE of the k=2 component, d=2 (bivariate) discrete mixture of Laplace.
def compute_mle_bivariate_discrete_mixture_laplace(data):
    #initial_params = np.random.rand(13)  # Locations, scales, Sigma terms, weight
    initial_params = np.array([
        1.2, 2.1,          # loc1 (mean of component 1)
        5.2, 6.9,          # loc2 (mean of component 2)
        1.0, 0.9,      # scale1, scale2
        1.0, 0.02, 1.0, # Sigma1 (3 values: [var1, cov, var2])
        1.4, -0.2, 0.6,# Sigma2 (3 values: [var1, cov, var2])
        0.5            # weight1 (mixture weight for component 1)
    ])

    # Minimize the negative log-likelihood using BFGS
    result = minimize(
        fun = negative_log_likelihood_two_component_mixture_bivariate_laplace,
        x0=initial_params,
        args=(data,),
        method='BFGS'
    )
    return result.x


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
    
    weight1 = params[12]  # Mixture weight for component 1
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
    
    weight1 = params[12]  # Mixture weight for component 1
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
    1.5, -0.3, 0.8,# Sigma2 (3 values: [var1, cov, var2])
    0.6            # weight1 (mixture weight for component 1)
])

# Generate synthetic data
T = 1000
data = generate_bivariate_discrete_laplace_with_sigma_param(T, params)

params1 = compute_mle_bivariate_discrete_mixture_laplace(data)

def print_mle_results(params):
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
print_mle_results(params1)







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

















