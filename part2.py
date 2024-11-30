from scipy.stats import norm, t, gaussian_kde, levy_stable
from scipy.special import beta, betainc, kv, gamma
from scipy.optimize import root_scalar, minimize, brentq
from scipy.fft import fftshift, ifft
from scipy.integrate import quad
from math import sqrt, exp
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
import seaborn as sns
import pandas as pd
from tqdm import tqdm
import random
import scipy.stats as stats
from scipy.linalg import cho_solve, cho_factor
from scipy.special import gammaln
from scipy.stats import gamma, multivariate_normal
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import minimize
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import inv, det

def bivariate_laplace_pdf(y, mu, Sigma, b):
    """
    Compute the PDF of the bivariate Laplace distribution.

    Parameters:
        y: ndarray
            Data point (1D array of shape (2,)).
        mu: ndarray
            Mean vector (1D array of shape (2,)).
        Sigma: ndarray
            Covariance matrix (2x2 array).
        b: float
            Scale parameter.

    Returns:
        pdf: float
            Probability density at y.
    """
    d = len(mu)
    diff = y - mu
    delta = diff.T @ inv(Sigma) @ diff  # Mahalanobis distance
    norm_const = (2 * np.pi)**(-d / 2) * det(Sigma)**(-1 / 2) * b**(-d)
    pdf = norm_const * np.exp(-np.sqrt(delta) / b)
    return pdf

def negative_log_likelihood_fixed_lambdas(params, x, lambdas):
    """
    Compute the negative log-likelihood for the 2-component bivariate Laplace mixture model
    with fixed mixture weights (lambdas).

    Parameters:
        params: ndarray
            Flattened parameter vector (remaining 12 parameters).
        x: ndarray
            Data points of shape (n_samples, 2).
        lambdas: list
            Fixed mixture weights [lambda1, lambda2].

    Returns:
        nll: float
            Negative log-likelihood.
    """
    n_samples, d = x.shape

    # Unpack parameters
    mu1 = params[0:2]
    mu2 = params[2:4]
    Sigma1 = np.array([[params[4], params[5]], [params[5], params[6]]])
    Sigma2 = np.array([[params[7], params[8]], [params[8], params[9]]])
    b1, b2 = params[10], params[11]

    # Ensure covariance matrices are positive definite
    if np.min(np.linalg.eigvals(Sigma1)) <= 1e-10 or np.min(np.linalg.eigvals(Sigma2)) <= 1e-10:
        return 1e10  # Large penalty for invalid covariance matrices

    # Compute the mixture PDF
    lambda1, lambda2 = lambdas
    pdf_values = np.zeros(n_samples)
    for i in range(n_samples):
        pdf1 = bivariate_laplace_pdf(x[i], mu1, Sigma1, b1)
        pdf2 = bivariate_laplace_pdf(x[i], mu2, Sigma2, b2)
        pdf_values[i] = lambda1 * pdf1 + lambda2 * pdf2

    # Compute the negative log-likelihood
    log_pdf = np.log(pdf_values + 1e-10)  # Add small value to avoid log(0)
    nll = -np.sum(log_pdf)
    return nll

def fit_bivariate_laplace_mixture_fixed_lambdas(x, lambdas, max_iters=300):
    """
    Fit a 2-component bivariate Laplace mixture model with fixed mixture weights.

    Parameters:
        x: ndarray
            Data points of shape (n_samples, 2).
        lambdas: list
            Fixed mixture weights [lambda1, lambda2].
        max_iters: int
            Maximum number of iterations for the optimizer.

    Returns:
        result: dict
            Estimated parameters.
    """
    n_samples, d = x.shape

    # Initial guess for remaining 12 parameters
    init_params = np.array([
        0, 0,  # Initial mu1
        1, 1,  # Initial mu2
        1, 0, 1,  # Initial Sigma1 (diag + off-diagonal)
        1, 0, 1,  # Initial Sigma2 (diag + off-diagonal)
        1, 1  # Initial scale parameters b1, b2
    ])

    # Bounds for parameters to ensure validity
    bounds = [
        (-10, 10), (-10, 10),  # Bounds for mu1
        (-10, 10), (-10, 10),  # Bounds for mu2
        (0.01, 10), (-5, 5), (0.01, 10),  # Bounds for Sigma1 (diag + off-diagonal)
        (0.01, 10), (-5, 5), (0.01, 10),  # Bounds for Sigma2 (diag + off-diagonal)
        (0.01, 10), (0.01, 10)  # Bounds for b1, b2
    ]

    # Optimize negative log-likelihood
    result = minimize(
        negative_log_likelihood_fixed_lambdas,
        init_params,
        args=(x, lambdas),
        method='L-BFGS-B',
        bounds=bounds,
        options={'disp': True, 'maxiter': max_iters}
    )

    # Extract optimized parameters
    optimized_params = result.x
    mu1 = optimized_params[0:2]
    mu2 = optimized_params[2:4]
    Sigma1 = np.array([[optimized_params[4], optimized_params[5]], [optimized_params[5], optimized_params[6]]])
    Sigma2 = np.array([[optimized_params[7], optimized_params[8]], [optimized_params[8], optimized_params[9]]])
    b1, b2 = optimized_params[10], optimized_params[11]

    return {
        'mu1': mu1,
        'mu2': mu2,
        'Sigma1': Sigma1,
        'Sigma2': Sigma2,
        'b1': b1,
        'b2': b2,
        'nll': result.fun,
        'success': result.success
    }





import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm  # For colormap
from scipy.integrate import quad
import scipy.special as sp
from scipy.optimize import minimize
from scipy.special import gammaln, gamma, kv
from scipy.linalg import cholesky

def bivariate_laplace_samples(n, mu, Sigma, b):
    """
    Generate samples from a bivariate Laplace distribution.
    
    Parameters:
        n: int
            Number of samples to generate.
        mu: ndarray
            Mean vector of shape (2,).
        Sigma: ndarray
            Covariance matrix of shape (2, 2).
        b: float
            Scale parameter.

    Returns:
        samples: ndarray
            Generated samples of shape (n, 2).
    """
    d = len(mu)
    z = np.random.laplace(size=(n, d))  # Laplace-distributed random variables
    eigvals, eigvecs = np.linalg.eigh(Sigma)
    A = eigvecs @ np.diag(np.sqrt(eigvals))  # Covariance transformation
    samples = mu + np.sqrt(b) * z @ A.T
    return samples

def mixture_bivariate_laplace(n, lambdas, mus, Sigmas, bs):
    """
    Generate samples from a mixture of bivariate Laplace distributions.
    
    Parameters:
        n: int
            Number of samples to generate.
        lambdas: list
            Mixture weights, e.g., [0.5, 0.5].
        mus: list of ndarrays
            List of mean vectors for each component.
        Sigmas: list of ndarrays
            List of covariance matrices for each component.
        bs: list of floats
            List of scale parameters for each component.

    Returns:
        samples: ndarray
            Generated samples of shape (n, 2).
    """
    k = len(lambdas)  # Number of components
    assert len(mus) == k and len(Sigmas) == k and len(bs) == k, "Mismatch in components"
    assert np.isclose(np.sum(lambdas), 1), "Mixture weights must sum to 1"

    # Step 1: Draw component labels based on the mixture weights
    component_labels = np.random.choice(k, size=n, p=lambdas)

    # Step 2: Generate samples for each component
    samples = np.zeros((n, 2))  # Storage for samples
    for c in range(k):
        n_c = np.sum(component_labels == c)  # Number of samples for component c
        if n_c > 0:
            samples[component_labels == c, :] = bivariate_laplace_samples(n_c, mus[c], Sigmas[c], bs[c])
    
    return samples

# Example usage
n_samples = 1000
lambdas = [0.7, 0.3]  # Mixture weights
mus = [np.array([0, 0]), np.array([0, 0])]  # Mean vectors
Sigmas = [np.array([[1, 0.5], [0.5, 1]]), np.array([[4, 2], [2, 4]])]  # Covariance matrices
bs = [10, 5]  # Scale parameters

# Generate samples
# samples = mixture_bivariate_laplace(n_samples, lambdas, mus, Sigmas, bs)
# print(samples)
# Fit the model
np.random.seed(42)
samples = np.random.randn(1000,2)
result = fit_bivariate_laplace_mixture_fixed_lambdas(samples, lambdas)

# Display results
print("Estimated Parameters:")
print("Means:", result['mu1'], result['mu2'])
print("Covariance Matrices:\n", result['Sigma1'], "\n", result['Sigma2'])
print("Scale Parameters:", result['b1'], result['b2'])
print("Negative Log-Likelihood:", result['nll'])

def calculate_aic_bic(nll, n_params, n_samples):
    """
    Calculate AIC and BIC for a given model.

    Parameters:
        nll: float
            Negative log-likelihood of the model.
        n_params: int
            Number of parameters in the model.
        n_samples: int
            Number of data samples.

    Returns:
        aic: float
            Akaike Information Criterion.
        bic: float
            Bayesian Information Criterion.
    """
    aic = 2 * n_params + 2 * nll
    bic = n_params * np.log(n_samples) + 2 * nll
    return aic, bic

# Number of parameters in the model (12 parameters)
n_params = 12
# Number of samples
n_samples = samples.shape[0]
# Negative log-likelihood
nll = result['nll']

# Calculate AIC and BIC
aic, bic = calculate_aic_bic(nll, n_params, n_samples)

# Display AIC and BIC
print("AIC:", aic)
print("BIC:", bic)

# # Create a 3D scatter plot
# fig = plt.figure(figsize=(10, 7))
# ax = fig.add_subplot(111, projection='3d')

# # Plot the samples
# ax.scatter(samples[:, 0], samples[:, 1], np.random.uniform(0, 1, size=len(samples)), c='blue', alpha=0.6)
# ax.set_title("3D Plot of 2-Component Bivariate Mixture of Laplace")
# ax.set_xlabel("X1")
# ax.set_ylabel("X2")
# ax.set_zlabel("Density")
# plt.show()

# # 2D Histogram (binning the data)
# bins = 20  # Number of bins in each dimension
# hist, xedges, yedges = np.histogram2d(samples[:, 0], samples[:, 1], bins=bins)

# # Coordinates of bin centers
# xpos, ypos = np.meshgrid(xedges[:-1] + np.diff(xedges) / 2, yedges[:-1] + np.diff(yedges) / 2, indexing="ij")
# xpos = xpos.ravel()
# ypos = ypos.ravel()
# zpos = np.zeros_like(xpos)

# # Heights of the bars (histogram counts)
# heights = hist.ravel()

# # Dimensions of the bars
# dx = dy = np.diff(xedges)[0]  # Width of each bin
# dz = heights  # Heights of the bars

# # Normalize heights for colormap
# norm = plt.Normalize(vmin=heights.min(), vmax=heights.max())
# colors = cm.viridis(norm(heights))  # Use the "viridis" colormap

# # Create the 3D bar plot
# fig = plt.figure(figsize=(10, 7))
# ax = fig.add_subplot(111, projection='3d')

# ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, alpha=0.7)

# # Set labels and title
# ax.set_title("3D Bar Plot of 2-Component Bivariate Mixture of Laplace")
# ax.set_xlabel("X1")
# ax.set_ylabel("X2")
# ax.set_zlabel("Frequency")

# # Add a colorbar for density
# sm = cm.ScalarMappable(cmap=cm.viridis, norm=norm)
# sm.set_array([])  # Dummy array for the colorbar
# plt.colorbar(sm, ax=ax, shrink=0.5, aspect=10, label="Density")

# plt.show()
# # Perform kernel density estimation
# kde = gaussian_kde(samples.T)

# # Create a grid of points in 2D space
# x1 = np.linspace(samples[:, 0].min(), samples[:, 0].max(), 100)
# x2 = np.linspace(samples[:, 1].min(), samples[:, 1].max(), 100)
# X1, X2 = np.meshgrid(x1, x2)
# positions = np.vstack([X1.ravel(), X2.ravel()])

# # Evaluate the density on the grid
# density = kde(positions).reshape(X1.shape)

# # Plot the density as a 3D surface plot
# fig = plt.figure(figsize=(10, 7))
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X1, X2, density, cmap='viridis', edgecolor='none')
# ax.set_title("3D Density Plot of 2-Component Bivariate Mixture of Laplace")
# ax.set_xlabel("X1")
# ax.set_ylabel("X2")
# ax.set_zlabel("Density")
# plt.show()


# # Compute the empirical CDF
# def empirical_cdf(samples, grid):
#     """
#     Compute the empirical CDF on a grid of points.
    
#     Parameters:
#         samples: ndarray
#             Sample points (n_samples, 2).
#         grid: ndarray
#             Grid points (n_grid_points, 2).

#     Returns:
#         cdf_values: ndarray
#             CDF values at the grid points.
#     """
#     n_samples = samples.shape[0]
#     cdf_values = np.zeros(grid.shape[0])
#     for i, point in enumerate(grid):
#         cdf_values[i] = np.sum(np.all(samples <= point, axis=1)) / n_samples
#     return cdf_values

# # Create a grid of points in 2D space
# x1 = np.linspace(samples[:, 0].min(), samples[:, 0].max(), 50)
# x2 = np.linspace(samples[:, 1].min(), samples[:, 1].max(), 50)
# X1, X2 = np.meshgrid(x1, x2)
# grid_points = np.vstack([X1.ravel(), X2.ravel()]).T

# # Evaluate the empirical CDF on the grid
# cdf_values = empirical_cdf(samples, grid_points).reshape(X1.shape)

# # Plot the CDF as a 3D surface plot
# fig = plt.figure(figsize=(10, 7))
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X1, X2, cdf_values, cmap='viridis', edgecolor='none')
# ax.set_title("3D CDF Plot of 2-Component Bivariate Mixture of Laplace")
# ax.set_xlabel("X1")
# ax.set_ylabel("X2")
# ax.set_zlabel("CDF")
# plt.show()





