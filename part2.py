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

def bivariate_laplace(mu, b, n):
    """
    Generate n samples from a bivariate Laplace distribution.
    Parameters:
        mu: [mu_x1, mu_x2] (location parameters)
        b: [b_x1, b_x2] (scale parameters)
        n: Number of samples
    Returns:
        samples: n x 2 array of samples
    """
    x1 = np.random.laplace(loc=mu[0], scale=b[0], size=n)
    x2 = np.random.laplace(loc=mu[1], scale=b[1], size=n)
    return np.column_stack((x1, x2))

def simulate_mixture(pi, mu1, b1, Sigma1, mu2, b2, Sigma2, n):
    """
    Simulate samples from a 2-component bivariate Laplace mixture.
    Parameters:
        pi: Mixing weight for the first component
        mu1: Location parameters for the first component
        b1: Scale parameters for the first component
        Sigma1: Covariance matrix for the first component
        mu2: Location parameters for the second component
        b2: Scale parameters for the second component
        Sigma2: Covariance matrix for the second component
        n: Number of samples
    Returns:
        samples: Simulated samples
    """
    # Generate samples for each component
    n1 = int(pi * n)
    n2 = n - n1
    
    # First component
    samples1 = bivariate_laplace(mu1, b1, n1)
    samples1 = samples1 @ np.linalg.cholesky(Sigma1).T  # Apply covariance
    
    # Second component
    samples2 = bivariate_laplace(mu2, b2, n2)
    samples2 = samples2 @ np.linalg.cholesky(Sigma2).T  # Apply covariance

    # Combine the samples
    return np.vstack((samples1, samples2))

pi = 0.7  # Mixing weight for the first component
mu1 = [0, 0]  # Location parameters for the first component
b1 = [10, 10]  # Scale parameters for the first component
Sigma1 = np.array([[1, 0.5], [0.5, 1]])  # Covariance matrix for the first component

mu2 = [0, 0]  # Location parameters for the second component
b2 = [5, 5]  # Scale parameters for the second component
Sigma2 = np.array([[4, 2], [2, 4]])  # Covariance matrix for the second component (more extreme returns)

n = 1000  # Number of samples

# Simulate the mixture
samples = simulate_mixture(pi, mu1, b1, Sigma1, mu2, b2, Sigma2, n)

print(samples)

# Create a 3D scatter plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Plot the samples
ax.scatter(samples[:, 0], samples[:, 1], np.random.uniform(0, 1, size=len(samples)), c='blue', alpha=0.6)
ax.set_title("3D Plot of 2-Component Bivariate Mixture of Laplace")
ax.set_xlabel("X1")
ax.set_ylabel("X2")
ax.set_zlabel("Density")
plt.show()


# a)

#Our definition
def modified_bessel_second_kind(nu, x):
    """
    Computes the modified Bessel function of the second kind K_nu(x) numerically.

    Parameters:
    - nu: Order of the Bessel function (ν > 0)
    - x: Argument of the Bessel function (x > 0)

    Returns:
    - K_nu(x): Value of the modified Bessel function of the second kind
    """
    if x < 0:
        raise ValueError("x must be positive.")
    if nu <= 0:
        raise ValueError("nu must be positive.")

    # Define the integrand
    def integrand(u):
        return 0.5 * (u**(nu - 1)) * np.exp(-x / 2 * (1 / u + u))

    # Compute the integral from 0 to infinity
    result, _ = quad(integrand, 0, np.inf, epsabs=1e-10, epsrel=1e-10)
    return result

# Example usage
nu = 1.5  # Order of the Bessel function
x = 2.0   # Argument of the Bessel function

K_nu_x = modified_bessel_second_kind(nu, x)
print(f"K_{nu}({x}) = {K_nu_x}", sp.kv(nu, x))


def graph_difference(nu, x_values):
    y_values = [abs(modified_bessel_second_kind(nu, x_val) - sp.kv(nu, x_val)) for x_val in tqdm(x_values)]
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, y_values, linestyle='-', color='blue', label='Difference')
    plt.title('Difference between bessel_k_int and scipy bessel_k')
    plt.xlabel('x values')
    plt.ylabel('Absolute Difference')
    plt.legend()
    plt.grid(True)
    plt.show()


# You can change the range, after 20 the difference is minimal
x_values = np.linspace(0, 10, 1000)
    
graph_difference(1.2, x_values)


import numpy as np
from scipy.special import gammaln
from scipy.optimize import minimize

def einschrk(pin, bound, Vin=None):
    lo = np.array(bound['low'])
    hi = np.array(bound['hi'])
    welche = np.array(bound['which'])
    pin = np.array(pin)
    if Vin is None or isinstance(Vin, (int, float)):
        # Vin is not given or is a placeholder
        trans = (hi + lo * pin ** 2) / (1 + pin ** 2)
        pout = (1 - welche) * pin + welche * trans
        Vout = None
    else:
        # Vin is given, adjust standard errors
        trans = (hi + lo * pin ** 2) / (1 + pin ** 2)
        pout = (1 - welche) * pin + welche * trans
        # Adjust the standard errors
        trans_derivative = 2 * pin * (lo - hi) / (1 + pin ** 2) ** 2
        d = (1 - welche) + welche * trans_derivative
        J = np.diag(d)
        Vout = J @ Vin @ J
    return pout, Vout

def mvtpdfmine(x, df, mu, Sigma):
    """
    Computes the density of a multivariate t-distribution.
    Input: x (d-dimensional array), df (degrees of freedom), mu (mean vector), Sigma (covariance matrix)
    """
    d = len(x)
    if mu is None or len(mu) == 0:
        mu = np.zeros(d)
    if Sigma is None or Sigma.size == 0:
        Sigma = np.eye(d)
    x = x.reshape(d)
    mu = mu.reshape(d)
    delta = x - mu
    term = delta.T @ np.linalg.solve(Sigma, delta)
    sign, log_det_Sigma = np.linalg.slogdet(Sigma)
    if sign <= 0:
        return 0
    log_pdf = gammaln((df + d)/2) - gammaln(df/2) - (d/2)*np.log(df * np.pi) - 0.5*log_det_Sigma \
              - ((df + d)/2)*np.log(1 + term / df)
    y = np.exp(log_pdf)
    return y

def MVTloglik(param, x, bound=None):
    if bound is not None:
        param, _ = einschrk(np.real(param), bound, Vin=999)
    nobs, d = x.shape
    Sig = np.zeros((d, d))
    k = param[0]
    mu = param[1:3]
    Sig[0, 0] = param[3]
    Sig[1, 1] = param[5]
    Sig[0, 1] = param[4]
    Sig[1, 0] = Sig[0, 1]
    if np.min(np.linalg.eigvals(Sig)) <= 0:
        return 1e5
    pdf = np.zeros(nobs)
    for i in range(nobs):
        pdf[i] = mvtpdfmine(x[i, :], k, mu, Sig)
        if pdf[i] == 0 or np.isnan(pdf[i]):
            return 1e5
    llvec = np.log(pdf)
    ll = -np.mean(llvec)
    if np.isinf(ll) or np.isnan(ll):
        ll = 1e5
    return ll

def MVTestimation(x):
    nobs, d = x.shape
    if d != 2:
        raise NotImplementedError('Not done yet, use EM for dimensions other than 2.')
    # Initialize parameters
    k = 2
    mu1 = 0
    mu2 = 0
    Sigma_11 = 1
    Sigma_12 = 0.5
    Sigma_22 = 1
    param = [k, mu1, mu2, Sigma_11, Sigma_12, Sigma_22]
    # Bounds for parameters
    bound = {
        'low': [0, -2, -2, 0.01, -90, 0.01],
        'hi': [90, 2, 2, 90, 90, 90],
        'which': [1, 1, 1, 1, 1, 1]
    }
    initvec = [2, -0.8, 0.2, 20, 2, 10]
    # Optimization options
    maxiter = 300
    tol = 1e-7
    options = {'disp': True, 'maxiter': maxiter, 'ftol': tol}
    # Define bounds for optimizer
    bounds = []
    for i in range(len(initvec)):
        if bound['which'][i]:
            bounds.append((bound['low'][i], bound['hi'][i]))
        else:
            bounds.append((None, None))
    # Optimization using scipy.optimize.minimize
    res = minimize(MVTloglik, initvec, args=(x, bound), method='L-BFGS-B', bounds=bounds, options=options)
    param = res.x
    fval = res.fun
    iters = res.nit
    hess_inv = res.hess_inv if hasattr(res, 'hess_inv') else None
    stderr = None
    if hess_inv is not None:
        if hasattr(hess_inv, 'todense'):
            hess_inv = hess_inv.todense()
        stderr = np.sqrt(np.diag(hess_inv))
    loglik = -fval * nobs
    return param, stderr, iters, loglik, hess_inv

# Example usage:
# x = np.random.multivariate_normal([0, 0], [[1, 0.5], [0.5, 1]], size=100)
# param, stderr, iters, loglik, hess_inv = MVTestimation(x)
# print("Estimated parameters:", param)
# print("Standard errors:", stderr)
# print("Iterations:", iters)
# print("Log-likelihood:", loglik)
