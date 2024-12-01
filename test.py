import numpy as np
from scipy.optimize import minimize
from scipy.special import kv, gamma
from scipy.linalg import inv, det, eigh


def bivariate_discrete_laplace_pdf(location, scale, Sigma, y):
    """
    Compute the PDF of the bivariate Laplace distribution.
    Parameters:
    - y: ndarray, observation vector of size (d,)
    - location: ndarray, location vector of size (d,)
    - Sigma: ndarray, positive-definite covariance matrix of size (d, d)
    - scale: float, parameter of the gamma distribution (must be > 0)
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
    pdf = normalization * ((bessel_arg / 2)**(scale / 2 - d / 4)) * bessel_factor

    return pdf


def log_likelihood(params, x, weights):
    """
    Compute the negative log-likelihood for the 2-component bivariate Laplace mixture.
    Parameters:
        params: Model parameters (12-dimensional vector)
        x: Data points (Nx2 matrix)
        weights: Mixture weights [w1, w2]
    Returns:
        neg_loglik: Negative log-likelihood value
    """
    # Extract parameters
    mu1 = params[0:2]
    mu2 = params[2:4]
    b1, b2 = params[4:6]
    sigma1 = np.array([[params[6], params[7]], [params[7], params[8]]])
    sigma2 = np.array([[params[9], params[10]], [params[10], params[11]]])
    
    # Ensure positive-definiteness of covariance matrices
    if np.min(eigh(sigma1, eigvals_only=True)) <= 0 or np.min(eigh(sigma2, eigvals_only=True)) <= 0:
        return 1e10

    # Log-likelihood for each component
    epsilon = 1e-10  # Small constant to avoid log(0)
    ll_comp1 = np.log(weights[0]) + np.log(np.array([bivariate_discrete_laplace_pdf(mu1, b1, sigma1, y) for y in x]) + epsilon)
    ll_comp2 = np.log(weights[1]) + np.log(np.array([bivariate_discrete_laplace_pdf(mu2, b2, sigma2, y) for y in x]) + epsilon)

    # Mixture log-likelihood
    ll_total = np.logaddexp(ll_comp1, ll_comp2)
    

    return -np.sum(ll_total)  # Negative log-likelihood


def estimate_mle(x, weights, init_params):
    """
    Estimate MLE for the 2-component bivariate Laplace mixture.
    Parameters:
        x: Data points (Nx2 matrix)
        weights: Mixture weights [w1, w2]
        init_params: Initial guess for the parameters (12-dimensional vector)
    Returns:
        result: Optimization result (MLE estimates, log-likelihood, standard errors)
    """
    bounds = [
        (-10, 10), (-10, 10),  # mu1
        (-10, 10), (-10, 10),  # mu2
        (-10, 10), (-10, 10),  # b1, b2
        (0.01, 10), (-0.99, 0.99), (0.01, 10),  # sigma1 (diag, off-diag, diag)
        (0.01, 10), (-0.99, 0.99), (0.01, 10)   # sigma2 (diag, off-diag, diag)
    ]

    result = minimize(
        log_likelihood,
        init_params,
        args=(x, weights),
        method='L-BFGS-B',
        bounds=bounds,
        options={'disp': True}
    )
    
    # Compute standard errors using the inverse Hessian
    hess_inv = result.hess_inv.todense() if hasattr(result.hess_inv, "todense") else result.hess_inv
    std_err = np.sqrt(np.diag(hess_inv)) if hess_inv is not None else None

    return {
        'params': result.x,
        'log_likelihood': -result.fun,
        'success': result.success,
        'std_err': std_err
    }


# Example usage
# Simulated data (Nx2)
n_samples = 2500
mu1 = np.array([1, 2])
mu2 = np.array([1, 0])
sigma1 = np.array([[1, 0.5], [0.5, 1]])
sigma2 = np.array([[4, 0.3], [0.3, 4]])
b1, b2 = 6, 2

# Generate data from two components
z = np.random.choice([0, 1], size=n_samples, p=[0.6, 0.4])  # Mixture indicator
x = np.array([
    np.random.multivariate_normal(mu1, b1 * sigma1) if z_i == 0 else
    np.random.multivariate_normal(mu2, b2 * sigma2)
    for z_i in z
])

# Mixture weights
weights = [0.6, 0.4]

# Initial parameter guess (12-dimensional vector)
init_params = [0, 0, 5, 5, 5, 1, 1, 0, 1, 2, 0, 2]

# Estimate MLE
result = estimate_mle(x, weights, init_params)

# Print results
print("Estimated Parameters:", result['params'])
print("Log-Likelihood:", result['log_likelihood'])
print("Standard Errors:", result['std_err'])
