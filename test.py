#%% Part 2 Prep

import numpy as np
from scipy.linalg import cho_solve, cho_factor
from scipy.special import gammaln
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.special import gammaln


import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.special import gammaln

def slog(x):
    # Truncated log to avoid -Inf or +Inf
    realmin = np.finfo(float).tiny
    realmax = np.finfo(float).max
    x_clipped = np.clip(x, realmin, realmax)
    return np.log(x_clipped)

def mvnctpdfln(x, mu, gam, v, Sigma):
    """
    Compute the log pdf of the multivariate noncentral t-distribution.

    Parameters
    ----------
    x : ndarray of shape (d, T)
        Evaluation points.
    mu : ndarray of shape (d,)
        Location vector.
    gam : ndarray of shape (d,)
        Noncentrality vector.
    v : float
        Degrees of freedom.
    Sigma : ndarray of shape (d, d)
        Dispersion matrix.

    Returns
    -------
    pdfLn : ndarray of shape (T,)
        Log pdf values at x.
    """
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




# Parameters for the multivariate noncentral t-distribution
mu = np.array([0, 0])        # Location vector (mean)
gam = np.array([0, 1])       # Noncentrality vector
v = 4                        # Degrees of freedom
Sigma = np.array([[1, 0.5],    # Covariance matrix
                  [0.5, 1]])

# Define a grid of points in 2D space
x1 = np.linspace(-5, 5, 100)  # Range for the first dimension
x2 = np.linspace(-5, 5, 100)  # Range for the second dimension
X1, X2 = np.meshgrid(x1, x2)  # Create a grid
x = np.vstack([X1.ravel(), X2.ravel()])  # Stack grid points into an Nx2 matrix

# Evaluate the log-PDF at each grid point
# Ensure you have the mvnctpdfln function defined in Python
pdfLn = mvnctpdfln(x, mu, gam, v, Sigma)
pdf = np.exp(pdfLn)
print(pdf)

# Reshape the output PDF to match the grid for plotting
pdf_grid = pdf.reshape(X1.shape)

# Plot the PDF as a contour plot
plt.figure(figsize=(10, 8))
contour = plt.contour(X1, X2, pdf_grid, levels=20)  # 20 contour levels
plt.colorbar(contour, label='PDF')
plt.title('PDF of Multivariate Noncentral t-Distribution')
plt.xlabel('x1')
plt.ylabel('x2')
plt.xlim(x1.min(), x1.max())
plt.ylim(x2.min(), x2.max())
plt.gca().set_aspect('equal', adjustable='box')  # Ensure the aspect ratio is equal
plt.show()
