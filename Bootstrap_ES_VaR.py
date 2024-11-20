# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Mon Nov 11 18:58:34 2024

# """

# #%% Libraries
# from scipy.special import betainc
# import numpy as np
# from scipy.optimize import fsolve
# from scipy.integrate import quad

# #%% Preparation functions

# #PDF
# def f_GAt(z, d, nu, theta, K):
#     if any(param <= 0 for param in [d, nu, theta]):
#         return "d, nu, or theta must be positive."
#     if z < 0:
#         return K * (1 + (-z * theta) ** d / nu) ** -(nu + 1/d)
#     else:
#         return K * (1 + (z / theta) ** d / nu) ** -(nu + 1/d)

# #Auxiliary function for the ES function
# def calculate_L(c, nu, theta, d):
#     return nu / (nu + (-c * theta) ** d)

# #ES
# def S_r(c, r, d, nu, theta):
#     if c >= 0:
#         return "c must be negative."
#     L = calculate_L(c, nu, theta, d)
#     B_L_num = betainc(nu - r / d, (r + 1) / d, L)  
#     B_L_den = betainc(nu, 1 / d, L)                
#     return ((-1) ** r) * (nu ** (r / d)) * ((1 + theta ** 2) / (theta ** r + theta ** (r + 2))) * (B_L_num / B_L_den)

# #The ES would be calculated using VaR as c and r=1

# #CDF
# def GAt(z, d, nu, theta, K=1):
#     pdf = f_GAt(z, d, nu, theta, K)
#     #Underscore because quad returns the error as the second parameter
#     cdf, _ = quad(lambda t: f_GAt(t, d, nu, theta, K), -np.inf, z)
#     return pdf, cdf


# #%% I.1


# # Paolella's code suggestion (it is a kinda revised copy-paste from ChatGPT)

# def GAtsim(sim, d, v, theta):
#     """
#     Simulates data from the GAt distribution.
    
#     Parameters:
#     sim (int): Number of simulations.
#     d (float): Parameter d of the GAt distribution.
#     v (float): Parameter v of the GAt distribution.
#     theta (float): Parameter theta of the GAt distribution.
    
#     Returns:
#     np.ndarray: Simulated data from the GAt distribution.
#     """
#     x = np.zeros(sim)
#     lo = 1e-6
#     hi = 1 - lo
    
#     for i in range(sim):
#         u = np.random.rand()
#         u = max(u, lo)
#         u = min(u, hi)
#         x[i] = GAtquantile(u, d, v, theta)
        
#     return x

# def GAtquantile(p, d, v, theta):
#     """
#     Computes the quantile of the GAt distribution for a given probability p.
    
#     Parameters:
#     p (float): Probability (0 < p < 1).
#     d (float): Parameter d of the GAt distribution.
#     v (float): Parameter v of the GAt distribution.
#     theta (float): Parameter theta of the GAt distribution.
    
#     Returns:
#     float: Quantile corresponding to the probability p.
#     """
#     # Find lower bound for the quantile
#     lobound = 0
#     while ff(lobound, p, d, v, theta) >= 0:
#         lobound -= 4
    
#     # Find upper bound for the quantile
#     hibound = 0
#     while ff(hibound, p, d, v, theta) <= 0:
#         hibound += 4
    
#     # Use fsolve to find the root of the function, i.e., the quantile
#     tol = 1e-5
#     q = fsolve(lambda x: ff(x, p, d, v, theta), [lobound, hibound], xtol=tol)[0]
    
#     return q





# def ff(x, p, d, v, theta, K=1):
#     """
#     Helper function to compute the difference between the CDF at x and the probability p.
    
#     Parameters:
#     x (float or array-like): Point at which to evaluate the CDF.
#     p (float): Target probability.
#     d (float): Shape parameter.
#     v (float): Shape parameter.
#     theta (float): Skew parameter.
    
#     Returns:
#     float: Difference between CDF(x) and p.
#     """
#     # If x is an array, extract the scalar value (fsolve may pass a one-element array)
#     if isinstance(x, np.ndarray):
#         x = x[0]
        
#     # Calculate the PDF and CDF at x
#     _, cdf = GAt(x, d, v, theta, K)
    
#     # Return the difference between CDF and p
#     return cdf - p


# #Something is not working but it is related to the array the GAtsim function returns which is then used in ff
# #I was thinking on doing a for loop but i havent revised Paolellas code in detail.


# # Set parameters
# sim = 5000
# d = 2.0
# v = 1.5
# theta = 0.9

# # Simulate data
# simulated_data = GAtsim(sim, d, v, theta)


# #%% Second request to ChatGPT, this one worked.


# import numpy as np
# from scipy.optimize import fsolve

# def GAtsim(sim, d, v, theta):
#     """
#     Simulates data from the GAt distribution.
    
#     Parameters:
#     sim (int): Number of simulations.
#     d (float): Parameter d of the GAt distribution.
#     v (float): Parameter v of the GAt distribution.
#     theta (float): Parameter theta of the GAt distribution.
    
#     Returns:
#     np.ndarray: Simulated data from the GAt distribution.
#     """
#     x = np.zeros(sim)
#     lo = 1e-6
#     hi = 1 - lo
    
#     for i in range(sim):
#         u = np.random.rand()
#         u = max(u, lo)
#         u = min(u, hi)
#         x[i] = GAtquantile(u, d, v, theta)
        
#     return x

# def GAtquantile(p, d, v, theta):
#     """
#     Computes the quantile of the GAt distribution for a given probability p.
    
#     Parameters:
#     p (float): Probability (0 < p < 1).
#     d (float): Parameter d of the GAt distribution.
#     v (float): Parameter v of the GAt distribution.
#     theta (float): Parameter theta of the GAt distribution.
    
#     Returns:
#     float: Quantile corresponding to the probability p.
#     """
#     # Find lower bound for the quantile
#     lobound = 0
#     while ff(lobound, p, d, v, theta) >= 0:
#         lobound -= 4
    
#     # Find upper bound for the quantile
#     hibound = 0
#     while ff(hibound, p, d, v, theta) <= 0:
#         hibound += 4
    
#     # Use fsolve to find the root of the function, i.e., the quantile
#     tol = 1e-5
#     q = fsolve(lambda x: ff(x, p, d, v, theta), x0=(lobound + hibound) / 2, xtol=tol)[0]
    
#     return q

# def ff(x, p, d, v, theta):
#     """
#     Helper function to compute the difference between the CDF at x and the probability p.
    
#     Parameters:
#     x (float): Point at which to evaluate the CDF.
#     p (float): Target probability.
#     d (float): Parameter d of the GAt distribution.
#     v (float): Parameter v of the GAt distribution.
#     theta (float): Parameter theta of the GAt distribution.
    
#     Returns:
#     float: Difference between CDF(x) and p.
#     """
#     _, cdf = GAt(x, d, v, theta)
#     return cdf - p

# def GAt(z, d, v, theta, K=1):
#     """
#     Computes the PDF and approximates the CDF of the GAt distribution.
    
#     Parameters:
#     z (float): Point at which to evaluate the PDF and CDF.
#     d (float): Shape parameter.
#     v (float): Shape parameter.
#     theta (float): Skew parameter.
#     K (float): Normalizing constant, if known.
    
#     Returns:
#     tuple: (pdf, cdf) values at point z.
#     """
#     # Calculate PDF using the GAt distribution formula
#     pdf = f_GAt(z, d, v, theta, K)
    
#     # Approximate CDF by integrating the PDF from -∞ to z
#     from scipy.integrate import quad
#     cdf, _ = quad(lambda t: f_GAt(t, d, v, theta, K), -np.inf, z)
    
#     return pdf, cdf

# def f_GAt(z, d, nu, theta, K=1):
#     """
#     PDF of the GAt distribution.
    
#     Parameters:
#     z (float): Point at which to evaluate the PDF.
#     d (float): Shape parameter.
#     nu (float): Shape parameter.
#     theta (float): Skew parameter.
#     K (float): Normalizing constant, if known.
    
#     Returns:
#     float: PDF value at z.
#     """
#     if any(param <= 0 for param in [d, nu, theta]):
#         raise ValueError("d, nu, or theta must be positive.")
    
#     if z < 0:
#         return K * (1 + (-z * theta) ** d / nu) ** -(nu + 1 / d)
#     else:
#         return K * (1 + (z / theta) ** d / nu) ** -(nu + 1 / d)


# sim = 5000
# d = 2.0
# v = 1.5
# theta = 0.9

# simulated_data = GAtsim(sim, d, v, theta)


# %%
from scipy.stats import norm, t, gaussian_kde, levy_stable
from scipy.special import beta, betainc, kv, gamma
from scipy.optimize import root_scalar, minimize, brentq
from scipy.fft import fftshift, ifft
from scipy.integrate import quad
from math import sqrt, exp
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from tqdm import tqdm
import random
import scipy.stats as stats

# Question 1 The following is the computation for GAt_pdf and GAt_cdf
def GAt_pdf(z, d, v, theta, loc=0, scale=1):
    z = (z - loc) / scale
    C_inv = (theta**-1 + theta) * d**-1 * v**(1/d) * beta(1/d, v)
    C = 1 / C_inv
    if isinstance(z, (float, int)):
        if z < 0:
            pdf = C * (1 + ((-z * theta)**d) / v)**(-(v + 1/d))
        else:
            pdf = C * (1 + ((z / theta)**d) / v)**(-(v + 1/d))
    else:
        pdf = np.where(z < 0, C * (1 + (np.abs(-z * theta)**d) / v)**(-(v + 1/d)), C * (1 + (np.abs(z / theta)**d) / v)**(-(v + 1/d)))
    pdf /= scale

    return pdf

def GAt_cdf(z, d, v, theta, loc=0, scale=1):
    z = (z - loc) / scale

    def B_L(a, b):
        return betainc(a, b, L)

    def B_U(a, b):
        return betainc(a, b, U)

    if isinstance(z, (float, int)):
        if z <= 0:
            L = v / (v + (-z * theta)**d)
            return B_L(v, 1/d) / (1 + theta**2)
        else:
            U = (z / theta)**d / (v + (z / theta)**d)
            return B_U(1/d, v) / (1 + theta**-2) + (1 + theta**2)**-1
    else:
        cdf = np.zeros_like(z)
        L = v / (v + (-z[z < 0] * theta)**d)
        cdf[z <= 0] = B_L(v, 1/d) / (1 + theta**2)
        U = (z[z > 0] / theta)**d / (v + (z / theta)**d)
        cdf[z > 0] = B_U(1/d, v) / (1 + theta**2)
    return cdf
    
def ff(x, p, d, v, theta, loc, scale):
    cdf = GAt_cdf(x, d, v, theta, loc, scale)
    return cdf - p

def GAtquantile(p, d, v, theta, loc=0, scale=1):
    lobound = loc
    while ff(lobound, p, d, v, theta, loc, scale) >= 0:
        lobound -= 4 * scale
    hibound = loc
    while ff(hibound, p, d, v, theta, loc, scale) <= 0:
        hibound += 4 * scale
    result = root_scalar(ff, args=(p, d, v, theta, loc, scale), bracket=[lobound, hibound], method='bisect', xtol=1e-5)
    if result.converged:
        return result.root
    else:
        raise ValueError("Root-finding did not converge")
    
sim = 10000
p = 0.01
d = 1.71
v = 1.96
theta = 0.879
loc = 0.175
scale = 1.21

def GAtsim(sim, d, v, theta, loc, scale):
    lo = 1e-6
    hi = 1 - lo
    u = np.random.uniform(low=lo, high=hi, size=sim)
    x = [GAtquantile(ui, d, v, theta, loc, scale) for ui in u]
    return x

quantile = GAtquantile(p, d, v, theta, loc, scale)
samples = GAtsim(sim, d, v, theta, loc, scale)

def conditional_expectation_Zr(r=1, c=quantile, d=1.71, v=1.96, theta=0.879, loc=0.175, scale=1.21):
    c = (c - loc) / scale
    L = v / (v + (-c * theta)**d)
    numerator = betainc(v - r/d, (r + 1)/d, L) * beta(v - r/d, (r + 1)/d)
    denominator = betainc(v, 1/d, L) * beta(v, 1/d)
    term1 = (-1)**r * v**(r/d)
    term2 = (1 + theta**2) / (theta**r + theta**(r + 2))
    expectation = term1 * term2 * numerator / denominator
    expectation *= scale**r
    
    return expectation

def nonparametric_bootstrap(data, ESlevel=0.01, B=500, n=10000):
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

def expected_shortfall(data, ESlevel=0.01):
    sorted_data = np.sort(data)
    index = int(ESlevel * len(data))
    ES = np.mean(sorted_data[:index])
    return ES

def parametric_bootstrap_t(data, ESlevel=0.01, B=500, n=10000):
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
    std_estimate = np.std(data, ddof=1)
    return mean_estimate, std_estimate

def parametric_bootstrap_gaussian(data, ESlevel=0.01, B=500, n=10000):
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
    def log_likelihood(params):
        nu, delta, loc, scale = params
        if nu <= 0 or scale <= 0:
            return -np.inf
        ll = np.sum(stats.nct.logpdf(data, nu, delta, loc=loc, scale=scale))
        return ll

    nu_init = 5.0
    delta_init = 0.0
    loc_init = 0
    scale_init = 1
    params_init = [nu_init, delta_init, loc_init, scale_init]
    bounds = [(1e-6, None), (None, None), (None, None),  (1e-6, None)]

    result = minimize(lambda params: -log_likelihood(params), params_init, bounds=bounds, method='L-BFGS-B')

    nu_mle, delta_mle, loc_mle, scale_mle = result.x
    return nu_mle, delta_mle, loc_mle, scale_mle

def parametric_bootstrap_nonc(data, ESlevel=0.01, B=500, n=10000):
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

def hartley(y, tol=1e-6, maxit=int(5e4)):
    init = np.array([np.mean(y) - np.std(y), np.mean(y) + np.std(y), np.std(y), np.std(y), 0.5])
    old = init.copy()
    new = np.zeros(5)
    iter = 0
    crit = 0
    y = np.asarray(y)
    while True:
        iter += 1
        mu1 = old[0]
        mu2 = old[1]
        s1 = old[2]
        s2 = old[3]
        lam = old[4]
        s1 = max(s1, 1e-6)
        s2 = max(s2, 1e-6)
        pdf1 = norm.pdf(y, mu1, s1)
        pdf2 = norm.pdf(y, mu2, s2)
        mixn = lam * pdf1 + (1 - lam) * pdf2
        epsilon = 1e-10
        mixn = np.maximum(mixn, epsilon)
        H1 = lam * pdf1 / mixn
        H2 = 1 - H1
        N1 = np.sum(H1)
        N2 = np.sum(H2)
        if N1 <= 0 or N2 <= 0:
            solvec = np.full(5, np.nan)
            break
        new[0] = np.sum(H1 * y) / N1
        new[1] = np.sum(H2 * y) / N2
        var1 = np.sum(H1 * (y - new[0])**2) / N1
        var2 = np.sum(H2 * (y - new[1])**2) / N2
        if var1 <= 0 or var2 <= 0:
            solvec = np.full(5, np.nan)
            break
        new[2] = np.sqrt(var1)
        new[3] = np.sqrt(var2)
        new[4] = np.mean(H1)
        crit = np.max(np.abs(old - new))
        solvec = new.copy()
        if np.isnan(solvec).any():
            break
        if crit < tol or iter >= maxit:
            break
        old = new.copy()
    mu_1, mu_2, sigma1, sigma2, pi1 = solvec
    return mu_1, mu_2, sigma1, sigma2, pi1

def generate_mixed_normal_samples(mu1, mu2, sigma1, sigma2, pi1, size=500):
    component_labels = np.random.choice([0, 1], size=size, p=[pi1, 1 - pi1])
    samples = np.where(component_labels == 0, np.random.normal(mu1, sigma1, size), np.random.normal(mu2, sigma2, size))
    
    return samples

def parametric_bootstrap_mixed(data, ESlevel=0.01, B=500, n=10000):
    mu1, mu2, sigma1, sigma2, pi1 = hartley(data)
    ES_samples = []
    VaR_samples = []
    for _ in range(B):
        bootstrap_sample = generate_mixed_normal_samples(mu1, mu2, sigma1, sigma2, pi1, size=n)
        ES = expected_shortfall(bootstrap_sample, ESlevel)
        VaR = np.percentile(bootstrap_sample, ESlevel * 100)
        ES_samples.append(ES)
        VaR_samples.append(VaR)
    return ES_samples, VaR_samples

def GAtestimation(x, vlo=None, initvec=None, fixd=None):
    if vlo is None:
        vlo = 0.01

    if initvec is None:
        initvec = []

    if fixd is None:
        fixd = []

    if not initvec:
        versuch = 3
        vhi = 4
        loglik = -np.inf
        vvec = np.linspace(vlo + 0.02, vhi, versuch)
        for vv in vvec:
            if not fixd:
                initvec = [2, vv, 0.98, 0, 3]
            else:
                initvec = [vv, 0.98, 0, 3]
            param0, stderr0, iters0, loglik0, Varcov0 = GAtestimation(x, vlo, initvec, fixd)
            if loglik0 > loglik:
                loglik = loglik0
                param = param0
                stderr = stderr0
                iters = iters0
                Varcov = Varcov0
        return param, stderr, iters, loglik, Varcov

    if not fixd:
        bound_lo = [0.1, vlo, 0.2, -1, 1e-4]
        bound_hi = [30, 100, 3, 2, 1e+4]
        bounds = [(lo, hi) for lo, hi in zip(bound_lo, bound_hi)]
    else:
        bound_lo = [vlo, 0.2, -1, 1e-4]
        bound_hi = [100, 3, 2, 1e+4]
        bounds = [(lo, hi) for lo, hi in zip(bound_lo, bound_hi)]
    nobs = len(x)
    maxiter = len(initvec) * 100
    tol = 1e-8
    opts = {'disp': False, 'maxiter': maxiter, 'gtol': tol}
    res = minimize(lambda param: GAtloglik(param, x, fixd), initvec, method='L-BFGS-B', bounds=bounds, options=opts)
    pout = res.x
    fval = res.fun
    theoutput = res
    hess_inv = res.hess_inv.todense() if hasattr(res.hess_inv, "todense") else res.hess_inv
    V = hess_inv / nobs
    param = pout
    Varcov = V
    stderr = np.sqrt(np.diag(V))
    loglik = -fval * nobs
    iters = theoutput.nit

    return param, stderr, iters, loglik, Varcov

def GAtloglik(param, x, fixd):
    if fixd is None or len(fixd) == 0:
        d, v, theta, mu, c = param
    else:
        d = fixd
        v, theta, mu, c = param
    z = (x - mu) / c
    pdf = GAt_pdf(z, d, v, theta) / c
    llvec = np.log(pdf + 1e-10)
    ll = -np.mean(llvec)
    if np.isinf(ll) or np.isnan(ll):
        ll = 1e5

    return ll

# print(GAtestimation(samples)[0])

def parametric_bootstrap_GAt(data, ESlevel=0.01, B=500, n=10000):
    param, _, _, _, _ = GAtestimation(data)
    ES_samples = []
    VaR_samples = []
    for k in tqdm(range(B)):
        bootstrap_sample = GAtsim(n, param[0], param[1], param[2], param[3], param[4])
        ES = expected_shortfall(bootstrap_sample, ESlevel)
        VaR = np.percentile(bootstrap_sample, ESlevel * 100)
        ES_samples.append(ES)
        VaR_samples.append(VaR)
    return ES_samples, VaR_samples


# Generate bootstrap samples
nonparametric_bootstrap_samples, var_nonpara = nonparametric_bootstrap(samples, 0.01, 100, 6667)
bootstrap_samples_t, var_t = parametric_bootstrap_t(samples, 0.01, 100, 10000)
bootstrap_samples_nonc, var_nonc = parametric_bootstrap_nonc(samples, 0.01, 100, 10000)
bootstrap_samples_gaussian, var_gaussian = parametric_bootstrap_gaussian(samples, 0.01, 100, 10000)
bootstrap_samples_mixed, var_mixed = parametric_bootstrap_mixed(samples, 0.01, 100, 10000)
parametric_bootstrap_GAt_samples, var_GAt = parametric_bootstrap_GAt(samples, 0.01, 100, 10000)

original_ES = conditional_expectation_Zr(1, quantile, d=1.71, v=1.96, theta=0.879, loc=0.175, scale=1.21)

data = {'Expected Shortfall': (nonparametric_bootstrap_samples + bootstrap_samples_t + bootstrap_samples_nonc + bootstrap_samples_gaussian + bootstrap_samples_mixed + parametric_bootstrap_GAt_samples), 'Type': (['Nonparametric'] * len(nonparametric_bootstrap_samples) + ['Student t'] * len(bootstrap_samples_t) + ['Noncentral t'] * len(bootstrap_samples_nonc) + ['Gaussian'] * len(bootstrap_samples_gaussian) + ['2-Comp Normal'] * len(bootstrap_samples_mixed) + ['GAt'] * len(parametric_bootstrap_GAt_samples))}
df = pd.DataFrame(data)

plt.figure(figsize=(10, 6))
sns.boxplot(x='Type', y='Expected Shortfall', data=df, flierprops=dict(marker='+', markeredgecolor='r', markersize=10), boxprops=dict(facecolor='none', edgecolor='blue'))
plt.axhline(y=original_ES, color='g', linestyle='--', label='Original ES')
plt.title('Boxplot of Expected Shortfall for DGP = GAt')
plt.xlabel('Type')
plt.ylabel('Expected Shortfall')
plt.legend()
plt.grid(True)
plt.show()

data_var = {'VaR': (var_nonpara + var_t + var_nonc + var_gaussian + var_mixed + var_GAt), 'Type': (['Nonparametric'] * len(var_nonpara) + ['Student t'] * len(var_t) + ['Noncentral t'] * len(var_nonc) + ['Gaussian'] * len(var_gaussian) + ['2-Comp Normal'] * len(var_mixed) + ['GAt'] * len(var_GAt))}
df_var = pd.DataFrame(data_var)

plt.figure(figsize=(10, 6))
sns.boxplot(x='Type', y='VaR', data=df_var, flierprops=dict(marker='+', markeredgecolor='r', markersize=10), boxprops=dict(facecolor='none', edgecolor='blue'))
plt.axhline(y=quantile, color='g', linestyle='--', label='True VaR')
plt.title('Boxplot of VaR for DGP = GAt')
plt.xlabel('Type')
plt.ylabel('VaR')
plt.legend()
plt.grid(True)
plt.show()

samples_t = np.random.standard_t(4, 10000)

def student_t_ES(df, loc=0, scale=1, ESlevel=0.01):
    return -stats.t.pdf(stats.t.ppf(ESlevel, df), df, loc, scale)/stats.t.cdf(stats.t.ppf(ESlevel, df), df, loc, scale)*(df+stats.t.ppf(ESlevel, df)**2)/(df-1)

def VaRt(df, alpha=0.01):
    return stats.t.ppf(alpha, df)

# Generate bootstrap samples
nonparametric_bootstrap_samples, var_nonpara = nonparametric_bootstrap(samples_t, 0.01, 100, 6667)
bootstrap_samples_t, var_t = parametric_bootstrap_t(samples_t, 0.01, 100, 10000)
bootstrap_samples_nonc, var_nonc = parametric_bootstrap_nonc(samples_t, 0.01, 100, 10000)
bootstrap_samples_gaussian, var_gaussian = parametric_bootstrap_gaussian(samples_t, 0.01, 100, 10000)
bootstrap_samples_mixed, var_mixed = parametric_bootstrap_mixed(samples_t, 0.01, 100, 10000)
parametric_bootstrap_GAt_samples, var_GAt = parametric_bootstrap_GAt(samples_t, 0.01, 100, 10000)

original_ES = student_t_ES(4)

data = {'Expected Shortfall': (nonparametric_bootstrap_samples + bootstrap_samples_t + bootstrap_samples_nonc + bootstrap_samples_gaussian + bootstrap_samples_mixed + parametric_bootstrap_GAt_samples), 'Type': (['Nonparametric'] * len(nonparametric_bootstrap_samples) + ['Student t'] * len(bootstrap_samples_t) + ['Noncentral t'] * len(bootstrap_samples_nonc) + ['Gaussian'] * len(bootstrap_samples_gaussian) + ['2-Comp Normal'] * len(bootstrap_samples_mixed) + ['GAt'] * len(parametric_bootstrap_GAt_samples))}
df = pd.DataFrame(data)

plt.figure(figsize=(10, 6))
sns.boxplot(x='Type', y='Expected Shortfall', data=df, flierprops=dict(marker='+', markeredgecolor='r', markersize=10), boxprops=dict(facecolor='none', edgecolor='blue'))
plt.axhline(y=original_ES, color='g', linestyle='--', label='Original ES')
plt.title('Boxplot of Expected Shortfall for DGP = Student t')
plt.xlabel('Type')
plt.ylabel('Expected Shortfall')
plt.legend()
plt.grid(True)
plt.show()

data_var = {'VaR': (var_nonpara + var_t + var_nonc + var_gaussian + var_mixed + var_GAt), 'Type': (['Nonparametric'] * len(var_nonpara) + ['Student t'] * len(var_t) + ['Noncentral t'] * len(var_nonc) + ['Gaussian'] * len(var_gaussian) + ['2-Comp Normal'] * len(var_mixed) + ['GAt'] * len(var_GAt))}
df_var = pd.DataFrame(data_var)

plt.figure(figsize=(10, 6))
sns.boxplot(x='Type', y='VaR', data=df_var, flierprops=dict(marker='+', markeredgecolor='r', markersize=10), boxprops=dict(facecolor='none', edgecolor='blue'))
plt.axhline(y=VaRt(4), color='g', linestyle='--', label='True VaR')
plt.title('Boxplot of VaR for DGP = Student t')
plt.xlabel('Type')
plt.ylabel('VaR')
plt.legend()
plt.grid(True)
plt.show()

def stabgen(nobs, a, b=0, c=1, d=0, seed=None):
    if seed is None:
        seed = np.random.randint(0, 9999999)
    rng_V = np.random.default_rng(seed)
    rng_W = np.random.default_rng(seed + 42)
    
    V = rng_V.uniform(-np.pi / 2, np.pi / 2, size=nobs)
    W = rng_W.exponential(scale=1.0, size=nobs)
    
    if a == 1:
        x = (2 / np.pi) * (((np.pi / 2) + b * V) * np.tan(V) - b * np.log((W * np.cos(V)) / ((np.pi / 2) + b * V)))
        x = c * x + d - (2 / np.pi) * d * np.log(d) * c + b
    else:
        Cab = np.arctan(b * np.tan(np.pi * a / 2)) / a
        Sab = (1 + b**2 * (np.tan(np.pi * a / 2))**2)**(1 / (2 * a))
        A = np.sin(a * (V + Cab)) / (np.cos(V))**(1 / a)
        B0 = np.cos(V - a * (V + Cab)) / W
        B = np.abs(B0)**((1 - a) / a)
        x = Sab * A * B * np.sign(B0)
        x = c * x + d
    return x

stable_sample = stabgen(10000, 1.7, -0.3)
nonparametric_bootstrap_samples, var_nonpara = nonparametric_bootstrap(stable_sample, 0.01, 100, 6667)
bootstrap_samples_t, var_t = parametric_bootstrap_t(stable_sample, 0.01, 100, 10000)
bootstrap_samples_nonc, var_nonc = parametric_bootstrap_nonc(stable_sample, 0.01, 100, 10000)
bootstrap_samples_gaussian, var_gaussian = parametric_bootstrap_gaussian(stable_sample, 0.01, 100, 10000)
bootstrap_samples_mixed, var_mixed = parametric_bootstrap_mixed(stable_sample, 0.01, 100, 10000)
parametric_bootstrap_GAt_samples, var_GAt = parametric_bootstrap_GAt(stable_sample, 0.01, 100, 10000)
#Check if this is correct - ES and VaR
True_ES = expected_shortfall(stable_sample, 0.01)
True_VaR = np.percentile(stable_sample, 1)

data = {'Expected Shortfall': (nonparametric_bootstrap_samples + bootstrap_samples_t + bootstrap_samples_nonc + bootstrap_samples_gaussian + bootstrap_samples_mixed + parametric_bootstrap_GAt_samples), 'Type': (['Nonparametric'] * len(nonparametric_bootstrap_samples) + ['Student t'] * len(bootstrap_samples_t) + ['Noncentral t'] * len(bootstrap_samples_nonc) + ['Gaussian'] * len(bootstrap_samples_gaussian) + ['2-Comp Normal'] * len(bootstrap_samples_mixed) + ['GAt'] * len(parametric_bootstrap_GAt_samples))}
df = pd.DataFrame(data)

plt.figure(figsize=(10, 6))
sns.boxplot(x='Type', y='Expected Shortfall', data=df, flierprops=dict(marker='+', markeredgecolor='r', markersize=10), boxprops=dict(facecolor='none', edgecolor='blue'))
plt.axhline(y=True_ES, color='g', linestyle='--', label='Original ES')
plt.title('Boxplot of Expected Shortfall for DGP = Stable Paretian')
plt.xlabel('Type')
plt.ylabel('Expected Shortfall')
plt.legend()
plt.grid(True)
plt.show()

data_var = {'VaR': (var_nonpara + var_t + var_nonc + var_gaussian + var_mixed + var_GAt), 'Type': (['Nonparametric'] * len(var_nonpara) + ['Student t'] * len(var_t) + ['Noncentral t'] * len(var_nonc) + ['Gaussian'] * len(var_gaussian) + ['2-Comp Normal'] * len(var_mixed) + ['GAt'] * len(var_GAt))}
df_var = pd.DataFrame(data_var)

plt.figure(figsize=(10, 6))
sns.boxplot(x='Type', y='VaR', data=df_var, flierprops=dict(marker='+', markeredgecolor='r', markersize=10), boxprops=dict(facecolor='none', edgecolor='blue'))
plt.axhline(y=True_VaR, color='g', linestyle='--', label='True VaR')
plt.title('Boxplot of VaR for DGP = Stable Paretian')
plt.xlabel('Type')
plt.ylabel('VaR')
plt.legend()
plt.grid(True)
plt.show()

def asymstableES(xi, a, b=0, mu=0, scale=1, method=1):
    alpha = a
    beta = b
    def stabcdfroot(x):
        F = levy_stable.cdf(x, alpha, beta)
        return F - xi
    try:
        q = brentq(stabcdfroot, -1e6, 1e6, xtol=1e-6)
    except ValueError as e:
        raise ValueError(f"Error finding quantile q: {e}")

    VaR = mu + scale * q
    if q == 0:
        t0 = (1 / a) * np.arctan(b * np.tan(np.pi * a / 2))
        ES = ((2 * gamma((a - 1) / a)) / (np.pi - 2 * t0)) * \
             (np.cos(t0) / np.cos(a * t0) ** (1 / a))
        ES = ES * scale + mu
        return ES, VaR

    if method == 1:
        method = 2
    if method == 2:
        tailcomp = stabletailcomp(q, a, b)
        ES = (scale * tailcomp / xi) + mu
        return ES, VaR

def stabletailcomp(q, a, b):

    K = (a / np.pi) * np.sin(np.pi * a / 2) * gamma(a) * (1 - b)
    ell = -1e6
    if a == 1:
        term1 = 0
    else:
        term1 = K * (-ell) ** (1 - a) / (1 - a)
    def stableCVARint(x):
        den = asymstabpdf(x, a, b)
        return x * den
    try:
        term3, _ = quad(stableCVARint, ell, q, limit=100, epsabs=1e-5)
    except Exception as e:
        return None

    tailcomp = term1 + term3
    return tailcomp

def asymstabpdf(x, a, b):
    pdf = levy_stable.pdf(x, a, b)
    return pdf

def generate_asymstab_samples(n, alpha, beta=0, mu=0, scale=1):
    u = np.random.uniform(0, 1, n)
    samples = np.zeros(n)
    for i in tqdm(range(n)):
        samples[i] = brentq(lambda x: levy_stable.cdf(x, alpha, beta) - u[i], -1e6, 1e6)
    return mu + scale * samples

def mle_asymmetric_stable(data):
    def neg_log_likelihood(params):
        alpha, beta, loc, scale = params
        if not (0 < alpha <= 2 and -1 <= beta <= 1 and scale > 0):
            return np.inf
        return -np.sum(levy_stable.logpdf(data, alpha, beta, loc=loc, scale=scale))

    initial_params = [1.7, -0.3, np.mean(data), np.std(data)]
    bounds = [(0.01, 2), (-1, 1), (None, None), (1e-6, None)]
    result = minimize(neg_log_likelihood, initial_params, bounds=bounds, method='L-BFGS-B')
    return result.x

def parametric_bootstrap_asymstab(data, ESlevel=0.01, B=100, n=2500):
    alpha_mle, beta_mle, loc_mle, scale_mle = mle_asymmetric_stable(data)
    ES_samples = []
    VaR_samples = []
    for k in tqdm(range(B)):
        bootstrap_sample = generate_asymstab_samples(n, alpha_mle, beta_mle, loc_mle, scale_mle)
        ES = expected_shortfall(bootstrap_sample, ESlevel)
        VaR = np.percentile(bootstrap_sample, ESlevel * 100)
        ES_samples.append(ES)
        VaR_samples.append(VaR)
    return ES_samples, VaR_samples

asymstab_samples = generate_asymstab_samples(10000, 1.7, -0.3)
True_ES, True_VaR = asymstableES(0.01, 1.7, -0.3)

nonparametric_bootstrap_samples, var_nonpara = nonparametric_bootstrap(asymstab_samples, 0.01, 100, 6667)
bootstrap_samples_t, var_t = parametric_bootstrap_t(asymstab_samples, 0.01, 100, 10000)
bootstrap_samples_nonc, var_nonc = parametric_bootstrap_nonc(asymstab_samples, 0.01, 100, 10000)
bootstrap_samples_gaussian, var_gaussian = parametric_bootstrap_gaussian(asymstab_samples, 0.01, 100, 10000)
bootstrap_samples_mixed, var_mixed = parametric_bootstrap_mixed(asymstab_samples, 0.01, 100, 10000)
parametric_bootstrap_GAt_samples, var_GAt = parametric_bootstrap_GAt(asymstab_samples, 0.01, 100, 10000)
boostrap_sample_asymstab, var_asymstab = parametric_bootstrap_asymstab(asymstab_samples, 0.01, 100, 10000)

data = {'Expected Shortfall': (nonparametric_bootstrap_samples + bootstrap_samples_t + bootstrap_samples_nonc + bootstrap_samples_gaussian + bootstrap_samples_mixed + parametric_bootstrap_GAt_samples + boostrap_sample_asymstab), 'Type': (['Nonparametric'] * len(nonparametric_bootstrap_samples) + ['Student t'] * len(bootstrap_samples_t) + ['Noncentral t'] * len(bootstrap_samples_nonc) + ['Gaussian'] * len(bootstrap_samples_gaussian) + ['2-Comp Normal'] * len(bootstrap_samples_mixed) + ['GAt'] * len(parametric_bootstrap_GAt_samples) + ['Asymmetric Stable'] * len(boostrap_sample_asymstab))}
df = pd.DataFrame(data)

plt.figure(figsize=(10, 6))
sns.boxplot(x='Type', y='Expected Shortfall', data=df, flierprops=dict(marker='+', markeredgecolor='r', markersize=10), boxprops=dict(facecolor='none', edgecolor='blue'))
plt.axhline(y=True_ES, color='g', linestyle='--', label='Original ES')
plt.title('Boxplot of Expected Shortfall for DGP = Asymmetric Stable')
plt.xlabel('Type')
plt.ylabel('Expected Shortfall')
plt.legend()
plt.grid(True)
plt.show()

data_var = {'VaR': (var_nonpara + var_t + var_nonc + var_gaussian + var_mixed + var_GAt + var_asymstab), 'Type': (['Nonparametric'] * len(var_nonpara) + ['Student t'] * len(var_t) + ['Noncentral t'] * len(var_nonc) + ['Gaussian'] * len(var_gaussian) + ['2-Comp Normal'] * len(var_mixed) + ['GAt'] * len(var_GAt) + ['Asymmetric Stable'] * len(var_asymstab))}
df_var = pd.DataFrame(data_var)

plt.figure(figsize=(10, 6))
sns.boxplot(x='Type', y='VaR', data=df_var, flierprops=dict(marker='+', markeredgecolor='r', markersize=10), boxprops=dict(facecolor='none', edgecolor='blue'))
plt.axhline(y=True_VaR, color='g', linestyle='--', label='True VaR')
plt.title('Boxplot of VaR for DGP = Asymmetric Stable')
plt.xlabel('Type')
plt.ylabel('VaR')
plt.legend()
plt.grid(True)
plt.show()




#%% Part 2

import numpy as np
from scipy.linalg import cho_solve, cho_factor
from scipy.special import gammaln
import matplotlib.pyplot as plt

def mvnctpdf_ln(x, mu, gam, v, Sigma):
    """
    Direct density approximation (d.d.a.) for the log of the d-variate canonical MVNCT density.

    Parameters:
    x : ndarray
        d x T matrix of evaluation points (each column is a point to evaluate).
    mu : ndarray
        d-length location vector.
    gam : ndarray
        d-length noncentrality vector.
    v : float
        Degrees of freedom.
    Sigma : ndarray
        Dispersion matrix (d x d covariance-like matrix).

    Returns:
    pdf_ln : ndarray
        Logarithm of the probability density function at the evaluation points.
    """
    # Dimensions
    d, T = x.shape

    # Cholesky decomposition of Sigma
    try:
        C, lower = cho_factor(Sigma, lower=True, check_finite=True)
    except np.linalg.LinAlgError:
        raise ValueError("Sigma must be positive definite.")
    
    # Reshape inputs
    mu = mu.reshape(-1, 1) if mu.ndim == 1 else mu
    gam = gam.reshape(-1, 1) if gam.ndim == 1 else gam
    vn2 = (v + d) / 2

    # Center x and compute normalized x
    xm = x - mu
    xm = np.linalg.solve(C, xm)  # Cholesky solve

    # Compute initial terms
    rho = np.sum((xm - gam) ** 2, axis=0)
    pdf_ln = gammaln(vn2) - (d / 2) * np.log(np.pi * v) - gammaln(v / 2) \
             - np.sum(np.log(np.diag(C))) - 0.5 * vn2 * np.log(1 + rho / v)

    # Return if gam is zero
    if np.all(gam == 0):
        return pdf_ln

    # Initialize variables
    idx = (pdf_ln > -37)
    maxiter = 10000
    k = 0
    logsumk = np.zeros_like(pdf_ln)

    # Iterate for approximation
    while np.any(idx) and k < maxiter:
        # Compute terms
        gcg = np.sum((np.linalg.solve(C, gam)) ** 2)
        term = 0.5 * np.log(2) + np.log(v + k) - 0.5 * np.log(v + rho) - 0.5 * gcg
        logterms = gammaln((v + d + k) / 2) - gammaln((k + 1) / 2) - gammaln(vn2) + k * term
        ff = np.exp(logterms[idx])
        
        # Update sum
        logsumk[idx] = np.logaddexp(logsumk[idx], np.log(ff))
        
        # Check convergence
        if np.all(np.abs(ff) < 1e-4):
            break

        k += 1

    # Add the final sum to the PDF log values
    pdf_ln += logsumk

    return pdf_ln

# Truncated logarithm function
def slog(x):
    """
    Truncated logarithm to avoid numerical issues with -inf or +inf.
    """
    realmin = np.finfo(np.float64).tiny
    realmax = np.finfo(np.float64).max
    return np.log(np.clip(x, realmin, realmax))



#%% Test for the plots

import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import quad
import scipy.special as sp
from scipy.optimize import minimize

def plot_mvnct_density(mu, gam, v, Sigma, title):
    # Define grid for x1 and x2
    x1 = np.linspace(-10, 10, 100)
    x2 = np.linspace(-10, 10, 100)
    X1, X2 = np.meshgrid(x1, x2)
    grid = np.array([X1.ravel(), X2.ravel()])

    # Reshape inputs to column vectors
    mu = np.reshape(mu, (-1, 1))
    gam = np.reshape(gam, (-1, 1))

    # Debugging: Print dimensions of inputs
    print(f"mu shape: {mu.shape}, gam shape: {gam.shape}, Sigma shape: {Sigma.shape}")

    # Ensure Sigma is positive definite
    if not np.all(np.linalg.eigvals(Sigma) > 0):
        raise ValueError("Sigma must be positive definite.")

    # Compute the log density at each grid point
    log_density = mvnctpdf_ln(grid, mu, gam, v, Sigma)
    density = np.exp(log_density).reshape(X1.shape)  # Convert log density to density

    # Plot contour
    plt.figure(figsize=(6, 6))
    plt.contour(X1, X2, density, levels=15, cmap="viridis")
    plt.title(title)
    plt.xlabel(r"$X_1$")
    plt.ylabel(r"$X_2$")
    plt.colorbar(label="Density")
    plt.grid(True)
    plt.show()

# Parameters for each plot
v = 4
Sigma_identity = np.eye(2)

# First plot: mu = [0, 0], gam = [0, 0], Sigma = I2
mu1 = np.zeros(2)
gam1 = np.zeros(2)
plot_mvnct_density(mu1, gam1, v, Sigma_identity, "MVNCT: v=4, γ=[0, 0], Σ=I2")

# Second plot: mu = [0, 1], gam = [0, 1], Sigma = I2
mu2 = np.array([0, 1])
gam2 = np.array([0, 1])
plot_mvnct_density(mu2, gam2, v, Sigma_identity, "MVNCT: v=4, γ=[0, 1], Σ=I2")

# Third plot: mu = [0, 1], gam = [0, 1], R = 0.5
R = np.array([[1, 0.5], [0.5, 1]])  # Correlation matrix
plot_mvnct_density(mu2, gam2, v, R, "MVNCT: v=4, γ=[0, 1], R=0.5")


def mvnctpdf_ln(x, mu, gam, v, Sigma):
    d, T = x.shape
    C = np.linalg.cholesky(Sigma)

    mu = mu.reshape(-1, 1) if mu.ndim == 1 else mu
    gam = gam.reshape(-1, 1) if gam.ndim == 1 else gam
    vn2 = (v + d) / 2

    xm = x - mu
    xm = np.linalg.solve(C, xm)
    rho = np.sum((xm - gam) ** 2, axis=0)
    pdf_ln = gammaln(vn2) - (d / 2) * np.log(np.pi * v) - gammaln(v / 2) \
             - np.sum(np.log(np.diag(C))) - 0.5 * vn2 * np.log(1 + rho / v)

    if np.all(gam == 0):
        return pdf_ln

    idx = (pdf_ln > -37)
    maxiter = 1000
    k = 0
    logsumk = np.zeros_like(pdf_ln)

    while np.any(idx) and k < maxiter:
        gcg = np.sum((np.linalg.solve(C, gam)) ** 2)
        term = 0.5 * np.log(2) + np.log(v + k) - 0.5 * np.log(v + rho) - 0.5 * gcg
        logterms = gammaln((v + d + k) / 2) - gammaln((k + 1) / 2) - gammaln(vn2) + k * term
        logterms = np.clip(logterms, -700, 700)  # Avoid numerical overflow
        ff = np.exp(logterms[idx])
        if np.all(np.abs(ff) < 1e-4):
            break

        logsumk[idx] = np.logaddexp(logsumk[idx], np.log(ff))
        k += 1

    pdf_ln += logsumk
    print("Eigenvalues of Sigma:", np.linalg.eigvals(Sigma))

    print(f"Iteration {k}, max logterms: {np.max(logterms)}, min logterms: {np.min(logterms)}")
    return pdf_ln

x1 = np.linspace(-5, 5, 100)
x2 = np.linspace(-5, 5, 100)



# First plot: mu = [0, 0], gam = [0, 0], Sigma = I2
mu1 = np.zeros(2)
gam1 = np.zeros(2)
plot_mvnct_density(mu1, gam1, v, Sigma_identity, "MVNCT: v=4, γ=[0, 0], Σ=I2")

# Second plot: mu = [0, 1], gam = [0, 1], Sigma = I2
mu2 = np.array([0, 1])
gam2 = np.array([0, 1])
plot_mvnct_density(mu2, gam2, v, Sigma_identity, "MVNCT: v=4, γ=[0, 1], Σ=I2")

# Third plot: mu = [0, 1], gam = [0, 1], R = 0.5
R = np.array([[1, 0.5], [0.5, 1]])  # Correlation matrix
plot_mvnct_density(mu2, gam2, v, R, "MVNCT: v=4, γ=[0, 1], R=0.5")



#%% II.1

import numpy as np
import matplotlib.pyplot as plt

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


# Parameters
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

#%% II.2

# a)

#Our definition
def bessel_k_int(nu, x):
    integrand = lambda u: 0.5*u**(nu-1)*np.exp(-x/2*(1/u+u))
    result, _ = quad(integrand, 0, np.inf)
    return result

# sp.kv

# Test
nu = 1.5
x = 2.0
sp_result = sp.kv(nu, x)
our_result = bessel_k_int(nu, x)

def graph_difference(nu, x, x_values):
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
x_values = list(range(-100, 20))
    
graph_difference(1.2, 2.0, x_values)
    

# In this part idk if he wants us to report a table w the differences (part a)


# b) 


# I think we should check this definition of the function
def bivariate_discrete_laplace_logpdf(mean, scale, x):
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
def negative_log_likelihood(params, x):
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

# Generate synthetic data (mixture of two bivariate discrete Laplace distributions)
np.random.seed(42)
n_samples = 200
data1 = np.random.randint(-5, 5, size=(n_samples // 2, 2))
data2 = np.random.randint(5, 15, size=(n_samples // 2, 2))
data = np.vstack([data1, data2])

# Initial guesses for parameters
initial_guess = [
    0.5,  # w1 (weight of the first component)
    0, 0,  # mu1 (mean of the first component)
    1, 1,  # b1 (scale of the first component)
    10, 10,  # mu2 (mean of the second component)
    1, 1  # b2 (scale of the second component)
]

# Set parameter bounds (to mimic MATLAB constraints)
bounds = [
    (0.01, 0.99),  # w1 (weights must sum to 1)
    (-10, 10),     # mu1_1
    (-10, 10),     # mu1_2
    (0.1, 10),     # b1_1
    (0.1, 10),     # b1_2
    (-10, 20),     # mu2_1
    (-10, 20),     # mu2_2
    (0.1, 10),     # b2_1
    (0.1, 10)      # b2_2
]

# Perform optimization using BFGS
result = minimize(
    fun=negative_log_likelihood,
    x0=initial_guess,
    args=(data,),
    method='L-BFGS-B',
    bounds=bounds,
    options={'disp': True}
)

# Extract results
estimated_params = result.x
print("Estimated parameters:", estimated_params)

# Optional: Extract standard errors using the Hessian inverse
hessian_inv = result.hess_inv.todense() if hasattr(result.hess_inv, "todense") else result.hess_inv
standard_errors = np.sqrt(np.diag(hessian_inv)) if hessian_inv is not None else None
print("Standard errors:", standard_errors)


# This part is to test if the function is approximating correctly the distribution
def generate_bivariate_discrete_laplace(n, w1, mu1, b1, mu2, b2):
    """
    Generate synthetic data from a k=2 bivariate discrete mixture of Laplace.
    Args:
        n: Number of samples.
        w1: Weight of the first component.
        mu1: Mean vector of the first component.
        b1: Scale parameters of the first component.
        mu2: Mean vector of the second component.
        b2: Scale parameters of the second component.
    Returns:
        Synthetic bivariate data (n x 2 array).
    """
    data = np.zeros((n, 2))
    for i in range(n):
        # Randomly choose component based on weights
        if np.random.rand() < w1:
            # Component 1
            data[i, 0] = np.random.laplace(mu1[0], b1[0])
            data[i, 1] = np.random.laplace(mu1[1], b1[1])
        else:
            # Component 2
            data[i, 0] = np.random.laplace(mu2[0], b2[0])
            data[i, 1] = np.random.laplace(mu2[1], b2[1])
    return data

# True parameters
w1_true = 0.6
mu1_true = [1, 2]
b1_true = [1, 1.5]
mu2_true = [5, 6]
b2_true = [1.2, 0.8]

# Generate synthetic data
n_samples = 1000
data = generate_bivariate_discrete_laplace(
    n=n_samples, w1=w1_true, mu1=mu1_true, b1=b1_true, mu2=mu2_true, b2=b2_true
)

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


# This is the important optimization part
# Perform optimization using the earlier code

result = minimize(
    fun=negative_log_likelihood,
    x0=initial_guess,
    args=(data,),
    method='L-BFGS-B',
    bounds=bounds,
    options={'disp': True}
)

# Extract results
estimated_params = result.x
print("\nEstimated Parameters:")
print(f"w1: {estimated_params[0]}")
print(f"mu1: {estimated_params[1:3]}")
print(f"b1: {estimated_params[3:5]}")
print(f"mu2: {estimated_params[5:7]}")
print(f"b2: {estimated_params[7:9]}")

# Compare with true parameters
print("\nTrue Parameters vs Estimated Parameters:")
print(f"True w1: {w1_true}, Estimated w1: {estimated_params[0]}")
print(f"True mu1: {mu1_true}, Estimated mu1: {estimated_params[1:3]}")
print(f"True b1: {b1_true}, Estimated b1: {estimated_params[3:5]}")
print(f"True mu2: {mu2_true}, Estimated mu2: {estimated_params[5:7]}")
print(f"True b2: {b2_true}, Estimated b2: {estimated_params[7:9]}")

# Calculate differences
print("\nDifferences:")
print(f"w1 difference: {abs(w1_true - estimated_params[0])}")
print(f"mu1 difference: {np.abs(np.array(mu1_true) - estimated_params[1:3])}")
print(f"b1 difference: {np.abs(np.array(b1_true) - estimated_params[3:5])}")
print(f"mu2 difference: {np.abs(np.array(mu2_true) - estimated_params[5:7])}")
print(f"b2 difference: {np.abs(np.array(b2_true) - estimated_params[7:9])}")


# It is not a good estimation, something must be checked, maybe our x_0 should be closer to the real values?














