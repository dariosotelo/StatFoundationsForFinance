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
nonparametric_bootstrap_samples, var_nonpara = nonparametric_bootstrap(samples, 0.01, 100, 7500)
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
nonparametric_bootstrap_samples, var_nonpara = nonparametric_bootstrap(samples_t, 0.01, 100, 7500)
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
        seed = np.random.randint(0, 1234567)
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

stable_sample = stabgen(10000, 1.5)
nonparametric_bootstrap_samples, var_nonpara = nonparametric_bootstrap(stable_sample, 0.01, 100, 7500)
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

asymstab_samples = generate_asymstab_samples(10000, 1.5)
True_ES, True_VaR = asymstableES(0.01, 1.5)

nonparametric_bootstrap_samples, var_nonpara = nonparametric_bootstrap(asymstab_samples, 0.01, 100, 7500)
bootstrap_samples_t, var_t = parametric_bootstrap_t(asymstab_samples, 0.01, 100, 10000)
bootstrap_samples_nonc, var_nonc = parametric_bootstrap_nonc(asymstab_samples, 0.01, 100, 10000)
bootstrap_samples_gaussian, var_gaussian = parametric_bootstrap_gaussian(asymstab_samples, 0.01, 100, 10000)
bootstrap_samples_mixed, var_mixed = parametric_bootstrap_mixed(asymstab_samples, 0.01, 100, 10000)
parametric_bootstrap_GAt_samples, var_GAt = parametric_bootstrap_GAt(asymstab_samples, 0.01, 100, 10000)

data = {'Expected Shortfall': (nonparametric_bootstrap_samples + bootstrap_samples_t + bootstrap_samples_nonc + bootstrap_samples_gaussian + bootstrap_samples_mixed + parametric_bootstrap_GAt_samples), 'Type': (['Nonparametric'] * len(nonparametric_bootstrap_samples) + ['Student t'] * len(bootstrap_samples_t) + ['Noncentral t'] * len(bootstrap_samples_nonc) + ['Gaussian'] * len(bootstrap_samples_gaussian) + ['2-Comp Normal'] * len(bootstrap_samples_mixed) + ['GAt'] * len(parametric_bootstrap_GAt_samples))}
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

data_var = {'VaR': (var_nonpara + var_t + var_nonc + var_gaussian + var_mixed + var_GAt), 'Type': (['Nonparametric'] * len(var_nonpara) + ['Student t'] * len(var_t) + ['Noncentral t'] * len(var_nonc) + ['Gaussian'] * len(var_gaussian) + ['2-Comp Normal'] * len(var_mixed) + ['GAt'] * len(var_GAt))}
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