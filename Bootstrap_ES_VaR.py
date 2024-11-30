
# %% Part 1
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

stable_sample = stabgen(10000, 0, -0.3)
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
    samples = levy_stable.ppf(u, alpha, beta, mu, scale)
    return samples

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

def parametric_bootstrap_asymstab(data, ESlevel=0.01, B=100, n=10000):
    alpha_mle, beta_mle, loc_mle, scale_mle = mle_asymmetric_stable(data)
    print(alpha_mle, beta_mle, loc_mle, scale_mle)
    ES_samples = []
    VaR_samples = []
    for _ in tqdm(range(B)):
        bootstrap_sample = stabgen(n, alpha_mle, beta_mle, scale_mle, loc_mle)
        ES = expected_shortfall(bootstrap_sample, ESlevel)
        VaR = np.percentile(bootstrap_sample, ESlevel * 100)
        ES_samples.append(ES)
        VaR_samples.append(VaR)
    return ES_samples, VaR_samples


asymstab_samples = stabgen(10000, 1.7, -0.3)
True_ES, True_VaR = asymstableES(0.01, 1.7, -0.3)
boostrap_sample_asymstab, var_asymstab = parametric_bootstrap_asymstab(asymstab_samples, 0.01, 100, 10000)

data = {'Expected Shortfall': (boostrap_sample_asymstab), 'Type': (['AsymStab'] * len(boostrap_sample_asymstab))}
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

data_var = {'VaR': (var_asymstab), 'Type': (['AsymStab'] * len(var_asymstab))}
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


nonparametric_bootstrap_samples, var_nonpara = nonparametric_bootstrap(asymstab_samples, 0.01, 100, 6667)
bootstrap_samples_t, var_t = parametric_bootstrap_t(asymstab_samples, 0.01, 100, 10000)
bootstrap_samples_nonc, var_nonc = parametric_bootstrap_nonc(asymstab_samples, 0.01, 100, 10000)
bootstrap_samples_gaussian, var_gaussian = parametric_bootstrap_gaussian(asymstab_samples, 0.01, 100, 10000)
bootstrap_samples_mixed, var_mixed = parametric_bootstrap_mixed(asymstab_samples, 0.01, 100, 10000)
parametric_bootstrap_GAt_samples, var_GAt = parametric_bootstrap_GAt(asymstab_samples, 0.01, 100, 10000)



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


#%% Part 2 Prep
import numpy as np
from scipy.linalg import cho_solve, cho_factor
from scipy.special import gammaln
import matplotlib.pyplot as plt
from scipy.linalg import cholesky, solve_triangular

def slog(x):
    # Truncated log to avoid -Inf or +Inf
    realmin = np.finfo(float).tiny
    realmax = np.finfo(float).max
    x_clipped = np.clip(x, realmin, realmax)
    return np.log(x_clipped)

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
plt.title('PDF of Multivariate Noncentral t-Distribution (mu=[0,0], gam=[0,1], v=4, Sigma=[[1,0.5],[0.5,1]])')
plt.xlabel('x1')
plt.ylabel('x2')
plt.xlim(x1.min(), x1.max())
plt.ylim(x2.min(), x2.max())
plt.gca().set_aspect('equal', adjustable='box')  # Ensure the aspect ratio is equal
plt.show()
# Parameters for the first plot
mu1 = [0, 0]
gam1 = [0, 0]
v1 = 4
Sigma1 = [[1, 0], [0, 1]]

# Evaluate the log-PDF at each grid point for the first set of parameters
pdfLn1 = mvnctpdfln(x, mu1, gam1, v1, Sigma1)
pdf1 = np.exp(pdfLn1)

# Reshape the output PDF to match the grid for plotting
pdf_grid1 = pdf1.reshape(X1.shape)

# Plot the PDF as a contour plot for the first set of parameters
plt.figure(figsize=(10, 8))
contour1 = plt.contour(X1, X2, pdf_grid1, levels=20)  # 20 contour levels
plt.colorbar(contour1, label='PDF')
plt.title('PDF of Multivariate Noncentral t-Distribution (mu=[0,0], gam=[0,0], v=4, Sigma=[[1,0],[0,1]])')
plt.xlabel('x1')
plt.ylabel('x2')
plt.xlim(x1.min(), x1.max())
plt.ylim(x2.min(), x2.max())
plt.gca().set_aspect('equal', adjustable='box')  # Ensure the aspect ratio is equal
plt.show()

# Parameters for the second plot
mu2 = [0, 0]
gam2 = [0, 1]
v2 = 4
Sigma2 = [[1, 0], [0, 1]]

# Evaluate the log-PDF at each grid point for the second set of parameters
pdfLn2 = mvnctpdfln(x, mu2, gam2, v2, Sigma2)
pdf2 = np.exp(pdfLn2)

# Reshape the output PDF to match the grid for plotting
pdf_grid2 = pdf2.reshape(X1.shape)

# Plot the PDF as a contour plot for the second set of parameters
plt.figure(figsize=(10, 8))
contour2 = plt.contour(X1, X2, pdf_grid2, levels=20)  # 20 contour levels
plt.colorbar(contour2, label='PDF')
plt.title('PDF of Multivariate Noncentral t-Distribution (mu=[0,0], gam=[0,1], v=4, Sigma=[[1,0],[0,1]])')
plt.xlabel('x1')
plt.ylabel('x2')
plt.xlim(x1.min(), x1.max())
plt.ylim(x2.min(), x2.max())
plt.gca().set_aspect('equal', adjustable='box')  # Ensure the aspect ratio is equal
plt.show()

#%% II.1
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm  # For colormap
from scipy.integrate import quad
import scipy.special as sp
from scipy.optimize import minimize
from scipy.special import gammaln, gamma, kv
from scipy.linalg import cholesky

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

def simulate_mixture_bivariate_laplace(pi, mu1, b1, Sigma1, mu2, b2, Sigma2, n):
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
samples = simulate_mixture_bivariate_laplace(pi, mu1, b1, Sigma1, mu2, b2, Sigma2, n)

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
# Perform kernel density estimation
kde = gaussian_kde(samples.T)

# Create a grid of points in 2D space
x1 = np.linspace(samples[:, 0].min(), samples[:, 0].max(), 100)
x2 = np.linspace(samples[:, 1].min(), samples[:, 1].max(), 100)
X1, X2 = np.meshgrid(x1, x2)
positions = np.vstack([X1.ravel(), X2.ravel()])

# Evaluate the density on the grid
density = kde(positions).reshape(X1.shape)

# Plot the density as a 3D surface plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X1, X2, density, cmap='viridis', edgecolor='none')
ax.set_title("3D Density Plot of 2-Component Bivariate Mixture of Laplace")
ax.set_xlabel("X1")
ax.set_ylabel("X2")
ax.set_zlabel("Density")
plt.show()
# CDF Plot for the Mixture of Bivariate Laplace and Bivariate Noncentral t-Distribution

# Define a function to compute the CDF from the PDF
def compute_cdf(pdf, x1, x2):
    cdf = np.zeros_like(pdf)
    for i in range(pdf.shape[0]):
        for j in range(pdf.shape[1]):
            cdf[i, j] = np.sum(pdf[:i+1, :j+1])
    cdf /= cdf[-1, -1]  # Normalize to ensure the CDF ranges from 0 to 1
    return cdf

# Compute the CDF for the mixture of bivariate Laplace
cdf_mixture_laplace = compute_cdf(pdf_grid, x1, x2)

# Plot the CDF as a 3D surface plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X1, X2, cdf_mixture_laplace, cmap='viridis', edgecolor='none')
ax.set_title("3D CDF Plot of 2-Component Bivariate Mixture of Laplace")
ax.set_xlabel("X1")
ax.set_ylabel("X2")
ax.set_zlabel("CDF")
plt.show()

#%% II.2

# a)

#Our definition
def bessel_k_int(nu, x):
    integrand = lambda u: 0.5*u**(nu-1)*np.exp(-x/2*(1/u+u))
    result, _ = quad(integrand, 0, np.inf)
    return result

# # Test
# nu = 1.5
# x = 2.0
# sp_result = sp.kv(nu, x)
# our_result = bessel_k_int(nu, x)

def graph_difference(nu, x_values):
    y_values = [abs(bessel_k_int(nu, x_val) - sp.kv(nu, x_val)) for x_val in x_values]
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, y_values, linestyle='-', color='blue', label='Difference')
    plt.title('Difference between bessel_k_int and scipy bessel_k')
    plt.xlabel('x values')
    plt.ylabel('Absolute Difference')
    plt.legend()
    plt.grid(True)
    plt.show()


# You can change the range, after 20 the difference is minimal
x_values = np.linspace(0, 20, 200)
    
graph_difference(1.2, x_values)
    

# In this part idk if he wants us to report a table w the differences (part a)


# b) 

# This part i am using it later on
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
    print("X:",x)
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
    fun=negative_log_likelihood_bvlp,
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
    fun=negative_log_likelihood_bvlp,
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
# this is not finished.

#%% II.3

# # Im not using this one afaik
# def negative_log_likelihood(params, x):
#     """
#     Negative log-likelihood for the bivariate NCT distribution.
    
#     Parameters:
#     params : ndarray
#         Parameter vector: [mu1, mu2, gam1, gam2, v, Sigma_11, Sigma_12, Sigma_22].
#     x : ndarray
#         T x 2 dataset (bivariate observations).
        
#     Returns:
#     float
#         Negative log-likelihood value.
#     """
#     # Extract parameters
#     mu = np.array([params[0], params[1]])  # Location vector
#     gam = np.array([params[2], params[3]])  # Noncentrality vector
#     v = params[4]  # Degrees of freedom
#     Sigma = np.array([
#         [params[5], params[6]],  # Covariance matrix
#         [params[6], params[7]]
#     ])

#     """  linear algebra requirements check, there might be something wrong w this  
#     if not np.all(np.isfinite(Sigma)):
#         print("Invalid Sigma (contains NaNs or infs):", Sigma)
#         return np.inf  # Penalize invalid Sigma

#     # Ensure Sigma is positive definite
#     if np.any(np.linalg.eigvals(Sigma) <= 0):
#         print("Non-positive definite Sigma:", Sigma)
#         return np.inf
    

#      """   
#     # Call the provided log-density function
    
#     try:
#         log_pdf = mvnctpdfln(x, mu, gam, v, Sigma)
#     except ValueError as e:
#         # Handle any errors (e.g., non-positive-definite Sigma)
#         return np.inf
    
#     # Return the negative log-likelihood
#     return -np.sum(log_pdf)


# def compute_mle(data_set, x_0):
    
#     # Minimize the negative log-likelihood
#     result = minimize(
#         fun=negative_log_likelihood,
#         x0=x_0,
#         args=(data_set,),
#         method='L-BFGS-B',
#         options={'disp': False}  # Display optimization details
#     )
    
#     # Check for successful optimization
#     if not result.success:
#         raise RuntimeError(f"Optimization failed: {result.message}")
    
#     # Return the estimated parameters
#     return result.x

# x_0 = [
#     0.0, 0.0,  # mu1, mu2
#     0.0, 0.0,  # gam1, gam2
#     5.0,        # v (degrees of freedom)
#     1.0, 0.0, 1.0  # Sigma_11, Sigma_12, Sigma_22
# ]

# np.random.seed(42)
# x = np.random.multivariate_normal(
#     mean=[1.0, 2.0],
#     cov=[[1.0, 0.5], [0.5, 1.5]],
#     size=200
# )

# mle_params = compute_mle(x.T, x_0)
# print("sadsfljkajdsflkajdfs")
# print("Estimated Parameters:", mle_params)
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import inv, eigvals
from scipy.stats import chi2
from scipy.linalg import cholesky
from tqdm import tqdm

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
    """
    Estimates the parameters of a bivariate noncentral t-distribution.

    Parameters
    ----------
    x : ndarray of shape (T, 2)
        Input data, where T is the number of observations.

    Returns
    -------
    param : ndarray
        Estimated parameters.
    stderr : ndarray
        Standard errors of the estimated parameters.
    iters : int
        Number of iterations performed by the optimizer.
    loglik : float
        Log-likelihood of the estimated model.
    Varcov : ndarray
        Variance-covariance matrix of the estimated parameters.
    """
    x = np.asarray(x)
    T, d = x.shape  # Now the data is (T, 2) instead of (2, T)
    if d != 2:
        raise ValueError('Not implemented for dimensions other than 2.')

    # Transpose the data for compatibility with the rest of the code
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





# # Parameters for the distribution
# mu = np.array([0, 0])           # Location vector
# gam = np.array([0, 2.5])          # Noncentrality vector
# v = 4                           # Degrees of freedom
# Sigma = np.array([[1, 0.5],
#                   [0.5, 1]])     # Covariance matrix
# n_samples = 500              # Number of samples

# x_min, x_max = -15, 15  # Range for both dimensions
# y_min, y_max = -15, 15

# # Precompute the maximum PDF value (optional for efficiency)
# x_test = np.linspace(x_min, x_max, 100)
# y_test = np.linspace(y_min, y_max, 100)
# X_test, Y_test = np.meshgrid(x_test, y_test)
# test_points = np.vstack([X_test.ravel(), Y_test.ravel()])
# log_pdf_test = mvnctpdfln(test_points, mu, gam, v, Sigma)
# pdf_test = np.exp(log_pdf_test)
# pdf_max = np.max(pdf_test)  # Maximum value of the PDF

# # Initialize storage for samples
# samples = np.zeros((2, n_samples))
# count = 0

# # Rejection sampling loop
# while count < (n_samples):
#     # Step 1: Generate a random point in the sampling space
#     x_rand = x_min + (x_max - x_min) * np.random.rand()
#     y_rand = y_min + (y_max - y_min) * np.random.rand()
#     candidate = np.array([x_rand, y_rand])

#     # Step 2: Evaluate the PDF at the candidate point
#     log_pdf_val = mvnctpdfln(candidate.reshape(2, 1), mu, gam, v, Sigma)
#     pdf_val = np.exp(log_pdf_val)[0]  # Since mvnctpdfln returns an array

#     # Step 3: Generate a uniform random number and accept/reject
#     u = np.random.rand() * pdf_max  # Scale uniform random number by maximum PDF value
#     if u <= pdf_val:
#         samples[:, count] = candidate
#         count += 1

# print(samples)

# actual = [v, mu[0], mu[1], Sigma[0][0], Sigma[1][1], Sigma[0][1], gam[0], gam[1]]
#print('Actual Parameters:')
# print(actual)

# Generate random 2 by 1000 data points
np.random.seed(42)
samples = np.random.randn(1000, 2)
print(samples)

# Assuming you have already implemented the mvnctpdfln and MVNCT2estimation functions
param, stderr, iters, loglik, Varcov = MVNCT2estimation(samples)
print('Estimated Parameters:')
print(param)
print('Standard Errors:')
print(stderr)

# Compute AIC and BIC for the MLE parameters
def compute_aic_bic(log_likelihood, num_params, num_data):
    aic = -2 * log_likelihood + 2 * num_params
    bic = -2 * log_likelihood + num_params * np.log(num_data)
    return aic, bic

# Number of parameters estimated
num_params = len(param)

# Number of data points
num_data = samples.shape[1]

# Compute AIC and BIC
aic, bic = compute_aic_bic(loglik, num_params, num_data)

print(f"AIC: {aic}")
print(f"BIC: {bic}")

#This version is a not robust one of the previous
#%%


# Program 1
def negative_log_likelihood_bvnct(params, x):
    mu = np.array([params[0], params[1]])  # Location vector
    gam = np.array([params[2], params[3]])  # Noncentrality vector
    v = params[4]  # Degrees of freedom
    Sigma = np.array([
        [params[5], params[6]],  # Covariance matrix
        [params[6], params[7]]
    ])    
    log_pdf = mvnctpdfln(x, mu, gam, v, Sigma)
    
    return -np.sum(log_pdf)


def compute_mle(data_set, x_0):    
    # Minimize the negative log-likelihood
    result = minimize(
        fun=negative_log_likelihood_bvnct,
        x0=x_0,
        args=(data_set,),
        method='L-BFGS-B',
        options={'disp': False}
    )
    
    if not result.success:
        raise RuntimeError(f"Optimization failed: {result.message}")
    
    return result.x




x_0 = [
    0.0, 0.0,  # mu1, mu2
    0.0, 0.0,  # gam1, gam2
    5.0,        # v (degrees of freedom)
    1.0, 0.0, 1.0  # Sigma_11, Sigma_12, Sigma_22
    
]

# I tried using a generated sample but it doesnt work, im so fkn tired

np.random.seed(42)
x = np.random.multivariate_normal(
    mean=[1.0, 2.0],
    cov=[[1.0, 0.5], [0.5, 1.5]],
    size=200
)

mle_params = compute_mle(x.T, x_0)
print("sadsfljkajdsflkajdfs")
print("Estimated Parameters:", mle_params)


#%% Test of the previous code but for a Gaussian sample

import numpy as np
from scipy.optimize import minimize
from scipy.stats import multivariate_normal

# Define the negative log-likelihood for the bivariate Gaussian distribution
def negative_log_likelihood(params, x):
    mu = np.array([params[0], params[1]])  # Mean vector
    Sigma = np.array([
        [params[2], params[3]],  # Covariance matrix
        [params[3], params[4]]
    ])
    
    # Check if the covariance matrix is valid (positive semi-definite)
    # This part is not included in the previous code
    if np.linalg.det(Sigma) <= 0:
        return np.inf  # Return a high value to penalize invalid covariance matrices

    # Compute the log likelihood using scipy's multivariate_normal
    log_pdf = multivariate_normal.logpdf(x, mean=mu, cov=Sigma)
    return -np.sum(log_pdf)  # Negative log-likelihood

# Define a function to compute the MLE
def compute_mle_example(data_set, x_0):    
    # Minimize the negative log-likelihood
    result = minimize(
        fun=negative_log_likelihood,
        x0=x_0,
        args=(data_set,),
        method='L-BFGS-B',
        options={'disp': False}
    )
    
    if not result.success:
        raise RuntimeError(f"Optimization failed: {result.message}")
    
    return result.x  # Return the estimated parameters

# Initial parameters for Gaussian MLE
# mu1, mu2, Sigma_11, Sigma_12, Sigma_22
x_0 = [
    0.0, 0.0,  # Initial guess for mean
    1.0, 0.0, 1.0  # Initial guess for covariance matrix
]

# Generate synthetic bivariate Gaussian data
np.random.seed(42)
x = np.random.multivariate_normal(
    mean=[1.0, 2.0],
    cov=[[1.0, 0.5], [0.5, 1.5]],
    size=200
)

# Compute the MLE
mle_params = compute_mle(x, x_0)

# Display the results
print("Estimated Parameters:")
print("Mean (mu1, mu2):", mle_params[:2])
print("Covariance Matrix:")
print(f"[[{mle_params[2]:.4f}, {mle_params[3]:.4f}]")
print(f" [{mle_params[3]:.4f}, {mle_params[4]:.4f}]]")

#%% Test for P1

#This is a function to sample a dataset of mvnct
def sample_mvnct(mu, Sigma, v, size=1):
    
    mu = np.asarray(mu)  # Ensure mu is a NumPy array
    d = len(mu)  # Dimensionality of the distribution
    
    # Step 1: Sample from a standard multivariate normal distribution
    Z = np.random.multivariate_normal(mean=np.zeros(d), cov=Sigma, size=size)
    
    # Step 2: Sample from a chi-squared distribution with v degrees of freedom
    W = np.random.chisquare(df=v, size=size)
    
    # Step 3: Transform Z and W to generate MVNCT samples
    # Each row of Z corresponds to a sample from MVNCT
    samples = mu + Z / np.sqrt(W[:, None] / v)
    
    return samples
x_0 = [
    1.0, 0.0,  # mu1, mu2
    0.0, 1.0,  # gam1, gam2
    6.0,        # v (degrees of freedom)
    1.0, 0.0, 1.0 ] # Sigma_11, Sigma_12, Sigma_22
    


    
mu = [0, 0]  # Mean vector (non-centrality)
Sigma = [[1, 0.5], [0.5, 1]]  # Covariance matrix
v = 5  # Degrees of freedom
n_samples = 1000  # Number of samples to generate

# Generate samples
np.random.seed(42)
samples = sample_mvnct(mu, Sigma, v, size=n_samples)

mle_params = compute_mle(samples, x_0)

# Display the results
print("Estimated Parameters:")
print("Mean (mu1, mu2):", mle_params[:2])
print("Covariance Matrix:")
print(f"[[{mle_params[2]:.4f}, {mle_params[3]:.4f}]")
print(f" [{mle_params[3]:.4f}, {mle_params[4]:.4f}]]")



#%% 

# Program 2 This is not program 2, the good code is in the next cell.

import numpy as np
from scipy.optimize import minimize
from scipy.stats import laplace, nct
from scipy.special import gammaln


# Safe log to handle numerical issues
def safe_log(x):
    return np.log(np.maximum(x, np.finfo(float).eps))


# Mixture Laplace Model
# I dont think i usedthis function
def neg_log_likelihood_mixture_laplace(params, data):
    pi, mu1_x, mu1_y, mu2_x, mu2_y, b1, b2 = params[:7]
    pi = 1 / (1 + np.exp(-pi))  # Map pi to (0, 1) using sigmoid
    mu1 = np.array([mu1_x, mu1_y])
    mu2 = np.array([mu2_x, mu2_y])
    b1, b2 = np.exp(b1), np.exp(b2)  # Ensure scales are positive

    # Log-likelihood for each component
    ll1 = laplace.logpdf(data, loc=mu1, scale=b1).sum(axis=1)
    ll2 = laplace.logpdf(data, loc=mu2, scale=b2).sum(axis=1)

    # Mixture log-likelihood
    ll = safe_log(pi * np.exp(ll1) + (1 - pi) * np.exp(ll2))
    return -np.sum(ll)


# Mixture Laplace from above


def bivariate_laplace(mu, b, n):
    x1 = np.random.laplace(loc=mu[0], scale=b[0], size=n)
    x2 = np.random.laplace(loc=mu[1], scale=b[1], size=n)
    return np.column_stack((x1, x2))

def simulate_mixture(pi, mu1, b1, Sigma1, mu2, b2, Sigma2, n):
    n1 = int(pi * n)
    n2 = n - n1
    
    samples1 = bivariate_laplace(mu1, b1, n1)
    samples1 = samples1 @ np.linalg.cholesky(Sigma1).T  # Apply covariance
    
    samples2 = bivariate_laplace(mu2, b2, n2)
    samples2 = samples2 @ np.linalg.cholesky(Sigma2).T  # Apply covariance

    return np.vstack((samples1, samples2))


# Bivariate NCT Model
# Neither this one
def neg_log_likelihood_bivariate_nct(params, data):
    mu_x, mu_y, gam_x, gam_y, v, sigma_x, sigma_y, rho = params
    mu = np.array([mu_x, mu_y])
    gam = np.array([gam_x, gam_y])
    v = np.exp(v)  # Ensure degrees of freedom are positive
    Sigma = np.array([[sigma_x**2, rho * sigma_x * sigma_y],
                      [rho * sigma_x * sigma_y, sigma_y**2]])

    # Log-likelihood using the multivariate NCT
    def mnvtcpdfln(x, mu, gam, v, Sigma):
        d, t = x.shape
        C = Sigma
        R = np.linalg.cholesky(C).T
        vn2 = (v + d) / 2
        xm = x - mu.reshape(-1, 1)
        rho = np.sum(np.linalg.solve(R.T, xm)**2, axis=0)
        pdfln = (gammaln(vn2) 
                 - (d / 2) * np.log(np.pi * v) 
                 - np.sum(np.log(np.diag(R))) 
                 - vn2 * safe_log(v + rho))
        return pdfln.sum()

    log_likelihood = mnvtcpdfln(data.T, mu, gam, v, Sigma)
    return -log_likelihood


# AIC and BIC Calculation
def compute_aic_bic(log_likelihood, num_params, num_data):
    aic = -2 * log_likelihood + 2 * num_params
    bic = -2 * log_likelihood + num_params * np.log(num_data)
    return aic, bic


# MLE Computation Function
def compute_mle(data, model, initial_guess):
    # Choose the likelihood function
    if model == 'mixture_laplace':
        neg_log_likelihood = neg_log_likelihood_mixture_laplace
    elif model == 'bivariate_nct':
        neg_log_likelihood = neg_log_likelihood_bivariate_nct
    else:
        raise ValueError("Invalid model type.")

    # Minimize the negative log-likelihood
    result = minimize(
        fun=neg_log_likelihood,
        x0=initial_guess,
        args=(data,),
        method='L-BFGS-B',
        options={'disp': True}
    )

    # Check for optimization success
    if not result.success:
        raise RuntimeError(f"Optimization failed: {result.message}")

    # Return the estimated parameters and log-likelihood
    return result.x, -result.fun


#%%

# This part should work when we fix the neg_log_likelihood of the bvnct
# Q3 part 2
def aic_bic_bvlp_bvnct(data, initial_guess):
         
    initial_guess_mlp = initial_guess[0]
    result_mixture_laplace = minimize(
    fun=negative_log_likelihood_bvlp,
    x0=initial_guess_mlp,
    args=(data,),  # Pass the dataset to the function
    method="L-BFGS-B"  # Robust method for bounds and large problems
    )
    
    initial_guess_bvnct = initial_guess[1]
    result_bvnct = minimize(
    fun=negative_log_likelihood_bvnct,
    x0=initial_guess_bvnct,
    args=(data,),  # Pass the dataset to the function
    method="L-BFGS-B"  # Robust method for bounds and large problems
    )
    
    log_likelihood_mixture_laplace = -result_mixture_laplace.fun
    log_likelihood_bvnct = -result_bvnct.fun
    
    # Number of parameters (dimensionality of initial_guess)
    k_mlp = len(initial_guess_mlp)
    k_bvnct = len(initial_guess_bvnct)
    
    # Number of data points
    n = len(data)
    
    aic_mixture_laplace = 2 * k_mlp - 2 * log_likelihood_mixture_laplace
    bic_mixture_laplace = k_mlp * np.log(n) - 2 * log_likelihood_mixture_laplace
    
    aic_bvnct = 2 * k_bvnct - 2 * log_likelihood_bvnct
    bic_bvnct = k_bvnct * np.log(n) - 2 * log_likelihood_bvnct

    return aic_mixture_laplace, bic_mixture_laplace, aic_bvnct, bic_bvnct

# Test of q3part2 with some simulated data
# This code part is the same as in question 1

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
samples = simulate_mixture_bivariate_laplace(pi, mu1, b1, Sigma1, mu2, b2, Sigma2, n)

# Initial guesses

initial_guess_bvlp = np.array([
    0.5,        # w1 (weight of the first component)
    0, 0,       # mu1_x, mu1_y (mean of the first component)
    np.log(10), np.log(10),  # b1_x, b1_y (scale of the first component)
    3, 3,       # mu2_x, mu2_y (mean of the second component)
    np.log(5), np.log(5)     # b2_x, b2_y (scale of the second component)
])

initial_guess_bvnct = np.array([
    0, 0,       # mu_x, mu_y (location vector)
    0.5, 0.5,   # gamma_x, gamma_y (noncentrality vector)
    10,         # v (degrees of freedom)
    1, 0, 1     # Sigma_xx, Sigma_xy, Sigma_yy (covariance matrix elements)
])

initial_guess = [
    initial_guess_bvlp,  
    initial_guess_bvnct   
]

# The order of printing is this one: aic_mixture_laplace, bic_mixture_laplace, aic_bvnct, bic_bvnct
print(aic_bic_bvlp_bvnct(samples, initial_guess))

#%% II.4
import pandas as pd
from scipy.stats import gaussian_kde

# Reading the data.
df = pd.read_csv("DJIA30stockreturns.csv")

def aic_bic_model_comparison_mlaplace_bvnct(data, rep, initial_guess):
    results_matrix_laplace_mixture = np.zeros((rep, 4))
    results_matrix_bvnct = np.zeros((rep, 4))
    for i in range(rep):
        random_stock1 = np.random.randint(0, 25)
        random_stock2 = np.random.randint(0, 25)
        stock_1 = data.iloc[:, random_stock1].to_numpy()
        stock_2 = data.iloc[:, random_stock2].to_numpy()
        aic_mixture_laplace1, bic_mixture_laplace1, aic_bvnct1, bic_bvnct1 = aic_bic_bvlp_bvnct(stock_1, initial_guess)
        aic_mixture_laplace2, bic_mixture_laplace2, aic_bvnct2, bic_bvnct2 = aic_bic_bvlp_bvnct(stock_2, initial_guess)
        results_matrix_laplace_mixture[i, :] = [aic_mixture_laplace1, bic_mixture_laplace1, aic_mixture_laplace2, bic_mixture_laplace2]
        results_matrix_bvnct[i, :] = [aic_bvnct1, bic_bvnct1, aic_bvnct2, bic_bvnct2]
    return results_matrix_laplace_mixture, results_matrix_bvnct



# Initial guesses

initial_guess_bvlp = np.array([
    0.5,        # w1 (weight of the first component)
    0, 0,       # mu1_x, mu1_y (mean of the first component)
    np.log(10), np.log(10),  # b1_x, b1_y (scale of the first component)
    3, 3,       # mu2_x, mu2_y (mean of the second component)
    np.log(5), np.log(5)     # b2_x, b2_y (scale of the second component)
])

initial_guess_bvnct = np.array([
    0, 0,       # mu_x, mu_y (location vector)
    0.5, 0.5,   # gamma_x, gamma_y (noncentrality vector)
    10,         # v (degrees of freedom)
    1, 0, 1     # Sigma_xx, Sigma_xy, Sigma_yy (covariance matrix elements)
])

initial_guess = [
    initial_guess_bvlp,  
    initial_guess_bvnct   
]

laplace_aic_bic, bvnct_aic_bic = aic_bic_model_comparison_mlaplace_bvnct(df, 50, initial_guess)


# This is not working because we need some way to give the data to the mixture laplace from one array to a dx2 sized matrix
# also there is a problem with the translation to python of the nct (again)


















