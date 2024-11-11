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


