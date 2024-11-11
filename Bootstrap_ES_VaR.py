#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 18:58:34 2024

"""

#%% Libraries
from scipy.special import betainc

#%% Preparation functions

def f_GAt(z, d, nu, theta, K):
    if any(param <= 0 for param in [d, nu, theta]):
        return "d, nu, or theta must be positive."
    if z < 0:
        return K * (1 + (-z * theta) ** d / nu) ** -(nu + 1/d)
    else:
        return K * (1 + (z / theta) ** d / nu) ** -(nu + 1/d)

def calculate_L(c, nu, theta, d):
    return nu / (nu + (-c * theta) ** d)

def S_r(c, r, d, nu, theta):
    if c >= 0:
        return "c must be negative."
    L = calculate_L(c, nu, theta, d)
    B_L_num = betainc(nu - r / d, (r + 1) / d, L)  
    B_L_den = betainc(nu, 1 / d, L)                
    return ((-1) ** r) * (nu ** (r / d)) * ((1 + theta ** 2) / (theta ** r + theta ** (r + 2))) * (B_L_num / B_L_den)

#The ES would be calculated using VaR as c and r=1



