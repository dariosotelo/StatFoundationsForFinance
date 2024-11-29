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
from scipy.special import gammaln
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
    
graph_difference(nu, x, x_values)



































