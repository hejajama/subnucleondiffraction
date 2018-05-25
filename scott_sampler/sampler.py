#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

# Maximum a posteriori parameters
radius = 0.87 # [fm]
structure_parameter = 0.3
constituent_number = 5
sigma_fluct = 0.81

# Gaussian width of each constituent
constituent_width = 0.2 + structure_parameter*(radius - 0.2)

# Gaussian sampling radius for the center of each constituent
sampling_radius = np.sqrt(radius**2 - constituent_width**2)

# Each constituent is weighted by a Gamma random variable with unit mean.
# The variance is controlled by the Gamma shape parameter k.
gamma_shape = 1./(constituent_number*sigma_fluct**2)

# Sample several example protons with constituent substructure
samples = 9
positions = np.random.normal(
    scale=sampling_radius,
    size=2*constituent_number*samples
).reshape(constituent_number, 2, -1)

# Create 2D grid
l = np.linspace(-3, 3, 100)
xx, yy = np.meshgrid(l, l)

# Plot several example nucleon densities
fig, axes = plt.subplots(
    nrows=3, ncols=3,
    sharex=True, sharey=True
)

# Generate figure
for ax, pos in zip(axes.flat, positions.T):
    rho = np.zeros_like(xx)

    for (x0, y0) in pos.T:
        weight = np.random.gamma(shape=gamma_shape, scale=1./gamma_shape)
        norm = weight/(2*np.pi*constituent_width**2)
        print x0,y0,norm
        rho += norm*np.exp(-((xx-x0)**2 + (yy-y0)**2)/(2*constituent_width**2))

    ax.imshow(rho, extent=(-3, 3, -3, 3))

    if ax.is_last_row():
        ax.set_xlabel('x [fm]')
    if ax.is_first_col():
        ax.set_ylabel('y [fm]')

plt.tight_layout(w_pad=0, h_pad=0)
plt.savefig('proton_substructure.pdf')
