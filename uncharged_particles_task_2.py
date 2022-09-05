#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 16:22:50 2022

@author: alema80
"""

import math
import numpy as np
from matplotlib import pyplot as plt


def sigma_KN(E):
    # Total Klein-Nishina cross section
    #
    # E   energy of the photon in MeV
    
    re = 2.818e-13  # classical electron radius in cm
    mec2 = 0.511    # rest energy of electron in MeV
    eps = E / mec2
    return(2 * math.pi * re**2 * ((1 + eps) / eps**2 * (2 * (1 + eps) / (1 + 2 * eps) - math.log(1 + 2 * eps) / eps) + math.log(1 + 2 * eps) / (2 * eps) - (1 + 3 * eps) / (1 + 2 * eps)**2))

# Total Thomson cross section
re = 2.818e-13  # classical electron radius in cm
sigma_Th = 8/3 * math.pi * re**2


# Compute plot points
nPoints = 128
energy_a = np.logspace(-3, 2, nPoints)
sigma_KN_a = np.empty(nPoints)
sigma_Th_a = np.empty(nPoints)
for i in range(nPoints):
    sigma_KN_a[i] = sigma_KN(energy_a[i])
    sigma_Th_a[i] = sigma_Th
     
# Plot the figure
plt.plot(energy_a, sigma_KN_a, 'r', label = 'Klein-Nishina')
plt.plot(energy_a, sigma_Th_a, 'b', label = 'Thomson')
plt.semilogx()
plt.semilogy()
plt.xlabel(r'energy (MeV)')
plt.ylabel(r'$\sigma_{KN}$  (cm$^2$)')
plt.legend(loc = 'lower left',)
plt.show()