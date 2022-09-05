#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 19:36:14 2022

@author: alema80
"""

import math
import numpy as np
from matplotlib import pyplot as plt


def dcs_R(theta, beta, z, m0c2, Z):
  # Relativistic unscreened Rutherford differential cross section
  #
  # theta  scattering angle in rad
  # beta   velocity of the particle relative to the speed of light c
  # z      charge of the particle in e
  # m0c2   rest energy of the particle
  # Z      atomic number of the medium
  
  r_e = 2.818e-13  # Classical electron radius in cm
  mec2 = 0.511     # Rest energy of electron in MeV
  return( (r_e * z * Z)**2 *(mec2 / m0c2)**2 * (1 - beta**2) / beta**2 * 1 / (1 - math.cos(theta)**2))


def dcs_R_screened(theta, beta, z, m0c2, Z):
  # Relativistic screened Rutherford differential cross section
  #
  # theta  scattering angle
  # beta   velocity of the particle relative to the speed of light c
  # z      charge of the particle in e
  # m0c2   rest energy of the particle
  # Z      atomic number of the medium
  
  r_e = 2.818e-13  # Classical electron radius in cm
  mec2 = 0.511     # Rest energy of electron in MeV
  return( (r_e * z * Z)**2 *(mec2 / m0c2)**2 * (1 - beta**2) / beta**2 * 1 / (1 - math.cos(theta) + 0.5 * chi_a**2)**2)


def dcs_M(theta, beta, z, m0c2, Z):
  # Relativistic unscreened Mott differential cross section
  #
  # theta  scattering angle
  # beta   velocity of the particle relative to the speed of light c
  # z      charge of the particle in e
  # m0c2   rest energy of the particle in MeV
  # Z      atomic number of the medium

  alpha = 0.0072973525693  # Fine structure constant (approx. 1/137)
  dcs_R_u = dcs_R(theta, beta, z, m0c2, Z)
  return( dcs_R_u * (1 - beta**2 * (math.sin(theta/2))**2 + math.pi * alpha * beta * Z * math.sin(theta/2) * (1 - math.sin(theta/2) ) ))


# Define the incident particle
beta = 0.95   # velocity relative to c
z = 1         # electron
m0c2 = 0.511  # rest energy in MeV

# Define the material
Z = 6        # carbon

# Calculate and print screening angles chi_0 and chi_a
mec2 = 0.511  # Rest energy of electron in MeV
chi_0 = 4.2121e-3 * math.sqrt(1 - beta**2) / (mec2 * beta) * z**(1/3)
chi_a = chi_0 * math.sqrt(1.13 + 3.76 * ( (z * Z) / (137 * beta) )**2)
print('Chi_0 = {:f} rad'.format(chi_0))
print('Chi_a = {:f} rad'.format(chi_a))

# Compute plot points
nPoints = 128
theta_a = np.linspace(1e-5, (math.pi - 0.001), nPoints)
dcs_R_u_a = np.empty(nPoints)
dcs_R_s_a = np.empty(nPoints)
dcs_M_u_a = np.empty(nPoints)
for i in range(nPoints):
     dcs_R_u_a[i] = dcs_R(theta_a[i], beta, z, m0c2, Z)
     dcs_R_s_a[i] = dcs_R_screened(theta_a[i], beta, z, m0c2, Z)
     dcs_M_u_a[i] = dcs_M(theta_a[i], beta, z, m0c2, Z)
     

# Plot the figure
plt.plot(theta_a, dcs_R_u_a, 'r', label = 'Rutherford_unscreened')
plt.plot(theta_a, dcs_R_s_a, 'g', label = 'Rutherford_screened')
plt.plot(theta_a, dcs_M_u_a, 'b', label = 'Mott_unscreened')
plt.semilogy()
plt.xlabel(r'$\theta$ (rad)')
plt.ylabel(r'$\frac{d \sigma_R}{d \Omega}$  (cm$^2$ / rad)')
plt.legend(loc = 'upper center',)
plt.show()