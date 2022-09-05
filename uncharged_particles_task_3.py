#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 16:22:50 2022

@author: alema80
"""

import math
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

def sigma_KN(E):
    # Total Klein-Nishina cross section
    #
    # E   energy of the photon in MeV
    
    re = 2.818e-13  # classical electron radius in cm
    mec2 = 0.511    # rest energy of electron in MeV
    eps = E / mec2
    return(2 * math.pi * re**2 * ((1 + eps) / eps**2 * (2 * (1 + eps) / (1 + 2 * eps) - math.log(1 + 2 * eps) / eps) + math.log(1 + 2 * eps) / (2 * eps) - (1 + 3 * eps) / (1 + 2 * eps)**2))

def mu_m_KN(E, Z, A):
    # Compton mass attenuation coefficient according to Klein-Nishina
    #
    # E  # energy of the photon in MeV
    # Z  # atomic number of the material
    # A  # mass number of the material
    
    N_A = 6.022e23  # Avogadro number
    return(Z * N_A / A * sigma_KN(E))


# Material
Z = 6
A = 12

# Read cross section data from a file
url = 'https://raw.githubusercontent.com/AlexandrMalusek/MSF_course/main/xcom_carbon.csv'
df = pd.read_csv(url)
#df = pd.read_csv('xcom_carbon.csv')


# Compute plot points
nPoints = 128
energy_a = np.logspace(-3, 2, nPoints)
mu_m_KN_a = np.empty(nPoints)
for i in range(nPoints):
    mu_m_KN_a[i] = mu_m_KN(energy_a[i], Z, A)
     
# Plot the figure
plt.plot(df['energy'], df['mac_co'], 'r', label = 'Coherent scattering')
plt.plot(df['energy'], df['mac_in'], 'b', label = 'Inoherent scattering')
plt.plot(df['energy'], df['mac_co'] + df['mac_in'], 'g', label = 'Coh. + Incoh. scattering')
plt.plot(energy_a, mu_m_KN_a, 'y', label = 'Klein-Nishina')
plt.xlim([1e-3, 1])
plt.ylim(1e-3, 10)
plt.semilogx()
plt.semilogy()
plt.xlabel(r'energy (MeV)')
plt.ylabel(r'MAC  (cm$^2$/g)')
plt.legend(loc = 'upper right',)
plt.show()