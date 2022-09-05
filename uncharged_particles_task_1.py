#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 16:22:50 2022

@author: alema80
"""

import math

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


# Incident particle
E = 1.0   # energy of the photon in MeV

# Material
Z = 29
A = 63.546

# Calculate the Klein-Nishina cross section
sigma_KN_v =  sigma_KN(E)

# Calculate the Compton mass attenuation coefficient for copper according to Klein-Nishina
mu_m_KN_Cu = mu_m_KN(E, Z, A)

# Print the results
print('Klein-Nishina cross section = {:g} cm^2'.format(sigma_KN_v))
print('Compton mass attenuation coefficient = {:g} cm^2 g^-1'.format(mu_m_KN_Cu))