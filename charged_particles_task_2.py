#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 19:36:14 2022

@author: alema80
"""

import math


def mSel_e(beta, Z, A, I, delta):
    # Mass electronic stopping power for electrons
    #
    # beta   velocity of the electron relative to c
    # Z      atomic number
    # A      mass number
    # I      mean exitation potential in eV
    # delta  density effect correction
    
    c1 = 0.153537  # MeV cm^2
    mec2 = 0.511  # Electron rest energy in MeV
    tau = 1 / math.sqrt(1 - beta**2) - 1
    E = tau * mec2
    return(c1 * Z / A * 1 / beta**2 * (math.log((E/(1e-6 * I))**2) + math.log(1 + tau / 2) + F_el(tau) - delta)) 


def mSel_p(beta, Z, A, I, delta):
    # Mass electronic stopping power for positrons
    # beta   velocity of the electron relative to c
    # Z      atomic number
    # A      mass number
    # I      mean exitation potential in eV
    # delta  density effect correction
    
    c1 = 0.153537  # MeV cm^2
    mec2 = 0.511  # Electron rest energy in MeV
    tau = 1 / math.sqrt(1 - beta**2) - 1
    E = tau * mec2
    return(c1 * Z / A * 1 / beta**2 * (math.log((E / (1e-6 * I))**2) + math.log(1 + tau / 2) + F_po(tau) - delta)) 


def F_el(tau):
    # The F-function for stopping power of electrons
    #
    # tau  kinetic energy relative to rest energy
    
    beta = math.sqrt(1 - 1 / (tau + 1)**2)
    return((1 - beta**2) * (1 + tau**2 / 8 - (2 * tau + 1) * math.log(2)))


def F_po(tau):
    # The F-function for stopping power of positrons
    #
    # tau  kinetic energy relative to rest energy
    
    beta = math.sqrt(1 - 1 / (tau + 1)**2)
    return(2 * math.log(2) - beta**2 / 12 * (23 + 14 / (tau + 2) + 10 / (tau + 2)**2 + 4 / (tau + 2)**3))


def delta(beta):
    # The density effect correction using Sternhaimer's model
    #
    # beta   velocity of the particle relative to c
    
    # The following numbers are specific for Al
    a = 0.08024
    X0 = 0.1708
    X1 = 3.0127
    C = -4.2395
    m = 3.6345
    
    tau = 1 / math.sqrt(1 - beta**2) - 1
    X = 0.5 * math.log10(tau * (tau + 2))
    if X > X1:
        return(4.6052 * X + C)
    if X > X0:
        return(4.6052 * X + a * (X1 - X)**m + C)
    return(0.0)


# Define the incident particle    
E = 50.0  # kinetic energy of the particle in MeV

# Define the material
Z = 13     # Atomic number
A = 26.98  # Mass number
I = 166    # Mean excitation energy in eV

# Compute beta corresponding to the kinetic energy E
mec2 = 0.511   # Electron rest energy in MeV
tau = 50.0 / mec2 
beta = math.sqrt(1 - 1 / (tau + 1)**2)

# Compute electronic stopping powers
mSel_electron_nodelta = mSel_e(beta, Z, A, I, 0.0)
mSel_positron_nodelta = mSel_p(beta, Z, A, I, 0.0)

deltaX = delta(beta) 
mSel_electron = mSel_e(beta, Z, A, I, deltaX)
mSel_positron = mSel_p(beta, Z, A, I, deltaX)

# Print the results
print('mass electronic stopping power for electron (no density effect correction) = {:g} MeV cm^2 g^-1'.format(mSel_electron_nodelta))
print('mass electronic stopping power for positron (no density effect correction) = {:g} MeV cm^2 g^-1'.format(mSel_positron_nodelta))

print('mass electronic stopping power for electron = {:g} MeV cm^2 g^-1'.format(mSel_electron))
print('mass electronic stopping power for positron = {:g} MeV cm^2 g^-1'.format(mSel_positron))
