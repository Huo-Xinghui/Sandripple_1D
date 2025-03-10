"""Based on the paper
'Unified Model of Sediment Transport Threshold and Rate Across Weak and Intense Subaqueous Bedload, Windblown Sand, and Windblown Snow',
Pahtz et al. 2021, JGR"""

import numpy as np
from scipy.special import lambertw

def calculate_vz_down(vz_e):
    """Calculate the settling velocity of a particle in a fluid"""
    val_in_lambertw = -(1 + vz_e)*np.exp(-1 - vz_e)
    w = lambertw(val_in_lambertw)
    vz_s = -1 - w
    return vz_s

def calculate_e(d1, d2, epsilon, nu, vz_s):
    """Calculate the restitution coefficients"""
    mu = ((d2 / d1)**3 + 1 / epsilon)**(-1)
    alpha = (1 + epsilon) / (1 + mu) - 1
    beta = 1 - (2 / 7) * (1 - nu) / (1 + mu)
    return e, e_z