import numpy as np
from scipy.stats import truncnorm

def generate_truncated_lognormal(mu, sigma, dmin, dmax, size):
    """Generate truncated log-normal samples."""
    a = (np.log(dmin) - mu) / sigma
    b = (np.log(dmax) - mu) / sigma
    samples = truncnorm.rvs(a, b, loc=mu, scale=sigma, size=size)
    return np.exp(samples)

def prepare_data_for_CCA(Phi_array, E_array, phi_array, e_array):
    """ Prepare normalized data for CCA analysis."""
    # Mean and standard deviation
    Phi_mean = np.mean(Phi_array)
    E_mean = np.mean(E_array)
    Phi_std = np.std(Phi_array, ddof=1)
    E_std = np.std(E_array, ddof=1)
    # Normalize the data
    exp_data = np.column_stack([
        (Phi_array - Phi_mean) / Phi_std,
        (E_array - E_mean) / E_std
    ])
    model_data = np.column_stack([
        (phi_array - Phi_mean) / Phi_std,
        (e_array - E_mean) / E_std
    ])
    return exp_data, model_data