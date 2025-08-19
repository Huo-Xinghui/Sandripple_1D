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

def normalize_data(Phi_array, E_array, phi_array, e_array):
    """Normalize data to zero mean and unit variance."""
    Phi_mean = np.mean(Phi_array)
    E_mean = np.mean(E_array)
    Phi_std = np.std(Phi_array, ddof=1)
    E_std = np.std(E_array, ddof=1)
    Phi_norm = (Phi_array - Phi_mean) / Phi_std
    E_norm = (E_array - E_mean) / E_std
    phi_norm = (phi_array - Phi_mean) / Phi_std
    e_norm = (e_array - E_mean) / E_std
    return Phi_norm, E_norm, phi_norm, e_norm

def print_time(ct, lt, st, cs, ts):
    """Printing remaining minutes and seconds from a time in seconds."""
    print_dict={
        1: 'thd=15, simple model data vs sigma calculated',
        2: 'thd=30, simple model data vs sigma calculated',
        3: 'thd=45, simple model data vs sigma calculated',
        4: 'thd=15, 3D model data vs sigma calculated',
        5: 'thd=30, 3D model data vs sigma calculated',
        6: 'thd=45, 3D model data vs sigma calculated',
        7: 'thd=15, monodisperse model data vs sigma calculated',
        8: 'thd=30, monodisperse model data vs sigma calculated',
        9: 'thd=45, monodisperse model data vs sigma calculated',
    }
    cost_time = ct - lt
    total_time = ct - st
    min = cost_time // 60
    sec = cost_time % 60
    remaining_time = total_time/cs * (ts - cs)
    rem_min = remaining_time // 60
    rem_sec = remaining_time % 60
    tot_min = total_time // 60
    tot_sec = total_time % 60
    print(print_dict[cs])
    print(f'{cs}/{ts} Cost: {min:.0f} min {sec:.0f} s. Remaining: {rem_min:.0f} min {rem_sec:.0f} s. Total: {tot_min:.0f} min {tot_sec:.0f} s.')