import numpy as np
from sklearn.cross_decomposition import CCA
from data_process import prepare_data_for_CCA, normalize_data
from get_data import get_model_data_array
from sklearn.covariance import LedoitWolf, OAS

def cca_objective(params, d1_array, th_array, Phi_array, E_array, dist_params, bed_type):
    """Objective function for CCA: minimize 1-rho_1 (maximal canonical correlation coefficient)."""
    epsilon, nu = params
    phi_array, e_array = get_model_data_array(epsilon, nu, d1_array, th_array, dist_params, bed_type)
    exp_data, model_data = prepare_data_for_CCA(Phi_array, E_array, phi_array, e_array)
    # Check if there are NaN values in the data
    if np.isnan(exp_data).any() or np.isnan(model_data).any():
        return np.inf
    # Perform CCA
    try:
        cca = CCA(n_components=1)
        cca.fit(exp_data, model_data)
        exp_c, model_c = cca.transform(exp_data, model_data)
        rho_1 = np.corrcoef(exp_c[:, 0], model_c[:, 0])[0, 1]
        return 1 - rho_1  # We want to maximize rho_1, so we minimize 1 - rho_1
    except Exception as e:
        print(f"CCA failed: {e}")
        return np.inf

def corr_objective(params, d1_array, th_array, Phi_array, E_array, dist_params, bed_type):
    """Objective function for correlation: minimize 1-rho_1 (maximal correlation coefficient)."""
    epsilon, nu = params
    phi_array, e_array = get_model_data_array(epsilon, nu, d1_array, th_array, dist_params, bed_type)
    exp_data, model_data = prepare_data_for_CCA(Phi_array, E_array, phi_array, e_array)
    # Check if there are NaN values in the data
    if np.isnan(exp_data).any() or np.isnan(model_data).any():
        return np.inf
    # Calculate correlation coefficient of the rebound angle
    corr_phi = np.corrcoef(exp_data[:, 0], model_data[:, 0])[0, 1]
    # Calculate correlation coefficient of the restitution coefficient
    corr_e = np.corrcoef(exp_data[:, 1], model_data[:, 1])[0, 1]
    # Calculate the overall correlation coefficient
    if np.isnan(corr_phi) or np.isnan(corr_e):
        return np.inf
    corr = np.sqrt(corr_phi*corr_e)
    if np.isnan(corr):
        return np.inf
    return 1 - corr_e

def Eucl_distance_objective(params, d1_array, th_array, Phi_array, E_array, dist_params, bed_type):
    """Objective function for Euclidean distance: minimize the Euclidean distance between model and experimental data."""
    epsilon, nu = params
    phi_array, e_array = get_model_data_array(epsilon, nu, d1_array, th_array, dist_params, bed_type)
    Phi_norm, E_norm, phi_norm, e_norm = normalize_data(Phi_array, E_array, phi_array, e_array)
    w_phi = 0.5  # Weight for rebound angle
    w_e = 0.5    # Weight for restitution coefficient
    distance = np.sqrt(
        w_phi * np.mean((Phi_norm - phi_norm) ** 2) +
        w_e * np.mean((E_norm - e_norm) ** 2)
    )
    return distance

def Maha_distance_objective(params, d1_array, th_array, Phi_array, E_array, dist_params, bed_type):
    """Objective function for Mahalanobis distance: minimize the Mahalanobis distance between model and experimental data."""
    epsilon, nu = params
    phi_array, e_array = get_model_data_array(epsilon, nu, d1_array, th_array, dist_params, bed_type)
    lw = LedoitWolf()
    lw.fit(np.column_stack((Phi_array, E_array)))
    S_lw = lw.covariance_
    #oas = OAS()
    #oas.fit(np.column_stack((Phi_array, E_array)))
    #S_oas = oas.covariance_
    diff = np.column_stack((phi_array - Phi_array, e_array - E_array))
    inv_S = np.linalg.inv(S_lw)
    return np.mean([d @ inv_S @ d for d in diff])  # Mahalanobis distance