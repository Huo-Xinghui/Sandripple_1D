import numpy as np
from sklearn.cross_decomposition import CCA
from data_process import prepare_data_for_CCA
from get_data import get_model_data_array

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