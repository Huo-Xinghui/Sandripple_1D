import numpy as np
import sys
from scipy.optimize import minimize
from data_process import generate_truncated_lognormal, prepare_data_for_CCA
from get_data import get_exp_data_array, get_model_data_array
from cca_calculation import cca_objective

# Optimize parameters
opt_method = 'Nelder-Mead'
parameter_tolerance = 1e-6
function_tolerance = 1e-6
max_iterations = 500
disp_process = True # display the optimization process
# Boundaries of epsilon and nu (invalid for Nelder-Mead)
epsilon_min = 0.1
epsilon_max = 0.9
nu_min = -1.0
nu_max = 0.0
# Initial guess for optimization
initial_guess = [0.7, -0.8]  # [epsilon, nu]
# Bed PSD parameters
dist_params = {
    'd_min': 1.5e-4,
    'd_max': 6e-4,
    'mu': -8.30271,
    'sigma': 0.25778,
    'sampling_num': 100
}
# Average bed diameter
sampling_num = 10000
mu = dist_params['mu']
sigma = dist_params['sigma']
d_min = dist_params['d_min']
d_max = dist_params['d_max']
d2_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, sampling_num)
d2_mid = np.percentile(d2_array, 50)
# Impactor diameters
d1_dict = {
    'coarse': 1.4*d2_mid,
    'medium': d2_mid,
    'fine': 0.73*d2_mid
}
# Control parameters
test_normalization = True # test the normalization of data
bed_type = {
    'three_D': False,  # 3D bed
    'monodisperse': False,  # monodisperse bed
    'd50': d2_mid,  # average bed diameter
}

# Get experimental data
d1_array, th_array, Phi_array, E_array = get_exp_data_array(d1_dict['coarse'], d1_dict['medium'], d1_dict['fine'])

if test_normalization:
    # Test the normalization of data
    test_epsilon = 0.7
    test_nu = -0.8
    phi_array, e_array = get_model_data_array(test_epsilon, test_nu, d1_array, th_array, dist_params, bed_type)
    exp_data_n, model_data_n = prepare_data_for_CCA(Phi_array, E_array, phi_array, e_array) 
    print("Characteristic of the normalized data:")
    print(f"Phi_norm: Mean = {np.mean(exp_data_n[:, 0]):.4f}, Std = {np.std(exp_data_n[:, 0], ddof=1):.4f}")
    print(f"E_norm: Mean = {np.mean(exp_data_n[:, 1]):.4f}, Std = {np.std(exp_data_n[:, 1], ddof=1):.4f}")
    cov_matrix = np.cov(exp_data_n, rowvar=False)
    print("Covariance matrix of the normalized data:")
    print(cov_matrix)
    sys.exit(0)

bounds = [(epsilon_min, epsilon_max), (nu_min, nu_max)]
result = minimize(
    fun =cca_objective,
    x0 = initial_guess,
    args = (d1_array, th_array, Phi_array, E_array, dist_params, bed_type),
    method = opt_method,
    bounds = bounds,
    options = {
        'maxiter': max_iterations,
        'disp': disp_process,
        'xatol': parameter_tolerance,
        'fatol': function_tolerance
    }
)
epsilon_opt, nu_opt = result.x
min_error = result.fun
rho_1 = 1 - min_error

print(f"\n Optimization results:")
print(f"epsilon = {epsilon_opt:.4f}")
print(f"nu = {nu_opt:.4f}")
print(f"Maximal canonical correlation coefficient (rho_1) = {rho_1:.4f}")
print(f"explained variance ratio = {rho_1**2:.4f} ({rho_1**2*100:.2f}%)")