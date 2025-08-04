import numpy as np
from data_process import generate_truncated_lognormal
from rebound_model import calculate_rebound_2D, calculate_rebound_3D

def get_exp_data_array(coarse_d1, medium_d1, fine_d1):
    """Get experimental data array for rebound angle and restitution coefficient."""
    Beladjine07_v26= {
        'd_in': [0.006, 0.006, 0.006, 0.006, 0.006],
        'v_in': [26, 26, 26, 26, 26],
        'ang_in': [10, 20, 40, 60, 90],
        'ang_re': [21.21, 26.68, 33.87, 40.94, 90.01],
        'e': [0.77, 0.61, 0.43, 0.26, 0.22],
        'Nej': [9, 14, 17, 21, 19],
        'v_ej': [0.92, 0.86, 0.85, 0.9, 0.83]
    }
    Zhou06_v_many= {
        'ang_in': [8, 11.5],
        'ang_re_mean': [46.97, 47.53],
        'ang_re_std': [0.59, 0.75],
        'e_mean': [0.63, 0.61],
        'e_std': [0.006, 0.009]
    }
    Rice95_v_many_coarse= {
        'd_in': coarse_d1,
        'v_in': np.mean([2.7323, 2.7013, 2.7070, 2.7979]),
        'ang_in': np.mean([13.94, 14.75, 14.73, 15.04]),
        'ang_re': np.mean([20.95, 23.03, 22.55, 25.63]),
        'e': np.mean([0.58, 0.56, 0.56, 0.58]),
        'Nej': np.mean([5.7, 5.09, 5.12, 4.64]),
        'v_ej': np.mean([0.243, 0.2422, 0.2522, 0.2637])
    }
    Rice95_v_many_medium= {
        'd_in': medium_d1,
        'v_in': np.mean([3.1649, 3.3638, 3.3604, 3.2963]),
        'ang_in': np.mean([11.82, 11.47, 11.53, 11.36]),
        'ang_re': np.mean([29.95, 30.31, 31.06, 29.56]),
        'e': np.mean([0.57, 0.54, 0.59, 0.58]),
        'Nej': np.mean([2.68, 3.43, 2.67, 2.76]),
        'v_ej': np.mean([0.2544, 0.252, 0.2607, 0.295])
    }
    Rice95_v_many_fine= {
        'd_in': fine_d1,
        'v_in': np.mean([3.8493, 3.7801, 3.7411]),
        'ang_in': np.mean([10.85, 10.24, 10.46]),
        'ang_re': np.mean([44.63, 38.03, 37.85]),
        'e': np.mean([0.52, 0.56, 0.58]),
        'Nej': np.mean([1.86, 1.71, 1.58]),
        'v_ej': np.mean([0.2757, 0.2566, 0.2773])
    }
    Chen18_v_many= {
        'ang_in': [23.2, 21.8, 21.9, 30.7, 30.6, 37.6, 47.1, 46.6, 46],
        'ang_re': [32.83, 35.59, 29.90, 43, 33.28, 44.16, 57.62, 50.97, 39.88],
        'e': [0.38, 0.4, 0.45, 0.32, 0.37, 0.27, 0.2, 0.19, 0.22]
    }
    Willetts89_v_many_coarse= {
        'd_in': [coarse_d1, coarse_d1, coarse_d1, coarse_d1],
        'v_in': [3.38, 3.43, 3.18, 3.50],
        'ang_in': [12.7, 17.8, 23.2, 27.7],
        'ang_re': [19.1, 25.2, 21.4, 27.2],
        'e': [0.63, 0.57, 0.54, 0.46],
        'v_ej': [0.3718, 0.3773, 0.4134, 0.42]
    }
    Willetts89_v_many_medium= {
        'd_in': [medium_d1, medium_d1, medium_d1, medium_d1],
        'v_in': [3.56, 3.99, 4.02, 4.39],
        'ang_in': [11.7, 18.2, 21.4, 26.3],
        'ang_re': [24.9, 33.4, 33.3, 44.7],
        'e': [0.61, 0.53, 0.50, 0.40],
        'v_ej': [0.356, 0.399, 0.3618, 0.439]
    }
    Willetts89_v_many_fine= {
        'd_in': [fine_d1, fine_d1, fine_d1, fine_d1],
        'v_in': [3.61, 4.41, 4.26, 4.35],
        'ang_in': [9.5, 15.4, 19.7, 24.9],
        'ang_re': [38.8, 42, 42.2, 42.5],
        'e': [0.57, 0.50, 0.48, 0.46],
        'v_ej': [0.3249, 0.3969, 0.3834, 0.348]
    }
    Rioual20_v_many= {
        'ang_in': 53.0,
        'e': 0.37813
    }
    Gordon21_v_many= {
        'ang_in': [16.7, 11.0],
        'e': [0.6, 0.62]
    }
    Gordon09_v_many= {
        'ang_in': [8.5, 9.5],
        'ang_re': [22.8, 18.0],
        'e': [0.79, 0.64]
    }

    # Experimental data of impactor
    d1_array = [Rice95_v_many_coarse['d_in'], Rice95_v_many_medium['d_in'], Rice95_v_many_fine['d_in']]
    d1_array = d1_array + Willetts89_v_many_coarse['d_in'] + Willetts89_v_many_medium['d_in'] + Willetts89_v_many_fine['d_in']
    d1_array = np.array(d1_array)
    theta1_array = [Rice95_v_many_coarse['ang_in'], Rice95_v_many_medium['ang_in'], Rice95_v_many_fine['ang_in']]
    theta1_array = theta1_array + Willetts89_v_many_coarse['ang_in'] + Willetts89_v_many_medium['ang_in'] + Willetts89_v_many_fine['ang_in']
    theta1_array = np.array(theta1_array)
    theta1_array = np.deg2rad(theta1_array)  # Convert angles to radians
    # Experimental data of rebound
    E_array = [Rice95_v_many_coarse['e'], Rice95_v_many_medium['e'], Rice95_v_many_fine['e']]
    E_array = E_array + Willetts89_v_many_coarse['e'] + Willetts89_v_many_medium['e'] + Willetts89_v_many_fine['e']
    E_array = np.array(E_array)
    Phi_array = [Rice95_v_many_coarse['ang_re'], Rice95_v_many_medium['ang_re'], Rice95_v_many_fine['ang_re']]
    Phi_array = Phi_array + Willetts89_v_many_coarse['ang_re'] + Willetts89_v_many_medium['ang_re'] + Willetts89_v_many_fine['ang_re']
    Phi_array = np.array(Phi_array)
    Phi_array = np.deg2rad(Phi_array)  # Convert angles to radians

    return d1_array, theta1_array, Phi_array, E_array

def sample_averaged_rebound(epsilon, nu, d1, th, dist_params, bed_type):
    """Sample average rebound angle and restitution coefficient (average over bed PSD)."""
    phi_list = []
    e_list = []
    # sampling
    mu = dist_params['mu']
    sigma = dist_params['sigma']
    d_min = dist_params['d_min']
    d_max = dist_params['d_max']
    sampling_num = dist_params['sampling_num']
    d2_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, sampling_num)
    d3_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, sampling_num)
    d4_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, sampling_num)
    # calculating
    for i in range(sampling_num):
        if bed_type['three_D']:
            phi, e = calculate_rebound_3D(d1, d2_array[i], d3_array[i], d4_array[i], th, epsilon, nu)
        else:
            phi, e = calculate_rebound_2D(d1, d2_array[i], d3_array[i], th, epsilon, nu)
        phi_list.append(phi)
        e_list.append(e)
    phi_avg = np.mean(phi_list)
    e_avg = np.mean(e_list)

    return phi_avg, e_avg

def get_model_data_array(epsilon, nu, d1_array, th_array, dist_params, bed_type):
    """Get model data array for rebound angle and restitution coefficient."""
    phi_list = []
    e_list = []
    for d1, th in zip(d1_array, th_array):
        if bed_type['monodisperse']:
            d2 = bed_type['d50']
            d3 = d2
            if bed_type['three_D']:
                d4 = d2
                phi, e = calculate_rebound_3D(d1, d2, d3, d4, th, epsilon, nu)
            else:
                phi, e = calculate_rebound_2D(d1, d2, d3, th, epsilon, nu)
            phi_list.append(phi)
            e_list.append(e)
        else:
            phi, e = sample_averaged_rebound(epsilon, nu, d1, th, dist_params, bed_type)
            phi_list.append(phi)
            e_list.append(e)

    return np.array(phi_list), np.array(e_list)