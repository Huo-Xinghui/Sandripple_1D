import numpy as np
from data_process import generate_truncated_lognormal
from rebound_model import calculate_rebound_2D, calculate_rebound_3D

def get_exp_data_array(coarse_d1, medium_d1, fine_d1):
    """Get experimental data array for rebound angle and restitution coefficient."""
    Beladjine07= {
        'd1': [0.006, 0.006, 0.006, 0.006, 0.006],
        'v1': [26, 26, 26, 26, 26],
        'thd': [10, 20, 40, 60, 90],
        'phid': [21.21, 26.68, 33.87, 40.94, 90.01],
        'e': [0.77, 0.61, 0.43, 0.26, 0.22],
        'Nej': [9, 14, 17, 21, 19],
        'vn': [0.92, 0.86, 0.85, 0.9, 0.83]
    }
    Zhou06= {
        'thd': [8, 11.5],
        'phid': [46.97, 47.53],
        'phid_std': [0.59, 0.75],
        'e': [0.63, 0.61],
        'e_std': [0.006, 0.009]
    }
    Rice95_coarse= {
        'd1': [coarse_d1],
        'v1': [np.mean([2.7323, 2.7013, 2.7070, 2.7979])],
        'v1_std': [np.std([2.7323, 2.7013, 2.7070, 2.7979], ddof=1)],
        'thd': [np.mean([13.94, 14.75, 14.73, 15.04])],
        'thd_std': [np.std([13.94, 14.75, 14.73, 15.04], ddof=1)],
        'phid': [np.mean([20.95, 23.03, 22.55, 25.63])],
        'phid_std': [np.std([20.95, 23.03, 22.55, 25.63], ddof=1)],
        'e': [np.mean([0.58, 0.56, 0.56, 0.58])],
        'e_std': [np.std([0.58, 0.56, 0.56, 0.58], ddof=1)],
        'Nej': [np.mean([5.7, 5.09, 5.12, 4.64])],
        'Nej_std': [np.std([5.7, 5.09, 5.12, 4.64], ddof=1)],
        'vn': [np.mean([0.243, 0.2422, 0.2522, 0.2637])],
        'vn_std': [np.std([0.243, 0.2422, 0.2522, 0.2637], ddof=1)]
    }
    Rice95_medium= {
        'd1': [medium_d1],
        'v1': [np.mean([3.1649, 3.3638, 3.3604, 3.2963])],
        'v1_std': [np.std([3.1649, 3.3638, 3.3604, 3.2963], ddof=1)],
        'thd': [np.mean([11.82, 11.47, 11.53, 11.36])],
        'thd_std': [np.std([11.82, 11.47, 11.53, 11.36], ddof=1)],
        'phid': [np.mean([29.95, 30.31, 31.06, 29.56])],
        'phid_std': [np.std([29.95, 30.31, 31.06, 29.56], ddof=1)],
        'e': [np.mean([0.57, 0.54, 0.59, 0.58])],
        'e_std': [np.std([0.57, 0.54, 0.59, 0.58], ddof=1)],
        'Nej': [np.mean([2.68, 3.43, 2.67, 2.76])],
        'Nej_std': [np.std([2.68, 3.43, 2.67, 2.76], ddof=1)],
        'vn': [np.mean([0.2544, 0.252, 0.2607, 0.295])],
        'vn_std': [np.std([0.2544, 0.252, 0.2607, 0.295], ddof=1)]
    }
    Rice95_fine= {
        'd1': [fine_d1],
        'v1': [np.mean([3.8493, 3.7801, 3.7411])],
        'v1_std': [np.std([3.8493, 3.7801, 3.7411], ddof=1)],
        'thd': [np.mean([10.85, 10.24, 10.46])],
        'thd_std': [np.std([10.85, 10.24, 10.46], ddof=1)],
        'phid': [np.mean([44.63, 38.03, 37.85])],
        'phid_std': [np.std([44.63, 38.03, 37.85], ddof=1)],
        'e': [np.mean([0.52, 0.56, 0.58])],
        'e_std': [np.std([0.52, 0.56, 0.58], ddof=1)],
        'Nej': [np.mean([1.86, 1.71, 1.58])],
        'Nej_std': [np.std([1.86, 1.71, 1.58], ddof=1)],
        'vn': [np.mean([0.2757, 0.2566, 0.2773])],
        'vn_std': [np.std([0.2757, 0.2566, 0.2773], ddof=1)]
    }
    Chen18= {
        'thd': [23.2, 21.8, 21.9, 30.7, 30.6, 37.6, 47.1, 46.6, 46],
        'phid': [32.83, 35.59, 29.90, 43, 33.28, 44.16, 57.62, 50.97, 39.88],
        'e': [0.38, 0.4, 0.45, 0.32, 0.37, 0.27, 0.2, 0.19, 0.22]
    }
    Willetts89_coarse= {
        'd1': [coarse_d1, coarse_d1, coarse_d1, coarse_d1],
        'v1': [3.38, 3.43, 3.18, 3.50],
        'thd': [12.7, 17.8, 23.2, 27.7],
        'phid': [19.1, 25.2, 21.4, 27.2],
        'e': [0.63, 0.57, 0.54, 0.46],
        'vn': [0.3718, 0.3773, 0.4134, 0.42]
    }
    Willetts89_medium= {
        'd1': [medium_d1, medium_d1, medium_d1, medium_d1],
        'v1': [3.56, 3.99, 4.02, 4.39],
        'thd': [11.7, 18.2, 21.4, 26.3],
        'phid': [24.9, 33.4, 33.3, 44.7],
        'e': [0.61, 0.53, 0.50, 0.40],
        'vn': [0.356, 0.399, 0.3618, 0.439]
    }
    Willetts89_fine= {
        'd1': [fine_d1, fine_d1, fine_d1, fine_d1],
        'v1': [3.61, 4.41, 4.26, 4.35],
        'thd': [9.5, 15.4, 19.7, 24.9],
        'phid': [38.8, 42, 42.2, 42.5],
        'e': [0.57, 0.50, 0.48, 0.46],
        'vn': [0.3249, 0.3969, 0.3834, 0.348]
    }
    Rioual20= {
        'thd': 53.0,
        'e': 0.37813
    }
    Gordon21= {
        'thd': [16.7, 11.0],
        'e': [0.6, 0.62]
    }
    Gordon09= {
        'thd': [8.5, 9.5],
        'phid': [22.8, 18.0],
        'e': [0.79, 0.64]
    }

    exp_data_dicts = {
        'Be07': Beladjine07,
        'Zh06': Zhou06,
        'Ri95_c': Rice95_coarse,
        'Ri95_m': Rice95_medium,
        'Ri95_f': Rice95_fine,
        'Ch18': Chen18,
        'Wi89_c': Willetts89_coarse,
        'Wi89_m': Willetts89_medium,
        'Wi89_f': Willetts89_fine,
        'Ri20': Rioual20,
        'Go21': Gordon21,
        'Go09': Gordon09
    }

    return exp_data_dicts

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