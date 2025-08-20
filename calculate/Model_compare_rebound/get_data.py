import numpy as np
from scipy.special import erfc
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

    v1_std_list = [0.5577, 0.5507, 0.5656, 0.5823]
    v2_std_list = [0.4673, 0.4732, 0.4890, 0.5326]
    v1_mean_list = [2.7323, 2.7013, 2.7070, 2.7979]
    v1_Ri95_c = np.mean(v1_mean_list)
    v2_mean_list = [1.5775, 1.5039, 1.5184, 1.5898]
    # e = v2/v1
    e_list = []
    e_std_list = []
    for i in range(len(v1_std_list)):
        e = v2_mean_list[i] / v1_mean_list[i]
        e_list.append(e)
        e_std = e * np.sqrt((v2_std_list[i] / v2_mean_list[i])**2 + (v1_std_list[i] / v1_mean_list[i])**2)
        e_std_list.append(e_std)
    Rice95_coarse= {
        'd1': [coarse_d1, coarse_d1, coarse_d1, coarse_d1],
        'v1': v1_mean_list,
        'v1_std': v1_std_list,
        'thd': [13.94, 14.75, 14.73, 15.04],
        'thd_std': [4.88, 4.66, 4.08, 5.76],
        'phid': [20.95, 23.03, 22.55, 25.63],
        'phid_std': [14.25, 14.45, 16.49, 16.06],
        'e': e_list,
        'e_std': e_std_list,
        'Nej': [5.7, 5.09, 5.12, 4.64],
        'Nej_succ': [5.81, 5.46, 5.55, 5.63], # based on the impacts with Nej!=0
        'vn': [0.243, 0.2422, 0.2522, 0.2637],
        'vn_std': [0.2093, 0.1930, 0.2142, 0.2121]
    }

    v1_std_list = [0.7499, 0.6975, 0.6191, 0.6330]
    v2_std_list = [0.6626, 0.6726, 0.6649, 0.6986]
    v1_mean_list = [3.1649, 3.3638, 3.3604, 3.2963]
    v1_Ri95_m = np.mean(v1_mean_list)
    v2_mean_list = [1.7851, 1.8085, 1.9645, 1.8794]
    # e = v2/v1
    e_list = []
    e_std_list = []
    for i in range(len(v1_std_list)):
        e = v2_mean_list[i] / v1_mean_list[i]
        e_list.append(e)
        e_std = e * np.sqrt((v2_std_list[i] / v2_mean_list[i])**2 + (v1_std_list[i] / v1_mean_list[i])**2)
        e_std_list.append(e_std)
    Rice95_medium= {
        'd1': [medium_d1, medium_d1, medium_d1, medium_d1],
        'v1': v1_mean_list,
        'v1_std': v1_std_list,
        'thd': [11.82, 11.47, 11.53, 11.36],
        'thd_std': [3.15, 3.13, 3.47, 3.42],
        'phid': [29.95, 30.31, 31.06, 29.56],
        'phid_std': [21.43, 22.12, 29.02, 22.06],
        'e': e_list,
        'e_std': e_std_list,
        'Nej': [2.68, 3.43, 2.67, 2.76],
        'Nej_succ': [3.22, 3.95, 3.30, 3.71],
        'vn': [0.2544, 0.252, 0.2607, 0.295],
        'vn_std': [0.2340, 0.1984, 0.2486, 0.2458]
    }

    v1_std_list = [0.6659, 0.6375, 0.6715]
    v2_std_list = [0.7859, 0.7966, 0.7919]
    v1_mean_list = [3.8493, 3.7801, 3.7411]
    v1_Ri95_f = np.mean(v1_mean_list)
    v2_mean_list = [1.9929, 2.1281, 2.156]
    # e = v2/v1
    e_list = []
    e_std_list = []
    for i in range(len(v1_std_list)):
        e = v2_mean_list[i] / v1_mean_list[i]
        e_list.append(e)
        e_std = e * np.sqrt((v2_std_list[i] / v2_mean_list[i])**2 + (v1_std_list[i] / v1_mean_list[i])**2)
        e_std_list.append(e_std)
    Rice95_fine= {
        'd1': [fine_d1, fine_d1, fine_d1],
        'v1': v1_mean_list,
        'v1_std': v1_std_list,
        'thd': [10.85, 10.24, 10.46],
        'thd_std': [2.8, 3.06, 2.47],
        'phid': [44.63, 38.03, 37.85],
        'phid_std': [31.31, 29.52, 31.89],
        'e': e_list,
        'e_std': e_std_list,
        'Nej': [1.86, 1.71, 1.58],
        'Nej_succ': [2.20, 2.24, 2.33],
        'vn': [0.2757, 0.2566, 0.2773],
        'vn_std': [0.1942, 0.2249, 0.2588]
    }

    Chen18= {
        'thd': [23.2, 21.8, 21.9, 30.7, 30.6, 37.6, 47.1, 46.6, 46],
        'phid': [32.83, 35.59, 29.90, 43, 33.28, 44.16, 57.62, 50.97, 39.88],
        'e': [0.38, 0.4, 0.45, 0.32, 0.37, 0.27, 0.2, 0.19, 0.22]
    }

    Willetts89_coarse= {
        'd1': [coarse_d1, coarse_d1, coarse_d1, coarse_d1],
        'v1': [3.38, 3.43, 3.18, 3.50],
        'v1_std': [0.89, 0.73, 0.56, 0.52],
        'thd': [12.7, 17.8, 23.2, 27.7],
        'thd_std': [5.1, 6.3, 4.9, 4.8],
        'phid': [19.1, 25.2, 21.4, 27.2],
        'phid_std': [13, 13.3, 13.8, 12.4],
        'e': [0.63, 0.57, 0.54, 0.46],
        'e_std': [0.16, 0.15, 0.16, 0.15],
        'vn': [0.3718, 0.3773, 0.4134, 0.42]
    }

    Willetts89_medium= {
        'd1': [medium_d1, medium_d1, medium_d1, medium_d1],
        'v1': [3.56, 3.99, 4.02, 4.39],
        'v1_std': [0.94, 0.70, 0.77, 0.71],
        'thd': [11.7, 18.2, 21.4, 26.3],
        'thd_std': [5.2, 4.0, 5.2, 4.6],
        'phid': [24.9, 33.4, 33.3, 44.7],
        'phid_std': [14.2, 22.9, 22.9, 24.2],
        'e': [0.61, 0.53, 0.50, 0.40],
        'e_std': [0.14, 0.16, 0.15, 0.15],
        'vn': [0.356, 0.399, 0.3618, 0.439]
    }

    Willetts89_fine= {
        'd1': [fine_d1, fine_d1, fine_d1, fine_d1],
        'v1': [3.61, 4.41, 4.26, 4.35],
        'v1_std': [1.22, 0.77, 0.76, 0.74],
        'thd': [9.5, 15.4, 19.7, 24.9],
        'thd_std': [4.6, 3.6, 3.7, 4.0],
        'phid': [38.8, 42, 42.2, 42.5],
        'phid_std': [25.8, 24.5, 23.4, 26.5],
        'e': [0.57, 0.50, 0.48, 0.46],
        'e_std': [0.19, 0.17, 0.16, 0.18],
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

    g = 9.8*(1.0 - 1.263/2650)
    a_list = [0.126, 0.119, 0.114, 0.138, 0.142]
    bed_list = [0, 5, 10, 12.5, 15]
    thd_list = [bed + 11.5 for bed in bed_list]
    d_bar = 3.2e-4

    Nej_list = []
    for i in range(len(a_list)):
        Nej = a_list[i] * v1_Ri95_c/np.sqrt(g*d_bar)*coarse_d1**3/d_bar**3
        Nej_list.append(Nej)
    Yin21_coarse_Nej_thd= {
        'thd': thd_list,
        'Nej': Nej_list
    }
    Nej_list = []
    for i in range(len(a_list)):
        Nej = a_list[i] * v1_Ri95_m/np.sqrt(g*d_bar)*medium_d1**3/d_bar**3
        Nej_list.append(Nej)
    Yin21_medium_Nej_thd= {
        'thd': thd_list,
        'Nej': Nej_list
    }
    Nej_list = []
    for i in range(len(a_list)):
        Nej = a_list[i] * v1_Ri95_f/np.sqrt(g*d_bar)*fine_d1**3/d_bar**3
        Nej_list.append(Nej)
    Yin21_fine_Nej_thd= {
        'thd': thd_list,
        'Nej': Nej_list
    }

    a = 0.126
    v1_list = np.linspace(1, 10, 100)

    Nej_list = []
    for i in range(len(v1_list)):
        Nej = a * v1_list[i]/np.sqrt(g*d_bar)*coarse_d1**3/d_bar**3
        Nej_list.append(Nej)
    Yin21_coarse_Nej_v1= {
        'v1': v1_list,
        'Nej': Nej_list
    }
    Nej_list = []
    for i in range(len(v1_list)):
        Nej = a * v1_list[i]/np.sqrt(g*d_bar)*medium_d1**3/d_bar**3
        Nej_list.append(Nej)
    Yin21_medium_Nej_v1= {
        'v1': v1_list,
        'Nej': Nej_list
    }
    Nej_list = []
    for i in range(len(v1_list)):
        Nej = a * v1_list[i]/np.sqrt(g*d_bar)*fine_d1**3/d_bar**3
        Nej_list.append(Nej)
    Yin21_fine_Nej_v1= {
        'v1': v1_list,
        'Nej': Nej_list
    }

    v1_hat_list = [17.73, 35.546, 53.43, 71.178, 89.13, 106.95, 124.69]
    vn_list = [0.219, 0.277, 0.308, 0.311, 0.330, 0.335, 0.349]
    v1_list = [v1_hat*np.sqrt(g*d_bar) for v1_hat in v1_hat_list]
    Yin21_vn_v1= {
        'v1': v1_list,
        'vn': vn_list
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
        'Go09': Gordon09,
        'Yi21_c_N_t': Yin21_coarse_Nej_thd,
        'Yi21_m_N_t': Yin21_medium_Nej_thd,
        'Yi21_f_N_t': Yin21_fine_Nej_thd,
        'Yi21_c_N_v': Yin21_coarse_Nej_v1,
        'Yi21_m_N_v': Yin21_medium_Nej_v1,
        'Yi21_f_N_v': Yin21_fine_Nej_v1,
        'Yi21_vn_v1': Yin21_vn_v1
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
    # calculating
    for i in range(sampling_num):
        if bed_type['three_D']:
            phi, e = calculate_rebound_3D(d1, d2_array[i], d3_array[i], th, epsilon, nu)
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
            d2 = bed_type['d_mean']
            d3 = d2
            #if th > 1:
            #    d3 = d2
            #else:
            #    sigma = dist_params['sigma']
            #    mu = dist_params['mu']
            #    dmin = dist_params['d_min']
            #    c = d1 + d2
            #    c1 = d1/c
            #    c2 = d2/c
            #    dc = (1.0 + c2*th - c1)*c/(1.0 - th)
            #    erfc1 = erfc(-(np.log(dc) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
            #    erfc2 = erfc(-(np.log(dc) - mu)/(np.sqrt(2.0)*sigma))
            #    erfc3 = erfc(-(np.log(dmin) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
            #    erfc4 = erfc(-(np.log(dmin) - mu)/(np.sqrt(2.0)*sigma))
            #    d3 = np.exp(mu + sigma**2/2)*(erfc1 - erfc3)/(erfc2 - erfc4)
            #sigma = dist_params['sigma']
            #mu = dist_params['mu']
            #c = d1 + d2
            #c1 = d1/c
            #c2 = d2/c
            #dc = (1.0 + c2*np.tan(th) - c1/np.cos(th))*c/(1.0/np.cos(th) - np.tan(th))
            #erfc1 = erfc(-(np.log(dc) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
            #erfc2 = erfc(-(np.log(dc) - mu)/(np.sqrt(2.0)*sigma))
            #d3 = np.exp(mu + sigma**2/2)*erfc1/erfc2
            if bed_type['three_D']:
                phi, e = calculate_rebound_3D(d1, d2, d3, th, epsilon, nu)
            else:
                phi, e = calculate_rebound_2D(d1, d2, d3, th, epsilon, nu)
            phi_list.append(phi)
            e_list.append(e)
        else:
            phi, e = sample_averaged_rebound(epsilon, nu, d1, th, dist_params, bed_type)
            phi_list.append(phi)
            e_list.append(e)

    return np.array(phi_list), np.array(e_list)