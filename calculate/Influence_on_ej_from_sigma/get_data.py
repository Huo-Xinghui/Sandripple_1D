import numpy as np
from data_process import generate_truncated_lognormal
from eject_model import calculate_eject

def sample_averaged_eject(th, v1, gamma, d1, physical_dict, bed_type, dist_params):
    """Sample average Nej and vn (average over bed PSD)."""
    Nej_list = []
    vn_list = []
    Eej_list = []
    E1_list = []
    # sampling
    mu = dist_params['mu']
    sigma = dist_params['sigma']
    d_min = dist_params['d_min']
    d_max = dist_params['d_max']
    sampling_num = dist_params['sampling_num']
    d2_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, sampling_num)
    d3_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, sampling_num)
    zc = bed_type['zc']
    # calculating
    for i in range(sampling_num):
        d_dict = {
            'd1': d1,
            'd2': d2_array[i],
            'd3': d3_array[i],
            'zc': zc
        }
        E1 = 0.5 * physical_dict['rho'] * (d1 ** 3 * np.pi / 6) * v1 ** 2
        Nej, vn, Eej, Ei = calculate_eject(th, v1, gamma, d_dict, physical_dict, bed_type, dist_params)
        Nej_list.append(Nej)
        vn_list.append(vn)
        Eej_list.append(Eej)
        E1_list.append(Ei)
    Nej_mean = np.mean(Nej_list)
    vn_mean = np.mean(vn_list)
    Eej_mean = np.mean(Eej_list)
    E1_mean = np.mean(E1_list)

    return Nej_mean, vn_mean, Eej_mean, E1_mean

def get_model_data_array(v1, th, sigma_array, dist_params, bed_type, physical_dict):
    N_list = []
    v_list = []
    E_list = []
    E1_list = []
    d_min = dist_params['d_min']
    d_max = dist_params['d_max']
    mu = dist_params['mu']
    gamma = physical_dict['gamma']
    g = physical_dict['g']
    v1hat = v1
    for sigma in sigma_array:
        dist_params['sigma'] = sigma
        d_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, 100000)
        d50 = np.percentile(d_array, 50)
        d90 = np.percentile(d_array, 90)
        d_mean = np.mean(d_array)
        v1 = v1hat * np.sqrt(g * d_mean)
        if bed_type['monodisperse']:
            if bed_type['d2'] == 50:
                d2 = d50
            elif bed_type['d2'] == 90:
                d2 = d90
            else:
                d2 = d_mean
            d_dict = {
                'd1': d_mean,
                'd2': d2,
                'd3': d2,
                'zc': d2
            }
            Nej, vej, Eej, Ei = calculate_eject(th, v1, gamma, d_dict, physical_dict, bed_type, dist_params)
            N_list.append(Nej)
            v_list.append(vej)
            E_list.append(Eej)
            E1_list.append(Ei)
        else:
            bed_type['zc'] = d90
            Nej, vej, Eej, Ei = sample_averaged_eject(th, v1, gamma, d_mean, physical_dict, bed_type, dist_params)
            N_list.append(Nej)
            v_list.append(vej)
            E_list.append(Eej)
            E1_list.append(Ei)

    return np.array(N_list), np.array(v_list), np.array(E_list), np.array(E1_list)