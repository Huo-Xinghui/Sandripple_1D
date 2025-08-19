import numpy as np
from scipy.special import erfc
from scipy import stats
from data_process import generate_truncated_lognormal
from rebound_model import calculate_rebound_2D, calculate_rebound_3D

def sample_averaged_rebound(epsilon, nu, th, sigma, dist_params, bed_type):
    """Sample average rebound angle and restitution coefficient (average over bed PSD)."""
    phi_list = []
    e_list = []
    # sampling
    mu = dist_params['mu']
    d_min = dist_params['d_min']
    d_max = dist_params['d_max']
    sampling_num = dist_params['sampling_num']
    d_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, 100000)
    d50 = np.percentile(d_array, 50)
    d10 = np.percentile(d_array, 10)
    d90 = np.percentile(d_array, 90)
    d_mean = np.mean(d_array)
    d1 = d_mean
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

def sample_averaged_rebound_simple(epsilon, nu, th, sigma, dist_params, bed_type):
    """Sample average rebound angle and restitution coefficient (average over bed PSD)."""
    phi_list = []
    e_list = []
    # sampling
    mu = dist_params['mu']
    d_min = dist_params['d_min']
    d_max = dist_params['d_max']
    sampling_num = dist_params['sampling_num']
    d_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, 100000)
    d50 = np.percentile(d_array, 50)
    d10 = np.percentile(d_array, 10)
    d90 = np.percentile(d_array, 90)
    d_mean = np.mean(d_array)
    d1 = d_mean
    d2_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, sampling_num)
    d3_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, sampling_num)
    # calculating
    for i in range(sampling_num):
        if bed_type['three_D']:
            d2 = d2_array[i]
            #if th > 1:
            #    d3 = d2
            #else:
            #    dc = d2*(1.0 + th)/(1.0 - th)
            #    erfc1 = erfc(-(np.log(dc) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
            #    erfc2 = erfc(-(np.log(dc) - mu)/(np.sqrt(2.0)*sigma))
            #    erfc3 = erfc(-(np.log(dmin) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
            #    erfc4 = erfc(-(np.log(dmin) - mu)/(np.sqrt(2.0)*sigma))
            #    d3 = np.exp(mu + sigma**2/2)*(erfc1 - erfc3)/(erfc2 - erfc4)
            if bed_type['monodisperse']:
                d3 = d2
            else:
                c = d1 + d2
                c1 = d1/c
                c2 = d2/c
                dc = (1.0 + c2*np.tan(th) - c1/np.cos(th))*c/(1.0/np.cos(th) - np.tan(th))
                erfc1 = erfc(-(np.log(dc) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                erfc2 = erfc(-(np.log(dc) - mu)/(np.sqrt(2.0)*sigma))
                erfc3 = erfc(-(np.log(d_min) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                erfc4 = erfc(-(np.log(d_min) - mu)/(np.sqrt(2.0)*sigma))
                d3 = np.exp(mu + sigma**2/2)*(erfc1 - erfc3)/(erfc2 - erfc4)
                #erfc1 = erfc(-(np.log(dc) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                #erfc2 = erfc(-(np.log(dc) - mu)/(np.sqrt(2.0)*sigma))
                #d3 = np.exp(mu + sigma**2/2)*erfc1/erfc2
            phi, e = calculate_rebound_3D(d1, d2, d3, th, epsilon, nu)
        else:
            d2 = d2_array[i]
            #if th > 1:
            #    d3 = d2
            #else:
            #    dc = d2*(1.0 + th)/(1.0 - th)
            #    erfc1 = erfc(-(np.log(dc) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
            #    erfc2 = erfc(-(np.log(dc) - mu)/(np.sqrt(2.0)*sigma))
            #    erfc3 = erfc(-(np.log(d_min) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
            #    erfc4 = erfc(-(np.log(d_min) - mu)/(np.sqrt(2.0)*sigma))
            #    d3 = np.exp(mu + sigma**2/2)*(erfc1 - erfc3)/(erfc2 - erfc4)
            if bed_type['monodisperse']:
                d3 = d2
            else:
                c = d1 + d2
                c1 = d1/c
                c2 = d2/c
                dc = (1.0 + c2*np.tan(th) - c1/np.cos(th))*c/(1.0/np.cos(th) - np.tan(th))
                erfc1 = erfc(-(np.log(dc) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                erfc2 = erfc(-(np.log(dc) - mu)/(np.sqrt(2.0)*sigma))
                erfc3 = erfc(-(np.log(d_min) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                erfc4 = erfc(-(np.log(d_min) - mu)/(np.sqrt(2.0)*sigma))
                d3 = np.exp(mu + sigma**2/2)*(erfc1 - erfc3)/(erfc2 - erfc4)
                #erfc1 = erfc(-(np.log(dc) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                #erfc2 = erfc(-(np.log(dc) - mu)/(np.sqrt(2.0)*sigma))
                #d3 = np.exp(mu + sigma**2/2)*erfc1/erfc2
            phi, e = calculate_rebound_2D(d1, d2, d3, th, epsilon, nu)
        phi_list.append(phi)
        e_list.append(e)
    phi_avg = np.mean(phi_list)
    e_avg = np.mean(e_list)

    return phi_avg, e_avg

def get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type):
    """Get model data array for rebound angle and restitution coefficient."""
    phi_list = []
    e_list = []
    d_min = dist_params['d_min']
    d_max = dist_params['d_max']
    mu = dist_params['mu']
    for sigma in sigma_array:
        if bed_type['simple']:
            d_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, 100000)
            d50 = np.percentile(d_array, 50)
            d10 = np.percentile(d_array, 10)
            d90 = np.percentile(d_array, 90)
            d_mean = np.mean(d_array)
            d1 = d_mean
            d2 = d_mean
            #if th > 1:
            #    d3 = d2
            #else:
            #    dc = d2*(1.0 + th)/(1.0 - th)
            #    erfc1 = erfc(-(np.log(dc) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
            #    erfc2 = erfc(-(np.log(dc) - mu)/(np.sqrt(2.0)*sigma))
            #    erfc3 = erfc(-(np.log(dmin) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
            #    erfc4 = erfc(-(np.log(dmin) - mu)/(np.sqrt(2.0)*sigma))
            #    d3 = np.exp(mu + sigma**2/2)*(erfc1 - erfc3)/(erfc2 - erfc4)
            if bed_type['monodisperse']:
                d3 = d2
            else:
                c = d1 + d2
                c1 = d1/c
                c2 = d2/c
                dc = (1.0 + c2*np.tan(th) - c1/np.cos(th))*c/(1.0/np.cos(th) - np.tan(th))
                erfc1 = erfc(-(np.log(dc) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                erfc2 = erfc(-(np.log(dc) - mu)/(np.sqrt(2.0)*sigma))
                erfc3 = erfc(-(np.log(d_min) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                erfc4 = erfc(-(np.log(d_min) - mu)/(np.sqrt(2.0)*sigma))
                d3 = np.exp(mu + sigma**2/2)*(erfc1 - erfc3)/(erfc2 - erfc4)
                #erfc1 = erfc(-(np.log(dc) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                #erfc2 = erfc(-(np.log(dc) - mu)/(np.sqrt(2.0)*sigma))
                #d3 = np.exp(mu + sigma**2/2)*erfc1/erfc2
            if bed_type['three_D']:
                phi, e = calculate_rebound_3D(d1, d2, d3, th, epsilon, nu)
            else:
                phi, e = calculate_rebound_2D(d1, d2, d3, th, epsilon, nu)
            phi_list.append(phi)
            e_list.append(e)
            #phi, e = sample_averaged_rebound_simple(epsilon, nu, th, sigma, dist_params, bed_type)
            #phi_list.append(phi)
            #e_list.append(e)
        else:
            phi, e = sample_averaged_rebound(epsilon, nu, th, sigma, dist_params, bed_type)
            phi_list.append(phi)
            e_list.append(e)

    return np.array(phi_list), np.array(e_list)