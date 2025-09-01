import numpy as np
from scipy.special import erfc
from scipy import stats
from data_process import generate_truncated_lognormal
from rebound_model import calculate_rebound_2D, calculate_rebound_3D, calculate_rebound_simple

def sample_averaged_rebound(epsilon, nu, th, sigma, dist_params, bed_type):
    """Sample average rebound angle and restitution coefficient (average over bed PSD)."""
    phi_list = []
    e_list = []
    ez_list = []
    ex_list = []
    ecx_list = []
    ecz_list = []
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
            phi, e, ecx, ecz, ez, ex = calculate_rebound_3D(d1, d2_array[i], d3_array[i], th, epsilon, nu)
        else:
            phi, e, ecx, ecz, ez, ex = calculate_rebound_2D(d1, d2_array[i], d3_array[i], th, epsilon, nu)
        phi_list.append(phi)
        e_list.append(e)
        ecx_list.append(ecx)
        ecz_list.append(ecz)
        ez_list.append(ez)
        ex_list.append(ex)
    phi_avg = np.mean(phi_list)
    e_avg = np.mean(e_list)
    ecx_avg = np.mean(ecx_list)
    ecz_avg = np.mean(ecz_list)
    ez_avg = np.mean(ez_list)
    ex_avg = np.mean(ex_list)

    return phi_avg, e_avg, ecx_avg, ecz_avg, ez_avg, ex_avg

def sample_averaged_rebound_simple(epsilon, nu, th, sigma, dist_params, bed_type):
    """Sample average rebound angle and restitution coefficient (average over bed PSD)."""
    phi_list = []
    e_list = []
    ecx_list = []
    ecz_list = []
    ez_list = []
    ex_list = []
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
    for d2, d3 in zip(d2_array, d3_array):
        if bed_type['three_D']:
            if bed_type['monodisperse']:
                d2 = d_mean
            else:
                d3 = d_mean
                #c = d1 + d2
                #c1 = d1/c
                #c2 = d2/c
                #dc = (1.0 + c2*np.tan(th) - c1/np.cos(th))*c/(1.0/np.cos(th) - np.tan(th))
                #erfc1 = erfc(-(np.log(dc) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                #erfc2 = erfc(-(np.log(dc) - mu)/(np.sqrt(2.0)*sigma))
                #erfc3 = erfc(-(np.log(d_min) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                #erfc4 = erfc(-(np.log(d_min) - mu)/(np.sqrt(2.0)*sigma))
                #d3 = np.exp(mu + sigma**2/2)*(erfc1 - erfc3)/(erfc2 - erfc4)
            phi, e, ecx, ecz, ez, ex = calculate_rebound_3D(d1, d2, d3, th, epsilon, nu)
        else:
            if bed_type['monodisperse']:
                d2 = d_mean
            else:
                d3 = d_mean
                #c = d1 + d2
                #c1 = d1/c
                #c2 = d2/c
                #dc = (1.0 + c2*np.tan(th) - c1/np.cos(th))*c/(1.0/np.cos(th) - np.tan(th))
                #erfc1 = erfc(-(np.log(dc) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                #erfc2 = erfc(-(np.log(dc) - mu)/(np.sqrt(2.0)*sigma))
                #erfc3 = erfc(-(np.log(d_min) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                #erfc4 = erfc(-(np.log(d_min) - mu)/(np.sqrt(2.0)*sigma))
                #d3 = np.exp(mu + sigma**2/2)*(erfc1 - erfc3)/(erfc2 - erfc4)
            phi, e, ecx, ecz, ez, ex = calculate_rebound_2D(d1, d2, d3, th, epsilon, nu)
        phi_list.append(phi)
        e_list.append(e)
        ecx_list.append(ecx)
        ecz_list.append(ecz)
        ez_list.append(ez)
        ex_list.append(ex)
    phi_avg = np.mean(phi_list)
    e_avg = np.mean(e_list)
    ecx_avg = np.mean(ecx_list)
    ecz_avg = np.mean(ecz_list)
    ez_avg = np.mean(ez_list)
    ex_avg = np.mean(ex_list)

    return phi_avg, e_avg, ecx_avg, ecz_avg, ez_avg, ex_avg

def get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type):
    """Get model data array for rebound angle and restitution coefficient."""
    phi_list = []
    e_list = []
    ecx_list = []
    ecz_list = []
    ez_list = []
    ex_list = []
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
            d_mode = stats.mode(d_array).mode
            d1 = d_mean
            if bed_type['monodisperse']:
                d2 = d_mean
                #d3 = d2
                dc = (d2*(np.cos(th) + np.sin(th)) - d1*(1.0 - np.cos(th)))/(1.0 - np.sin(th))
                dmin = d_min
                dmax = d_max
                m = mu
                s = sigma
                #d_max = min(d_max, dc)
                erfc1 = erfc((np.log(dmin) - m - s**2)/(np.sqrt(2.0)*s))
                erfc2 = erfc((np.log(dmin) - m)/(np.sqrt(2.0)*s))
                erfc3 = erfc((np.log(dmax) - m - s**2)/(np.sqrt(2.0)*s))
                erfc4 = erfc((np.log(dmax) - m)/(np.sqrt(2.0)*s))
                aa = (erfc1 - erfc3)/(erfc2 - erfc4)
                dmax = min(d_max, dc)
                erfc1 = erfc((-np.log(dmax) + m - s**2)/(np.sqrt(2.0)*s))
                erfc2 = erfc((np.log(dmin) - m)/(np.sqrt(2.0)*s))
                erfc3 = erfc((-np.log(dmin) + m - s**2)/(np.sqrt(2.0)*s))
                erfc4 = erfc((np.log(dmax) - m)/(np.sqrt(2.0)*s))
                bb = (erfc1 - erfc3)/(erfc2 - erfc4)
                #d2d3 = np.exp(m + s**2/2)*(erfc1 - erfc3)/(erfc2 - erfc4)
                d2d3 = np.exp(s**2)*aa*bb
                #d3 = d_mean
                d3 = d2 / d2d3
            else:
                d3 = d50
                #d3 = d2
                dc1 = (d3*(1.0 - np.sin(th)) + d1*(1.0 - np.cos(th)))/(np.cos(th) + np.sin(th))
                d_min = max(d_min, dc1)
                erfc1 = erfc(-(np.log(d_max) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                erfc2 = erfc(-(np.log(d_max) - mu)/(np.sqrt(2.0)*sigma))
                erfc3 = erfc(-(np.log(d_min) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
                erfc4 = erfc(-(np.log(d_min) - mu)/(np.sqrt(2.0)*sigma))
                d2 = np.exp(mu + sigma**2/2)*(erfc1 - erfc3)/(erfc2 - erfc4)
                #d2 = d_mean
            phi, e, ecx, ecz, ez, ex = calculate_rebound_2D(d1, d2, d3, th, epsilon, nu)
            #phi, e, ecx, ecz, ez, ex = calculate_rebound_simple(d1, d2, d3, th, epsilon, nu, sigma)
            #phi, e, ecx, ecz, ez, ex = sample_averaged_rebound_simple(epsilon, nu, th, sigma, dist_params, bed_type)
            phi_list.append(phi)
            e_list.append(e)
            ecx_list.append(ecx)
            ecz_list.append(ecz)
            ez_list.append(ez)
            ex_list.append(ex)
        else:
            phi, e, ecx, ecz, ez, ex = sample_averaged_rebound(epsilon, nu, th, sigma, dist_params, bed_type)
            phi_list.append(phi)
            e_list.append(e)
            ecx_list.append(ecx)
            ecz_list.append(ecz)
            ez_list.append(ez)
            ex_list.append(ex)

    return np.array(phi_list), np.array(e_list), np.array(ecx_list), np.array(ecz_list), np.array(ez_list), np.array(ex_list)

def get_d_data_array(th, sigma_array, dist_params):
    """Get model d2, d3 array."""
    d1_list = []
    d2_list = []
    d3_list = []
    d_min = dist_params['d_min']
    d_max = dist_params['d_max']
    mu = dist_params['mu']
    for sigma in sigma_array:
        d_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, 100000)
        d50 = np.percentile(d_array, 50)
        d10 = np.percentile(d_array, 10)
        d90 = np.percentile(d_array, 90)
        d_mean = np.mean(d_array)
        d_mode = stats.mode(d_array).mode
        d1 = d90
        d3 = d50
        dc1 = (d3*(1.0 - np.sin(th)) + d1*(1.0 - np.cos(th)))/(np.cos(th) + np.sin(th))
        d_min = max(d_min, dc1)
        erfc1 = erfc(-(np.log(d_max) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
        erfc2 = erfc(-(np.log(d_max) - mu)/(np.sqrt(2.0)*sigma))
        erfc3 = erfc(-(np.log(d_min) - mu - sigma**2)/(np.sqrt(2.0)*sigma))
        erfc4 = erfc(-(np.log(d_min) - mu)/(np.sqrt(2.0)*sigma))
        d2 = np.exp(mu + sigma**2/2)*(erfc1 - erfc3)/(erfc2 - erfc4)
        d1_list.append(d1)
        d2_list.append(d2)
        d3_list.append(d3)
    return np.array(d1_list), np.array(d2_list), np.array(d3_list)