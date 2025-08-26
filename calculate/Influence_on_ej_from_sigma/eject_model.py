import numpy as np
from scipy.integrate import quad
from scipy.stats import norm
from scipy.special import erfc
from data_process import generate_truncated_lognormal

def calculate_x_min(alpha, beta, d1, d2, d3, th):
    d13 = (d1 + d3)/2
    d23 = (d2 + d3)/2
    cos0 = (d13**2 + d23**2 - 1)/(2*d13*d23)
    a = np.arccos(cos0)
    theta_c = np.pi/2 - a
    l_r = (1 + d23**2 - d13**2)/(2*d23)
    if th <= theta_c:
        x_min = d13/np.sin(th) - d23
    else:
        x_min = np.sqrt(1 - l_r**2)/np.tan(th) - l_r
        x_max_try = 1/np.tan(th)
        x_min_try = x_min
        attempt_num = 100
        delta_x = (x_max_try - x_min_try)/attempt_num
        for n in range(0, attempt_num):
            x = x_min_try + n*delta_x
            x_sin_square = x**2*(np.sin(th))**2
            x_cos = x*np.cos(th)
            neq_LHS = 1/(1 + beta/alpha)
            neq_RHS = x_sin_square - x_cos*np.sqrt(1.0 - x_sin_square)
            if neq_LHS >= neq_RHS:
                break
        x_min = x

    return x_min

def calculate_x_min_3D(alpha, beta, d1, d2, d3, th):
    d13 = (d1 + d3)/2
    d23 = (d2 + d3)/2
    sin_zeta_max = (d3/2)/d23
    zeta_max = np.arcsin(sin_zeta_max)
    zeta = np.random.uniform(0, zeta_max)
    d13_temp = d13*np.sqrt(1.0 - (np.sin(zeta)*d23/d13)**2)
    d23_temp = (d2 + d3*np.sqrt(1.0 - (np.sin(zeta)*d23/d13)**2))/2
    while d23_temp - d13_temp >= 1.0:
        zeta = np.random.uniform(0, zeta_max)
        d13_temp = d13*np.sqrt(1.0 - (np.sin(zeta)*d23/d13)**2)
        d23_temp = (d2 + d3*np.sqrt(1.0 - (np.sin(zeta)*d23/d13)**2))/2
        print("d23-d13<1")
    d13 = d13_temp
    d23 = d23_temp
    cos0 = (d13**2 + d23**2 - 1)/(2*d13*d23)
    a = np.arccos(cos0)
    theta_c = np.pi/2 - a
    l_r = (1 + d23**2 - d13**2)/(2*d23)
    if th <= theta_c:
        x_min = d13/np.sin(th) - d23
    else:
        x_min = np.sqrt(1 - l_r**2)/np.tan(th) - l_r
        x_max_try = 1/np.tan(th)
        x_min_try = x_min
        attempt_num = 100
        delta_x = (x_max_try - x_min_try)/attempt_num
        for n in range(0, attempt_num):
            x = x_min_try + n*delta_x
            x_sin_square = x**2*(np.sin(th))**2
            x_cos = x*np.cos(th)
            neq_LHS = 1/(1 + beta/alpha)
            neq_RHS = x_sin_square - x_cos*np.sqrt(1.0 - x_sin_square)
            if neq_LHS >= neq_RHS:
                break
        x_min = x
    return x_min

def calculate_x_max(alpha, beta, th):
    x_max_try = 1/np.sin(th)
    x_min_try = 1/np.tan(th)
    attempt_num = 100
    delta_x = (x_max_try - x_min_try)/attempt_num
    xx = x_min_try
    for n in range(1, attempt_num + 1):
        x = x_min_try + n*delta_x
        x_sin_square = x**2*np.sin(th)**2
        x_cos = x*np.cos(th)
        neq_LHS = 1/(1 + beta/alpha)
        neq_RHS = x_sin_square - x_cos*np.sqrt(1 - x_sin_square)
        if neq_LHS < neq_RHS:
            break
        xx = x
    x_max = xx
    return x_max

def rebound_res_eq(x, alpha, beta, th):
    e_vx = -alpha*np.cos(th) + (alpha + beta)*x**2*np.sin(th)**2*np.cos(th) + (alpha + beta)*x*np.sin(th)**2*np.sqrt(1 - x**2*np.sin(th)**2)
    e_vz = alpha*np.sin(th) - (alpha + beta)*x**2*np.sin(th)**3 + (alpha + beta)*x*np.sin(th)*np.cos(th)*np.sqrt(1 - x**2*np.sin(th)**2)
    res = np.sqrt(e_vx**2 + e_vz**2)
    return res

def calculate_rebound(d1, d2, d3, th, epsilon, nu, mono):
    # normalize diameters
    D = (d1 + d2)/2
    D1 = d1/D
    D2 = d2/D
    D3 = d3/D
    # restitution coefficient
    mu = epsilon*D1**3/(D1**3 + epsilon*D2**3)
    alpha = (1 + epsilon)/(1 + mu) - 1
    beta = 1 - (2/7)*(1 - nu)/(1 + mu)
    if mono:
        x_min = calculate_x_min(alpha, beta, D1, D2, D3, th)
    else:
        x_min = calculate_x_min_3D(alpha, beta, D1, D2, D3, th)
    x_max = calculate_x_max(alpha, beta, th)
    if x_min < x_max:
        e, error = quad(rebound_res_eq, x_min, x_max, args=(alpha, beta, th))
        e /= (x_max - x_min)
    else:
        # exchange d2 and d3
        d_temp = d2
        d2 = d3
        d3 = d_temp
        D = (d1 + d2)/2
        D1 = d1/D
        D2 = d2/D
        D3 = d3/D
        # restitution coefficient
        mu = epsilon*D1**3/(D1**3 + epsilon*D2**3)
        alpha = (1 + epsilon)/(1 + mu) - 1
        beta = 1 - (2/7)*(1 - nu)/(1 + mu)
        if mono:
            x_min = calculate_x_min(alpha, beta, D1, D2, D3, th)
        else:
            x_min = calculate_x_min_3D(alpha, beta, D1, D2, D3, th)
        x_max = calculate_x_max(alpha, beta, th)
        if x_min > x_max:
            e = rebound_res_eq(x_max, alpha, beta, th)
        else:
            e, error = quad(rebound_res_eq, x_min, x_max, args=(alpha, beta, th))
            e /= (x_max - x_min)
    return e

def f_erfc1(mu, sigma, d_c):
    f1 = erfc((np.log(d_c) - mu) / (np.sqrt(2) * sigma))
    return f1

def f_erfc2(mu, sigma, d_c):
    f2 = erfc((np.log(d_c) - mu - sigma**2) / (np.sqrt(2) * sigma))
    return f2

def calculate_Ec(mu, sigma, d_min, d_max, C_Ec):
    a = np.log(d_min)
    b_max = np.log(d_max)
    x_min = (a - mu) / sigma
    x_max = (b_max - mu) / sigma
    b = 3.0*np.log(2.0)*(2.0 - np.log(2.0))
    erfc1 = erfc(-(x_max - b*sigma)/np.sqrt(2))
    erfc2 = erfc(-(x_min - b*sigma)/np.sqrt(2))
    erfc3 = erfc(-x_max/np.sqrt(2))
    erfc4 = erfc(-x_min/np.sqrt(2))
    erfc5 = erfc(-(x_max - sigma)/np.sqrt(2))
    erfc6 = erfc(-(x_min - sigma)/np.sqrt(2))
    A = np.exp((b - 1)*sigma**2/2)*((erfc1 - erfc2) / (erfc5 - erfc6))**(1/b)
    dc = np.exp(mu + b*sigma**2 / 2) * ((erfc1 - erfc2)/(erfc3 - erfc4))**(1/b)
    Ec = C_Ec * dc**3
    return Ec

def calculate_eject(th, v1, gamma, d_dict, physical_dict, bed_type, dist_params):
    """calculate ejection of one impact"""
    d_min = dist_params['d_min']
    d_max = dist_params['d_max']
    mu = dist_params['mu']
    sigma = dist_params['sigma']
    d1 = d_dict['d1']
    d2 = d_dict['d2']
    d3 = d_dict['d3']
    zc = d_dict['zc']
    epsilon = physical_dict['epsilon']
    nu = physical_dict['nu']
    rho = physical_dict['rho']
    g = physical_dict['g']
    mono = bed_type['monodisperse']

    e = calculate_rebound(d1, d2, d3, th, epsilon, nu, mono)
    m1 = rho*(d1**3*np.pi/6)
    m2 = rho*(d2**3*np.pi/6)
    if bed_type['monodisperse']:
        Ec = m2*g*zc
    else:
        C_Ec = rho*g*np.pi/6*zc
        #Ec_min = C_Ec*d_min**3
        #Ec_max = C_Ec*d_max**3
        #mu_Ec = np.log(C_Ec) + 3.0*mu
        #sigma_Ec = 3.0*sigma
        #z_min = (np.log(Ec_min) - mu_Ec - sigma_Ec**2) / sigma_Ec
        #z_max = (np.log(Ec_max) - mu_Ec - sigma_Ec**2) / sigma_Ec
        #P_Ec = norm.cdf((np.log(Ec_max) - mu_Ec) / sigma_Ec) - norm.cdf((np.log(Ec_min) - mu_Ec) / sigma_Ec)
        #Ec = (C_Ec*np.exp(3.0*mu + 9.0*sigma**2/2.0) * (norm.cdf(z_max) - norm.cdf(z_min))) / P_Ec
        Ec = calculate_Ec(mu, sigma, d_min, d_max, C_Ec)
    Ee1 = 0.5*m1*v1**2
    Ei = (1.0 - e**2)*Ee1
    k_max = Ei/Ec
    lambda_ej = 2.0*np.log(k_max)
    mu_ej = np.log((1.0 - e**2)*Ee1) - lambda_ej*np.log(2.0)
    sigma_ej = np.sqrt(lambda_ej)*np.log(2.0)
    En_bar = np.exp(mu_ej + sigma_ej**2/2)
    Ec_ej = Ec
    Nej = gamma*(1.0 - e**2)*Ee1/En_bar/2.0*erfc((np.log(Ec_ej) - mu_ej)/(np.sqrt(2.0)*sigma_ej))
    erfc2 = erfc((np.log(Ec_ej) - mu_ej)/(np.sqrt(2.0)*sigma_ej))
    erfc3 = erfc((np.log(Ec_ej) - mu_ej - sigma_ej**2)/(np.sqrt(2.0)*sigma_ej))
    E_mean = erfc3/erfc2 * np.exp(mu_ej + sigma_ej**2/2)
    if bed_type['monodisperse']:
        erfc1 = erfc((np.log(Ec_ej) - mu_ej - sigma_ej**2/2)/(np.sqrt(2.0)*sigma_ej))
        vn_mean = erfc1/erfc2 * np.sqrt(2.0/m2)*np.exp(mu_ej/2 + sigma_ej**2/8)
    else:
        #if Nej > 0:
        #    nej = int(np.round(Nej))
        #    vn_list = []
        #    for _ in range(nej):
        #        dn = generate_truncated_lognormal(mu, sigma, d_min, d_max, 1)[0]
        #        mn = rho*(dn**3*np.pi/6)
        #        En = 0.0
        #        while En < Ec_ej:
        #            En = np.random.lognormal(mu_ej, sigma_ej, size=1)[0]
        #        vn = np.sqrt(2.0*En/mn)
        #        vn_list.append(vn)
        #    vn_mean = np.mean(vn_list)
        #else:
        #    vn_mean = 0.0
        erfc1 = erfc((np.log(Ec_ej) - mu_ej - sigma_ej**2/2)/(np.sqrt(2.0)*sigma_ej))
        #vn_mean = erfc1/erfc2 * np.sqrt(2.0/mn)*np.exp(mu_ej/2 + sigma_ej**2/8)
        vn_mean = erfc1/erfc2 * np.sqrt(12.0/(rho*np.pi))*np.exp(mu_ej/2 - 3/2*mu + sigma_ej**2/8 + 9/8*sigma**2)
    return Nej, vn_mean, E_mean, Ei