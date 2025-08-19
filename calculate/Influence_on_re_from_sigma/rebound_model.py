import numpy as np
from scipy.integrate import quad
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
    #cos_xi = (d23**2 + d24**2 - d34**2)/(2*d23*d24)
    #xi = np.arccos(cos_xi)
    #K = d1*(d2 - d4)+d2*(d2 + d4)
    #if K < 0:
    #    zeta_max = 0
    #elif np.sqrt(K) <= 2*d24:
    #    zeta_max = xi - np.arccos(np.sqrt(K)/(2*d24))
    #else:
    #    zeta_max = xi
    sin_zeta_max = (d3/2)/d23
    zeta_max = np.arcsin(sin_zeta_max)
    zeta = np.random.uniform(0, zeta_max)
    #zeta = zeta_max/2
    d13_temp = d13*np.sqrt(1.0 - (np.sin(zeta)*d23/d13)**2)
    d23_temp = d23*np.cos(zeta)
    while d23_temp - d13_temp >= 1.0:
        zeta = np.random.uniform(0, zeta_max)
        #zeta = zeta_max/2
        d13_temp = d13*np.sqrt(1.0 - (np.sin(zeta)*d23/d13)**2)
        d23_temp = d23*np.cos(zeta)
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

def calculate_x_max_try(alpha, beta, d1, d2, d3, d4, th):
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
    x_max_try1 = xx
    d14 = (d1 + d4)/2
    d24 = (d2 + d4)/2
    a = (1 + d24**2 - d14**2)/(2*d24)
    angle = np.arccos(a)
    h = np.sin(angle)
    b = h/np.tan(th)
    x_max_try2 = a + b
    x_max = min(x_max_try1, x_max_try2)
    return x_max

def rebound_angle_eq(x, alpha, beta, th):
    e_vx = -alpha*np.cos(th) + (alpha + beta)*x**2*np.sin(th)**2*np.cos(th) + (alpha + beta)*x*np.sin(th)**2*np.sqrt(1 - x**2*np.sin(th)**2)
    e_vz = alpha*np.sin(th) - (alpha + beta)*x**2*np.sin(th)**3 + (alpha + beta)*x*np.sin(th)*np.cos(th)*np.sqrt(1 - x**2*np.sin(th)**2)
    theta2 = np.arctan2(e_vz, e_vx)
    return theta2

def rebound_res_eq(x, alpha, beta, th):
    e_vx = -alpha*np.cos(th) + (alpha + beta)*x**2*np.sin(th)**2*np.cos(th) + (alpha + beta)*x*np.sin(th)**2*np.sqrt(1 - x**2*np.sin(th)**2)
    e_vz = alpha*np.sin(th) - (alpha + beta)*x**2*np.sin(th)**3 + (alpha + beta)*x*np.sin(th)*np.cos(th)*np.sqrt(1 - x**2*np.sin(th)**2)
    res = np.sqrt(e_vx**2 + e_vz**2)
    return res

def calculate_rebound_2D(d1, d2, d3, th, epsilon, nu):
    # normalize diameters
    D = (d1 + d2)/2
    D1 = d1/D
    D2 = d2/D
    D3 = d3/D
    # restitution coefficient
    mu = epsilon*D1**3/(D1**3 + epsilon*D2**3)
    alpha = (1 + epsilon)/(1 + mu) - 1
    beta = 1 - (2/7)*(1 - nu)/(1 + mu)
    x_min = calculate_x_min(alpha, beta, D1, D2, D3, th)
    x_max = calculate_x_max(alpha, beta, th)
    if x_min < x_max:
        phi, error = quad(rebound_angle_eq, x_min, x_max, args=(alpha, beta, th))
        phi /= (x_max - x_min)
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
        x_min = calculate_x_min(alpha, beta, D1, D2, D3, th)
        x_max = calculate_x_max(alpha, beta, th)
        if x_min > x_max:
            phi = rebound_angle_eq(x_max, alpha, beta, th)
            e = rebound_res_eq(x_max, alpha, beta, th)
        else:
            phi, error = quad(rebound_angle_eq, x_min, x_max, args=(alpha, beta, th))
            phi /= (x_max - x_min)
            e, error = quad(rebound_res_eq, x_min, x_max, args=(alpha, beta, th))
            e /= (x_max - x_min)
    return phi, e

def calculate_rebound_3D(d1, d2, d3, th, epsilon, nu):
    # normalize diameters
    D = (d1 + d2)/2
    D1 = d1/D
    D2 = d2/D
    D3 = d3/D
    # restitution coefficient
    mu = epsilon*D1**3/(D1**3 + epsilon*D2**3)
    alpha = (1 + epsilon)/(1 + mu) - 1
    beta = 1 - (2/7)*(1 - nu)/(1 + mu)
    x_min = calculate_x_min_3D(alpha, beta, D1, D2, D3, th)
    x_max = calculate_x_max(alpha, beta, th)
    if x_min < x_max:
        phi, error = quad(rebound_angle_eq, x_min, x_max, args=(alpha, beta, th))
        phi /= (x_max - x_min)
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
        x_min = calculate_x_min_3D(alpha, beta, D1, D2, D3, th)
        x_max = calculate_x_max(alpha, beta, th)
        if x_min > x_max:
            phi = rebound_angle_eq(x_max, alpha, beta, th)
            e = rebound_res_eq(x_max, alpha, beta, th)
        else:
            phi, error = quad(rebound_angle_eq, x_min, x_max, args=(alpha, beta, th))
            phi /= (x_max - x_min)
            e, error = quad(rebound_res_eq, x_min, x_max, args=(alpha, beta, th))
            e /= (x_max - x_min)
    return phi, e