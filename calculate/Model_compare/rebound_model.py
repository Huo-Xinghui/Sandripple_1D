import numpy as np
from scipy.integrate import quad

def calculate_x_min(d1, d2, d3, th):
    d13 = (d1 + d3)/2
    d23 = (d2 + d3)/2
    cos0 = (d13**2 + d23**2 - 1)/(2*d13*d23)
    alpha = np.arccos(cos0)
    theta_c = np.pi/2 - alpha
    l_r = (1 + d23**2 - d13**2)/(2*d23)
    if th <= theta_c:
        x_min = d13/np.sin(th) - d23
    else:
        x_min = np.sqrt(1 - l_r**2)/np.tan(th) - l_r
    return x_min

def calculate_x_min_3D(d1, d2, d3, d4, th):
    d13 = (d1 + d3)/2
    d23 = (d2 + d3)/2
    d14 = (d1 + d4)/2
    d34 = (d3 + d4)/2
    smax = (d13**2 + d34**2 - d14**2)/(2*d23*d34)
    s = np.random.uniform(0, smax)
    #s = smax/2 # for simplicity, use the middle value
    d13 = d13*np.sqrt(1.0 - (s*d23/d13)**2)
    d23 = d23*np.sqrt(1.0 - s**2)
    cos0 = (d13**2 + d23**2 - 1)/(2*d13*d23)
    alpha = np.arccos(cos0)
    theta_c = np.pi/2 - alpha
    l_r = (1 + d23**2 - d13**2)/(2*d23)
    if th <= theta_c:
        x_min = d13/np.sin(th) - d23
    else:
        x_min = np.sqrt(1 - l_r**2)/np.tan(th) - l_r
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
    x_min = calculate_x_min(D1, D2, D3, th)
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
        x_min = calculate_x_min(D1, D2, D3, th)
        x_max = calculate_x_max(alpha, beta, th)
        phi, error = quad(rebound_angle_eq, x_min, x_max, args=(alpha, beta, th))
        phi /= (x_max - x_min)
        e, error = quad(rebound_res_eq, x_min, x_max, args=(alpha, beta, th))
        e /= (x_max - x_min)
    return phi, e

def calculate_rebound_3D(d1, d2, d3, d4, th, epsilon, nu):
    # normalize diameters
    D = (d1 + d2)/2
    D1 = d1/D
    D2 = d2/D
    D3 = d3/D
    D4 = d4/D
    # restitution coefficient
    mu = epsilon*D1**3/(D1**3 + epsilon*D2**3)
    alpha = (1 + epsilon)/(1 + mu) - 1
    beta = 1 - (2/7)*(1 - nu)/(1 + mu)
    x_min = calculate_x_min_3D(D1, D2, D3, D4, th)
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
        D4 = d4/D
        # restitution coefficient
        mu = epsilon*D1**3/(D1**3 + epsilon*D2**3)
        alpha = (1 + epsilon)/(1 + mu) - 1
        beta = 1 - (2/7)*(1 - nu)/(1 + mu)
        x_min = calculate_x_min_3D(D1, D2, D3, D4, th)
        x_max = calculate_x_max(alpha, beta, th)
        phi, error = quad(rebound_angle_eq, x_min, x_max, args=(alpha, beta, th))
        phi /= (x_max - x_min)
        e, error = quad(rebound_res_eq, x_min, x_max, args=(alpha, beta, th))
        e /= (x_max - x_min)
    return phi, e