import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import truncnorm
from scipy.stats import norm
from tqdm import tqdm
from scipy.special import erfc

# 根据对数正态分布的均值和标准差求正态分布的参数
def get_normal_params(log_mean, log_std):
    # 解关于sigma的方程
    def equation(sigma):
        mu = np.log(log_mean) - 0.5 * sigma ** 2
        return log_std ** 2 - np.exp(2 * mu + sigma ** 2) * (np.exp(sigma ** 2) - 1)
    # 简单的二分法求根
    left, right = 0.001, 10
    while right - left > 1e-6:
        mid = (left + right) / 2
        if equation(mid) > 0:
            left = mid
        else:
            right = mid
    sigma = mid
    mu = np.log(log_mean) - 0.5 * sigma ** 2
    return mu, sigma

def E_out(E_in, E_min_mean, e):
    lamb = 2*np.log((1 - e**2)*E_in/E_min_mean)
    sigma = np.sqrt(lamb)*np.log(2)
    mu = np.log((1 - e**2)*E_in) - lamb*np.log(2)
    E_out_mean = np.exp(mu + sigma**2/2)
    E_out_sigma = np.exp(2*mu + sigma**2)*(np.exp(sigma**2) - 1)
    return E_out_mean, E_out_sigma

def generate_truncated_lognormal(mu, sigma, dmin, dmax, size):
    a = (np.log(dmin) - mu) / sigma
    b = (np.log(dmax) - mu) / sigma
    samples = truncnorm.rvs(a, b, loc=mu, scale=sigma, size=size)
    return np.exp(samples)

def generate_bimodal(mu1, sigma1, mu2, sigma2, weight1, dmin, dmax, size):
    """生成截断的双峰分布样本（混合两个正态分布）"""
    size1 = int(size * weight1)
    size2 = size - size1
    # 生成第一个截断正态分布
    a1 = (dmin - mu1) / sigma1
    b1 = (dmax - mu1) / sigma1
    samples1 = truncnorm.rvs(a1, b1, loc=mu1, scale=sigma1, size=size1)
    # 生成第二个截断正态分布
    a2 = (dmin - mu2) / sigma2
    b2 = (dmax - mu2) / sigma2
    samples2 = truncnorm.rvs(a2, b2, loc=mu2, scale=sigma2, size=size2)
    # 合并并打乱
    combined = np.concatenate([samples1, samples2])
    np.random.shuffle(combined)
    return combined

def calculate_x_min(d1, d2, d3, theta1):
    d13 = (d1 + d3)/2
    d23 = (d2 + d3)/2
    cos1 = (1 + d13**2 - d23**2)/(2*d13)
    cos2 = (1 + d23**2 - d13**2)/(2*d23)
    gamma = np.arccos(cos1)
    delta = np.arccos(cos2)
    epsilon = np.pi - gamma - delta
    theta_c = gamma + delta - np.pi/2
    l_r = (1 + d23**2 - d13**2)/(2*d23)
    if theta1 <= theta_c:
        x_min = d13/np.sin(theta1) - d23
    else:
        x_min = np.sqrt(1 - l_r**2)/np.tan(theta1) - l_r
    return x_min, delta, epsilon

def calculate_x_max(alpha, beta, theta1):
    x_max_try = 1/np.sin(theta1)
    x_min_try = 1/np.tan(theta1)
    attempt_num = 100
    delta_x = (x_max_try - x_min_try)/attempt_num
    xx = x_min_try
    for n in range(1, attempt_num + 1):
        x = x_min_try + n*delta_x
        x_sin_square = x**2*np.sin(theta1)**2
        x_cos = x*np.cos(theta1)
        neq_LHS = 1/(1 + beta/alpha)
        neq_RHS = x_sin_square - x_cos*np.sqrt(1 - x_sin_square)
        if neq_LHS < neq_RHS:
            break
        xx = x
    x_max = xx
    return x_max

def calculate_x_max_shallow(theta1):
    x_max = 1.0/theta1
    return x_max

def calculate_e(alpha, beta, x, theta1):
    # 任意x位置上的e_vx可能出现小于0的情况，但其在xmin和xmax范围内的积分是正的
    e_vx = -alpha*np.cos(theta1) + (alpha + beta)*x**2*np.sin(theta1)**2*np.cos(theta1) + (alpha + beta)*x*np.sin(theta1)**2*np.sqrt(1 - x**2*np.sin(theta1)**2)
    e_vz = alpha*np.sin(theta1) - (alpha + beta)*x**2*np.sin(theta1)**3 + (alpha + beta)*x*np.sin(theta1)*np.cos(theta1)*np.sqrt(1 - x**2*np.sin(theta1)**2)
    e = np.sqrt(e_vx**2 + e_vz**2)
    return e, e_vx, e_vz

def calculate_e_bar(alpha, beta, x_min, x_max, theta1):
    xmin_sin_square = x_min**2*np.sin(theta1)**2
    xmax_sin_square = x_max**2*np.sin(theta1)**2
    xmin_cos = x_min*np.cos(theta1)
    xmax_cos = x_max*np.cos(theta1)
    eqx_xmin = -(alpha + beta)*(1 - xmin_sin_square)**(3/2) + xmin_cos*(-3*alpha + (alpha + beta)*xmin_sin_square)
    eqx_xmax = -(alpha + beta)*(1 - xmax_sin_square)**(3/2) + xmax_cos*(-3*alpha + (alpha + beta)*xmax_sin_square)
    e_vx = (eqx_xmax - eqx_xmin)/(x_max - x_min)/3
    xmin_3_sin_2 = x_min**3*np.sin(theta1)**2
    xmax_3_sin_2 = x_max**3*np.sin(theta1)**2
    csc_theta1 = 1/np.sin(theta1)
    cot_theta1 = 1/np.tan(theta1)
    eqz_xmin = alpha*x_min - 1/3*(alpha + beta)*(xmin_3_sin_2 + cot_theta1*csc_theta1*(1 - xmin_sin_square)**(3/2))
    eqz_xmax = alpha*x_max - 1/3*(alpha + beta)*(xmax_3_sin_2 + cot_theta1*csc_theta1*(1 - xmax_sin_square)**(3/2))
    e_vz = (eqz_xmax - eqz_xmin)/(x_max - x_min)*np.sin(theta1)
    e = np.sqrt(e_vx**2 + e_vz**2)
    return e, e_vx, e_vz

def calculate_e_bar_shallow_comp(alpha, beta, d1, d2, d3, theta1, x_min, x_max):
    d13 = (d1 + d3)/2
    d23 = (d2 + d3)/2
    if d13 > 1.0:
        d13 = 1.0
    #a = d13/theta1 - d23
    #b = 1.0/theta1
    a = d13 - d23*theta1
    #e_vx = -alpha + (alpha + beta)*theta1**2/3.0*(a**2 + a*b + b**2) + (alpha + beta)/3.0/(b - a)*(1.0 - a**2*theta1**2)**(3/2)
    #e_vz = alpha*theta1 - (alpha + beta)*theta1**3/3.0*(a**2 + a*b + b**2) + (alpha + beta)/3.0/theta1/(b - a)*(1.0 - a**2*theta1**2)**(3/2)
    #e_vx = (-3.0*alpha + 3.0*d13*alpha - 3.0*d23*alpha*theta1 + (alpha + beta)*(1.0 - (d13 - d23*theta1)**3 + theta1*(1.0 - (d13 - d23*theta1)**2)**(3/2)))/3.0/(1.0 - d13 + d23*theta1)
    e_vx = (-3.0*alpha + 3.0*alpha*a + (alpha + beta)*(1.0 - a**3 + theta1*(1.0 - a**2)**(3/2)))/3.0/(1.0 - a)
    e_vz = theta1*(alpha - alpha*a + 1.0/3.0*(alpha + beta)*(-1.0 + a**3 + (1.0 - a**2)**(3/2)/theta1))/(1.0 - a)
    e = np.sqrt(e_vx**2 + e_vz**2)
    return e, e_vx, e_vz

def calculate_e_bar_shallow_simp(alpha, beta, d1, d2, d3, theta1):
    d13 = (d1 + d3)/2
    d23 = (d2 + d3)/2
    if d2 == d3:
        e_vx = beta - d2*(alpha + beta)*theta1
        e_vz = 2.0*np.sqrt(2.0)/3.0*np.sqrt(d2)*(alpha + beta)*np.sqrt(theta1) - beta*theta1
        e = np.sqrt(beta**2) + (-2.0*d2*beta*(alpha + beta) + 8.0/9.0*d2*(alpha + beta)**2)*theta1/2.0/np.sqrt(beta**2)
    else:
        d13 = (d1 + d3)/2
        d23 = (d2 + d3)/2
        if d13 > 1.0:
            d13 = 1.0
        e_vx0 = 1.0/3.0*(-2.0*alpha + d13*alpha + d13**2*alpha + beta + d13*beta + d13**2*beta)
        e_vx1 = 1.0/3.0*(np.sqrt(1.0 - d13**2) + d13*(np.sqrt(1.0 - d13**2) - 2.0*d23) - d23)*(alpha + beta)*theta1
        e_vx = e_vx0 + e_vx1
        e_vz = 1.0/3.0*(1.0 + d13)*np.sqrt(1.0 - d13**2)*(alpha + beta)
        e = np.sqrt(e_vx**2 + e_vz**2)
    return e, e_vx, e_vz

def calculate_survive(d1, d2, d3, psi1, psi2, g, v2_x, v2_z, e):
    if v2_x > 0.0:
        zb = (1.0 - np.sin(psi1))*0.5*(d1 + d2)
        z_final = v2_z**2/(2.0*g)
        if z_final > zb:
            res_x = (np.sqrt((v2_z/g)**2 - (2.0*zb)/g) + v2_z/g)*v2_x - 0.5*(d1 + d2)*np.cos(psi1)
        else:
            res_x = -1
    elif v2_x < 0.0:
        zb = (1.0 - np.sin(psi2))*0.5*(d1 + d3)
        z_final = v2_z**2/(2.0*g)
        if z_final > zb:
            res_x = (np.sqrt((v2_z/g)**2 - (2.0*zb)/g) + v2_z/g)*v2_x - 0.5*(d1 + d3)*np.cos(psi2)
        else:
            res_x = -1
    else:
        res_x = -1
    if res_x <= 0.0 or e <= 0.0:
        is_survive = False
    else:
        is_survive = True
    return is_survive

def calculate_survive1(d1, d2, v1, theta1, x, g, v2_x, v2_z, e):
    v1_vec = np.array([v1*np.cos(theta1), v1*np.sin(theta1)])
    x_vec = np.array([x, 0])
    t = x_vec @ v1_vec / (v1**2) + (1/v1**2)*np.sqrt((1-x**2)*v1**2 + (x_vec @ v1_vec)**2)
    v1t = v1 * t
    zb = (1 - v1t*np.sin(theta1))*0.5*(d1 + d2)
    z_final = v2_z**2/(2.0*g)
    if z_final > zb:
        x_distance = (v1t*np.cos(theta1) - x)*0.5*(d1 + d2)
        x_min_try = 1/np.tan(theta1)*0.5*(d1 + d2)
        if x_distance < x_min_try:
            res_x = (np.sqrt((v2_z/g)**2 - (2.0*zb)/g) + v2_z/g)*v2_x - x_distance
        else:
            res_x = 1
    else:
        res_x = -1
    if res_x <= 0.0 or e <= 0.0:
        is_survive = False
    else:
        is_survive = True
    return is_survive

#-------------------------------------------------------------------
output_rebound_ratio = False # output rebound ratio
output_e = True # output e
output_ex = False # output e_x
output_ez = False # output e_z
output_theta2 = True # output rebound angle
output_vn = True # output eject velocity
output_Nej = True # output eject number

distribution = 1 # 0:uniform, 1:lognormal, 2:bidisperse, 3:polydisperse, 4:normal
shallow = False # shallow impact
simplify = False # first order approximation
lognormal_param = True # lognormal distribution parameters

d_min = 0.5e-4
d_max = 6e-4
normal_E = 2e-4
normal_D = 4e-4
mu = -8.30271
sigma = 0.25778
mu1 = (d_min + d_max) * 0.3  # 第一个峰靠左
mu2 = (d_min + d_max) * 0.7  # 第二个峰靠右
sigma1 = (d_max - d_min) * 0.1
sigma2 = (d_max - d_min) * 0.1
weight1 = 0.5 # 第一个峰的权重
gamma_ej = 0.05 # the fraction of remaining energy to eject particles
rho = 2650
g = 9.8*(1 - 1.263/rho)
epsilon = 0.15
nu = -0.9
num_samples = 100
#------------------------------------------------------------------
# average bed diameter
if distribution == 0:
    d2_array = np.random.uniform(d_min, d_max, size=num_samples)
    d2_mid = np.percentile(d2_array, 50)
elif distribution == 1:
    if not lognormal_param:
        mu, sigma = get_normal_params(normal_E, normal_D)
    d2_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, num_samples)
    d2_mid = np.percentile(d2_array, 50)
    d2_90 = np.percentile(d2_array, 90)
elif distribution == 2:
    d2_mid = (d_min + d_max) * 0.5
elif distribution == 3:
    d2_array = generate_bimodal(mu1, sigma1, mu2, sigma2, weight1, d_min, d_max, num_samples)
    d2_mid = np.percentile(d2_array, 50)
elif distribution == 4:
    d2_array = np.random.normal(loc=normal_E, scale=normal_D, size=num_samples)
    d2_mid = np.percentile(d2_array, 50)
# impactor
coarse_d1 = 1.4*d2_mid
medium_d1 = d2_mid
fine_d1 = 0.73*d2_mid

"""other's data"""
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
    'd_in': [coarse_d1, coarse_d1, coarse_d1, coarse_d1],
    'v_in': [2.7323, 2.7013, 2.7070, 2.7979],
    'ang_in': [13.94, 14.75, 14.73, 15.04],
    'ang_re': [20.95, 23.03, 22.55, 25.63],
    'e': [0.58, 0.56, 0.56, 0.58],
    'Nej': [5.7, 5.09, 5.12, 4.64],
    'v_ej': [0.243, 0.2422, 0.2522, 0.2637]
}
Rice95_v_many_medium= {
    'd_in': [medium_d1, medium_d1, medium_d1, medium_d1],
    'v_in': [3.1649, 3.3638, 3.3604, 3.2963],
    'ang_in': [11.82, 11.47, 11.53, 11.36],
    'ang_re': [29.95, 30.31, 31.06, 29.56],
    'e': [0.57, 0.54, 0.59, 0.58],
    'Nej': [2.68, 3.43, 2.67, 2.76],
    'v_ej': [0.2544, 0.252, 0.2607, 0.295]
}
Rice95_v_many_fine= {
    'd_in': [fine_d1, fine_d1, fine_d1],
    'v_in': [3.8493, 3.7801, 3.7411],
    'ang_in': [10.85, 10.24, 10.46],
    'ang_re': [44.63, 38.03, 37.85],
    'e': [0.52, 0.56, 0.58],
    'Nej': [1.86, 1.71, 1.58],
    'v_ej': [0.2757, 0.2566, 0.2773]
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

d1_array = Rice95_v_many_coarse['d_in'] + Rice95_v_many_medium['d_in'] + Rice95_v_many_fine['d_in']
d1_array = d1_array + Willetts89_v_many_coarse['d_in'] + Willetts89_v_many_medium['d_in'] + Willetts89_v_many_fine['d_in']
v1_array = Rice95_v_many_coarse['v_in'] + Rice95_v_many_medium['v_in'] + Rice95_v_many_fine['v_in']
v1_array = v1_array + Willetts89_v_many_coarse['v_in'] + Willetts89_v_many_medium['v_in'] + Willetts89_v_many_fine['v_in']
theta1_array = Rice95_v_many_coarse['ang_in'] + Rice95_v_many_medium['ang_in'] + Rice95_v_many_fine['ang_in']
theta1_array = theta1_array + Willetts89_v_many_coarse['ang_in'] + Willetts89_v_many_medium['ang_in'] + Willetts89_v_many_fine['ang_in']
iteration_num = len(d1_array)

rebound_ratio0_list_0 = []
rebound_ratio0_list_1 = []
rebound_ratio0_list_2 = []
e0_bar_list_0 = []
e0_bar_list_1 = []
e0_bar_list_2 = []
ez0_bar_list_0 = []
ez0_bar_list_1 = []
ez0_bar_list_2 = []
ex0_bar_list_0 = []
ex0_bar_list_1 = []
ex0_bar_list_2 = []
theta20_bar_list_0 = []
theta20_bar_list_1 = []
theta20_bar_list_2 = []
Nej_bar_list_0 = []
Nej_bar_list_1 = []
Nej_bar_list_2 = []
vn_bar_list_0 = []
vn_bar_list_1 = []
vn_bar_list_2 = []
for n in tqdm(range(iteration_num)):
    d1 = d1_array[n]
    v1 = v1_array[n]
    theta1 = np.radians(theta1_array[n])
    e0_list_0 = []
    e0_list_1 = []
    e0_list_2 = []
    ex0_list_0 = []
    ex0_list_1 = []
    ex0_list_2 = []
    ez0_list_0 = []
    ez0_list_1 = []
    ez0_list_2 = []
    theta20_list_0 = []
    theta20_list_1 = []
    theta20_list_2 = []
    Nej_list_0 = []
    Nej_list_1 = []
    Nej_list_2 = []
    vn_list_0 = []
    vn_list_1 = []
    vn_list_2 = []
    current_n = 0
    while current_n < num_samples:
        current_n += 1
        if distribution == 0:
            d_array = np.random.uniform(d_min, d_max, size=3)
            d2 = d_array[1]
            d3 = d_array[2]
        elif distribution == 1:
            if not lognormal_param:
                mu, sigma = get_normal_params(normal_E, normal_D)
            d_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, 3)
            d2 = d_array[1]
            d3 = d_array[2]
        elif distribution == 2:
            r = np.random.uniform(0.5,1.5,size=3)
            rr = np.floor(r)
            d_bin = d_max - d_min
            d2 = d_min + rr[1]*d_bin
            d3 = d_min + rr[2]*d_bin
        elif distribution == 3:
            d3 = np.random.choice(
                generate_bimodal(mu1, sigma1, mu2, sigma2, weight1, d_min, d_max, 1000),
                size=1, replace=False
            )
            d2_min = (-(d1+d3)+np.sqrt((d1+d3)**2+4*d1*d3))/2
            d2_min = max(d_min, d2_min)
            if d2_min > d_max:
                continue
            d2_candidates = generate_bimodal(mu1, sigma1, mu2, sigma2, weight1, d2_min, d_max, 100)
            if len(d2_candidates) == 0:
                continue
            d2 = np.random.choice(d2_candidates)
        elif distribution == 4:
            d_array = np.random.normal(loc=normal_E, scale=normal_D, size=3)
            d2 = d_array[1]
            d3 = d_array[2]
        # 三种情况:
        # d1, d2, d3分别为不同的直径
        # d1, d2=d3
        # d1=d2=d3
        for i in range(3):
            if i == 1:
                d3 = d2
            elif i == 2:
                d2 = d1
                d3 = d1

            d = (d1 + d2)/2
            d1_hat = d1/d
            d2_hat = d2/d
            d3_hat = d3/d
            # restitution coefficient
            mu_re = epsilon*d1_hat**3/(d1_hat**3 + epsilon*d2_hat**3)
            alpha = (1 + epsilon)/(1 + mu_re) - 1
            beta = 1 - (2/7)*(1 - nu)/(1 + mu_re)
            x_min, psi1, psi2 = calculate_x_min(d1_hat, d2_hat, d3_hat, theta1)

            x_max = calculate_x_max(alpha, beta, theta1)
            if x_min < x_max:
                dont_change_flag = True
            else:
                dont_change_flag = False

            if dont_change_flag:
                x0_hat = np.random.uniform(x_min, x_max)
                e0, evx0, evz0 = calculate_e(alpha, beta, x0_hat, theta1)
            else:
                d_temp = d2
                d2 = d3
                d3 = d_temp
                d = (d1 + d2)/2
                d1_hat = d1/d
                d2_hat = d2/d
                d3_hat = d3/d
                mu_re = epsilon*d1_hat**3/(d1_hat**3 + epsilon*d2_hat**3)
                alpha = (1 + epsilon)/(1 + mu_re) - 1
                beta = 1 - (2/7)*(1 - nu)/(1 + mu_re)
                x_min, psi1, psi2 = calculate_x_min(d1_hat, d2_hat, d3_hat, theta1)
                x_max = calculate_x_max(alpha, beta, theta1)
                x0_hat = np.random.uniform(x_min, x_max)
                e0, evx0, evz0 = calculate_e(alpha, beta, x0_hat, theta1)

            if not lognormal_param:
                mu, sigma = get_normal_params(normal_E, normal_D)
            dn = generate_truncated_lognormal(mu, sigma, d_min, d_max, 1)[0]
            m1 = rho*(d1**3*np.pi/6)
            mn = rho*(dn**3*np.pi/6)
            m_mid = rho*(d2_mid**3*np.pi/6)
            C_Ec = rho*g*np.pi/6*d2_90
            Ec_min = C_Ec*d_min**3
            Ec_max = C_Ec*d_max**3
            mu_Ec = np.log(C_Ec) + 3.0*mu
            sigma_Ec = 3.0*sigma
            z_min = (np.log(Ec_min) - mu_Ec - sigma_Ec**2) / sigma_Ec
            z_max = (np.log(Ec_max) - mu_Ec - sigma_Ec**2) / sigma_Ec
            P_Ec = norm.cdf((np.log(Ec_max) - mu_Ec) / sigma_Ec) - norm.cdf((np.log(Ec_min) - mu_Ec) / sigma_Ec)
            Ec = (C_Ec*np.exp(3.0*mu + 9.0*sigma**2/2.0) * (norm.cdf(z_max) - norm.cdf(z_min))) / P_Ec
            E1 = 0.5*m1*v1**2
            k_max = (1.0 - e0**2)*E1/Ec
            lambda_ej = 2.0*np.log(k_max)
            mu_ej = np.log((1.0 - e0**2)*E1) - lambda_ej*np.log(2.0)
            sigma_ej = np.sqrt(lambda_ej)*np.log(2.0)
            En_bar = np.exp(mu_ej + sigma_ej**2/2)
            #Ec_ej_min = mn*g*d_min
            #Ec_ej_max = mn*g*d_max
            #mu_Ec_ej = np.log(mn*g) + mu
            #sigma_Ec_ej = sigma
            #z_min = (np.log(Ec_ej_min) - mu_Ec_ej - sigma_Ec_ej**2) / sigma_Ec_ej
            #z_max = (np.log(Ec_ej_max) - mu_Ec_ej - sigma_Ec_ej**2) / sigma_Ec_ej
            #P_Ec_ej = norm.cdf((np.log(Ec_ej_max) - mu_Ec_ej) / sigma_Ec_ej) - norm.cdf((np.log(Ec_ej_min) - mu_Ec_ej) / sigma_Ec_ej)
            #Ec_ej = (mn*g*np.exp(mu + sigma**2/2.0) * (norm.cdf(z_max) - norm.cdf(z_min))) / P_Ec_ej
            Ec_ej = Ec
            Nej = gamma_ej*(1.0 - e0**2)*E1/En_bar/2.0*erfc((np.log(Ec_ej) - mu_ej)/(np.sqrt(2.0)*sigma_ej))
            Nej = max(Nej, 0.0)
            Nejed = 0
            vn_list = []
            while Nejed <= Nej:
                Nejed += 1
                dn = generate_truncated_lognormal(mu, sigma, d_min, d_max, 1)[0]
                mn = rho*(dn**3*np.pi/6)
                En = 0.0
                while En < Ec_ej:
                    En = np.random.lognormal(mu_ej, sigma_ej, size=1)[0]
                vn = np.sqrt(2.0*En/mn)
                vn_list.append(vn)
            vn = np.mean(vn_list)
            #C_ej = np.sqrt(2.0/mn)
            #Qn = (np.log(Ec_ej) - mu_ej - sigma_ej**2/2.0)/(np.sqrt(2.0)*sigma_ej)
            #Pn = (np.log(Ec_ej) - mu_ej)/(np.sqrt(2.0)*sigma_ej)
            #vn = erfc(Qn)/erfc(Pn)*C_ej*np.exp(mu_ej/2.0 + sigma_ej**2/8.0)
            if i == 0:
                vn_list_0.append(vn)
                Nej_list_0.append(Nej)
            elif i == 1:
                vn_list_1.append(vn)
                Nej_list_1.append(Nej)
            elif i == 2:
                vn_list_2.append(vn)
                Nej_list_2.append(Nej)

            ez0 = evz0 #/np.sin(theta1)
            ex0 = evx0 #/np.cos(theta1)
            theta20 = np.arctan(evz0/evx0)
            if theta20 < 0:
                theta20 += np.pi
            if i == 0:
                e0_list_0.append(e0)
                ez0_list_0.append(ez0)
                ex0_list_0.append(ex0)
                theta20_list_0.append(theta20)
            elif i == 1:
                e0_list_1.append(e0)
                ez0_list_1.append(ez0)
                ex0_list_1.append(ex0)
                theta20_list_1.append(theta20)
            elif i == 2:
                e0_list_2.append(e0)
                ez0_list_2.append(ez0)
                ex0_list_2.append(ex0)
                theta20_list_2.append(theta20)
            v20 = e0*v1
            v2_z0 = v1*evz0
            v2_x0 = v1*evx0
    """d1, d2, d3"""
    ez0_bar_0 = np.mean(ez0_list_0)
    ex0_bar_0 = np.mean(ex0_list_0)
    e0_bar_0 = np.sqrt(ex0_bar_0**2 + ez0_bar_0**2)
    theta20_bar_0 = np.mean(theta20_list_0)
    e0_bar_list_0.append(e0_bar_0)
    ez0_bar_list_0.append(ez0_bar_0)
    ex0_bar_list_0.append(ex0_bar_0)
    theta20_bar_list_0.append(theta20_bar_0)
    """d2=d3"""
    ez0_bar_1 = np.mean(ez0_list_1)
    ex0_bar_1 = np.mean(ex0_list_1)
    e0_bar_1 = np.sqrt(ex0_bar_1**2 + ez0_bar_1**2)
    theta20_bar_1 = np.mean(theta20_list_1)
    e0_bar_list_1.append(e0_bar_1)
    ez0_bar_list_1.append(ez0_bar_1)
    ex0_bar_list_1.append(ex0_bar_1)
    theta20_bar_list_1.append(theta20_bar_1)
    """d1=d2=d3"""
    ez0_bar_2 = np.mean(ez0_list_2)
    ex0_bar_2 = np.mean(ex0_list_2)
    e0_bar_2 = np.sqrt(ex0_bar_2**2 + ez0_bar_2**2)
    theta20_bar_2 = np.mean(theta20_list_2)
    e0_bar_list_2.append(e0_bar_2)
    ez0_bar_list_2.append(ez0_bar_2)
    ex0_bar_list_2.append(ex0_bar_2)
    theta20_bar_list_2.append(theta20_bar_2)
    """ejection velocity"""
    vn_bar_0 = np.mean(vn_list_0)
    vn_bar_1 = np.mean(vn_list_1)
    vn_bar_2 = np.mean(vn_list_2)
    vn_bar_list_0.append(vn_bar_0)
    vn_bar_list_1.append(vn_bar_1)
    vn_bar_list_2.append(vn_bar_2)
    """ejection number"""
    Nej_bar_0 = np.mean(Nej_list_0)
    Nej_bar_1 = np.mean(Nej_list_1)
    Nej_bar_2 = np.mean(Nej_list_2)
    Nej_bar_list_0.append(Nej_bar_0)
    Nej_bar_list_1.append(Nej_bar_1)
    Nej_bar_list_2.append(Nej_bar_2)

len_R95_c = len(Rice95_v_many_coarse['e'])
len_R95_m = len(Rice95_v_many_medium['e'])
len_R95_f = len(Rice95_v_many_fine['e'])
len_W89_c = len(Willetts89_v_many_coarse['e'])
len_W89_m = len(Willetts89_v_many_medium['e'])
len_W89_f = len(Willetts89_v_many_fine['e'])
color_map = 'jet'
size_c = 80
size_m = 50
size_f = 20
x_array_total = []
y_array_total = []
if output_e:
    plt.figure(1, figsize=(8, 6))
    colors = Rice95_v_many_coarse['ang_in']
    x_array = Rice95_v_many_coarse['e']
    y_array = e0_bar_list_0[0:len_R95_c]
    len0 = len_R95_c
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='o', label='Rice95 coarse')
    min_xy = min(x_array + y_array)
    max_xy = max(x_array + y_array)
    x_array_total += x_array
    y_array_total += y_array

    colors = Rice95_v_many_medium['ang_in']
    x_array = Rice95_v_many_medium['e']
    y_array = e0_bar_list_0[len0:len0 + len_R95_m]
    len0 += len_R95_m
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='o', label='Rice95 medium')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    colors = Rice95_v_many_fine['ang_in']
    x_array = Rice95_v_many_fine['e']
    y_array = e0_bar_list_0[len0:len0 + len_R95_f]
    len0 += len_R95_f
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='o', label='Rice95 fine')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    colors = Willetts89_v_many_coarse['ang_in']
    x_array = Willetts89_v_many_coarse['e']
    y_array = e0_bar_list_0[len0:len0 + len_W89_c]
    len0 += len_W89_c
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='^', label='Willetts89 coarse')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    colors = Willetts89_v_many_medium['ang_in']
    x_array = Willetts89_v_many_medium['e']
    y_array = e0_bar_list_0[len0:len0 + len_W89_m]
    len0 += len_W89_m
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='^', label='Willetts89 medium')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    colors = Willetts89_v_many_fine['ang_in']
    x_array = Willetts89_v_many_fine['e']
    y_array = e0_bar_list_0[len0:len0 + len_W89_f]
    len0 += len_W89_f
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='^', label='Willetts89 fine')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    line_points = [min_xy, max_xy]
    plt.plot(line_points, line_points, 'k--', label='y=x')

    plt.xlabel('$e_{exp}$')
    plt.ylabel('$e_{sim}$')
    #plt.xlim(min_xy, max_xy)
    #plt.ylim(min_xy, max_xy)
    plt.colorbar(label='Impact angle (degrees)')

    ms_c = np.sqrt(size_c)
    ms_m = np.sqrt(size_m)
    ms_f = np.sqrt(size_f)
    circle_c = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_c, label='Rice95 coarse')
    circle_m = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_m, label='Rice95 medium')
    circle_f = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_f, label='Rice95 fine')
    triangle_c = mlines.Line2D([], [], color='k', marker='^', linestyle='None', markersize=ms_c, label='Willetts89 coarse')
    triangle_m = mlines.Line2D([], [], color='k', marker='^', linestyle='None', markersize=ms_m, label='Willetts89 medium')
    triangle_f = mlines.Line2D([], [], color='k', marker='^', linestyle='None', markersize=ms_f, label='Willetts89 fine')
    plt.legend(handles=[circle_c, circle_m, circle_f, triangle_c, triangle_m, triangle_f], loc='best')

    plt.figure(2, figsize=(8, 6))
    x_array = Rice95_v_many_coarse['e'] + Rice95_v_many_medium['e'] + Rice95_v_many_fine['e']
    x_array = x_array + Willetts89_v_many_coarse['e'] + Willetts89_v_many_medium['e'] + Willetts89_v_many_fine['e']
    y_array1 = e0_bar_list_1
    y_array2 = e0_bar_list_2
    min_xy = min(x_array + y_array1 + y_array2)
    max_xy = max(x_array + y_array1 + y_array2)
    line_points = [min_xy, max_xy]
    plt.plot(line_points, line_points, 'k--', label='y=x')
    plt.scatter(x_array, y_array1, c='k', marker='s', label='d2=d3')
    plt.scatter(x_array, y_array2, c='k', marker='*', label='d1=d2=d3')
    plt.xlabel('$e_{exp}$')
    plt.ylabel('$e_{sim}$')
    plt.legend()

if output_theta2:
    plt.figure(3, figsize=(8, 6))
    colors = Rice95_v_many_coarse['ang_in']
    x_array = Rice95_v_many_coarse['ang_re']
    y_array = np.degrees(theta20_bar_list_0[0:len_R95_c])
    y_array = y_array.tolist()
    len0 = len_R95_c
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='o', label='Rice95 coarse')
    min_xy = min(x_array + y_array)
    max_xy = max(x_array + y_array)
    x_array_total += x_array
    y_array_total += y_array

    colors = Rice95_v_many_medium['ang_in']
    x_array = Rice95_v_many_medium['ang_re']
    y_array = np.degrees(theta20_bar_list_0[len0:len0 + len_R95_m])
    y_array = y_array.tolist()
    len0 += len_R95_m
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='o', label='Rice95 medium')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    colors = Rice95_v_many_fine['ang_in']
    x_array = Rice95_v_many_fine['ang_re']
    y_array = np.degrees(theta20_bar_list_0[len0:len0 + len_R95_f])
    y_array = y_array.tolist()
    len0 += len_R95_f
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='o', label='Rice95 fine')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    colors = Willetts89_v_many_coarse['ang_in']
    x_array = Willetts89_v_many_coarse['ang_re']
    y_array = np.degrees(theta20_bar_list_0[len0:len0 + len_W89_c])
    y_array = y_array.tolist()
    len0 += len_W89_c
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='^', label='Willetts89 coarse')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    colors = Willetts89_v_many_medium['ang_in']
    x_array = Willetts89_v_many_medium['ang_re']
    y_array = np.degrees(theta20_bar_list_0[len0:len0 + len_W89_m])
    y_array = y_array.tolist()
    len0 += len_W89_m
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='^', label='Willetts89 medium')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    colors = Willetts89_v_many_fine['ang_in']
    x_array = Willetts89_v_many_fine['ang_re']
    y_array = np.degrees(theta20_bar_list_0[len0:len0 + len_W89_f])
    y_array = y_array.tolist()
    len0 += len_W89_f
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='^', label='Willetts89 fine')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    line_points = [min_xy, max_xy]
    plt.plot(line_points, line_points, 'k--', label='y=x')

    plt.xlabel('$\\theta_{re,exp}$ (degrees)')
    plt.ylabel('$\\theta_{re,sim}$ (degrees)')
    #plt.xlim(min_xy, max_xy)
    #plt.ylim(min_xy, max_xy)
    plt.colorbar(label='Impact angle (degrees)')

    ms_c = np.sqrt(size_c)
    ms_m = np.sqrt(size_m)
    ms_f = np.sqrt(size_f)
    circle_c = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_c, label='Rice95 coarse')
    circle_m = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_m, label='Rice95 medium')
    circle_f = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_f, label='Rice95 fine')
    triangle_c = mlines.Line2D([], [], color='k', marker='^', linestyle='None', markersize=ms_c, label='Willetts89 coarse')
    triangle_m = mlines.Line2D([], [], color='k', marker='^', linestyle='None', markersize=ms_m, label='Willetts89 medium')
    triangle_f = mlines.Line2D([], [], color='k', marker='^', linestyle='None', markersize=ms_f, label='Willetts89 fine')
    plt.legend(handles=[circle_c, circle_m, circle_f, triangle_c, triangle_m, triangle_f], loc='best')

    plt.figure(4, figsize=(8, 6))
    x_array = Rice95_v_many_coarse['ang_re'] + Rice95_v_many_medium['ang_re'] + Rice95_v_many_fine['ang_re']
    x_array = x_array + Willetts89_v_many_coarse['ang_re'] + Willetts89_v_many_medium['ang_re'] + Willetts89_v_many_fine['ang_re']
    y_array1 = np.degrees(theta20_bar_list_1)
    y_array2 = np.degrees(theta20_bar_list_2)
    y_array1 = y_array1.tolist()
    y_array2 = y_array2.tolist()
    min_xy = min(x_array + y_array1 + y_array2)
    max_xy = max(x_array + y_array1 + y_array2)
    line_points = [min_xy, max_xy]
    plt.plot(line_points, line_points, 'k--', label='y=x')
    plt.scatter(x_array, y_array1, c='k', marker='s', label='d2=d3')
    plt.scatter(x_array, y_array2, c='k', marker='*', label='d1=d2=d3')
    plt.xlabel('$\\theta_{re,exp}$ (degrees)')
    plt.ylabel('$\\theta_{re,sim}$ (degrees)')
    plt.legend()

if output_Nej:
    plt.figure(5, figsize=(8, 6))
    colors = Rice95_v_many_coarse['ang_in']
    x_array = Rice95_v_many_coarse['Nej']
    y_array = Nej_bar_list_0[0:len_R95_c]
    len0 = len_R95_c
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='o', label='Rice95 coarse')
    min_xy = min(x_array + y_array)
    max_xy = max(x_array + y_array)
    x_array_total += x_array
    y_array_total += y_array

    colors = Rice95_v_many_medium['ang_in']
    x_array = Rice95_v_many_medium['Nej']
    y_array = Nej_bar_list_0[len0:len0 + len_R95_m]
    len0 += len_R95_m
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='o', label='Rice95 medium')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    colors = Rice95_v_many_fine['ang_in']
    x_array = Rice95_v_many_fine['Nej']
    y_array = Nej_bar_list_0[len0:len0 + len_R95_f]
    len0 += len_R95_f
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='o', label='Rice95 fine')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    line_points = [min_xy, max_xy]
    plt.plot(line_points, line_points, 'k--', label='y=x')

    plt.xlabel('$N_{ej,exp}$')
    plt.ylabel('$N_{ej,sim}$')
    #plt.xlim(min_xy, max_xy)
    #plt.ylim(min_xy, max_xy)
    plt.colorbar(label='Impact angle (degrees)')
    ms_c = np.sqrt(size_c)
    ms_m = np.sqrt(size_m)
    ms_f = np.sqrt(size_f)
    circle_c = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_c, label='Rice95 coarse')
    circle_m = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_m, label='Rice95 medium')
    circle_f = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_f, label='Rice95 fine')
    plt.legend(handles=[circle_c, circle_m, circle_f], loc='best')

if output_vn:
    plt.figure(6, figsize=(8, 6))
    colors = Rice95_v_many_coarse['ang_in']
    x_array = Rice95_v_many_coarse['v_ej']
    y_array = vn_bar_list_0[0:len_R95_c]
    len0 = len_R95_c
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='o', label='Rice95 coarse')
    min_xy = min(x_array + y_array)
    max_xy = max(x_array + y_array)
    x_array_total += x_array
    y_array_total += y_array

    colors = Rice95_v_many_medium['ang_in']
    x_array = Rice95_v_many_medium['v_ej']
    y_array = vn_bar_list_0[len0:len0 + len_R95_m]
    len0 += len_R95_m
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='o', label='Rice95 medium')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    colors = Rice95_v_many_fine['ang_in']
    x_array = Rice95_v_many_fine['v_ej']
    y_array = vn_bar_list_0[len0:len0 + len_R95_f]
    len0 += len_R95_f
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='o', label='Rice95 fine')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)
    x_array_total += x_array
    y_array_total += y_array

    #colors = Willetts89_v_many_coarse['ang_in']
    #x_array = Willetts89_v_many_coarse['v_ej']
    #y_array = vn_bar_list_0[len0:len0 + len_W89_c]
    #len0 += len_W89_c
    #plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='^', label='Willetts89 coarse')
    #min_xy = min(min(x_array + y_array), min_xy)
    #max_xy = max(max(x_array + y_array), max_xy)

    #colors = Willetts89_v_many_medium['ang_in']
    #x_array = Willetts89_v_many_medium['v_ej']
    #y_array = vn_bar_list_0[len0:len0 + len_W89_m]
    #len0 += len_W89_m
    #plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='^', label='Willetts89 medium')
    #min_xy = min(min(x_array + y_array), min_xy)
    #max_xy = max(max(x_array + y_array), max_xy)

    #colors = Willetts89_v_many_fine['ang_in']
    #x_array = Willetts89_v_many_fine['v_ej']
    #y_array = vn_bar_list_0[len0:len0 + len_W89_f]
    #len0 += len_W89_f
    #plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='^', label='Willetts89 fine')
    #min_xy = min(min(x_array + y_array), min_xy)
    #max_xy = max(max(x_array + y_array), max_xy)

    line_points = [min_xy, max_xy]
    plt.plot(line_points, line_points, 'k--', label='y=x')

    plt.xlabel('$v_{n,exp}$ (m/s)')
    plt.ylabel('$v_{n,sim}$ (m/s)')
    #plt.xlim(min_xy, max_xy)
    #plt.ylim(min_xy, max_xy)
    plt.colorbar(label='Impact angle (degrees)')

    ms_c = np.sqrt(size_c)
    ms_m = np.sqrt(size_m)
    ms_f = np.sqrt(size_f)
    circle_c = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_c, label='Rice95 coarse')
    circle_m = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_m, label='Rice95 medium')
    circle_f = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_f, label='Rice95 fine')
    plt.legend(handles=[circle_c, circle_m, circle_f], loc='best')

plt.show()

rmse = np.sqrt(np.mean((np.array(x_array_total) - np.array(y_array_total))**2))
print(f'RMSE: {rmse:.4f}')