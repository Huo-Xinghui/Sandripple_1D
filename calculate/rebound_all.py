import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import truncnorm
from tqdm import tqdm
from scipy.special import erfc
from scipy.integrate import quad
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.usetex'] = True
label_size = 20
ticks_size = 15

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
    cos0 = (d13**2 + d23**2 - 1)/(2*d13*d23)
    cos1 = (1 + d13**2 - d23**2)/(2*d13)
    cos2 = (1 + d23**2 - d13**2)/(2*d23)
    alpha = np.arccos(cos0)
    gamma = np.arccos(cos1)
    delta = np.arccos(cos2)
    epsilon = np.pi - gamma - delta
    #theta_c = gamma + delta - np.pi/2
    theta_c = np.pi/2 - alpha
    l_r = (1 + d23**2 - d13**2)/(2*d23)
    if theta1 <= theta_c:
        x_min = d13/np.sin(theta1) - d23
    else:
        x_min = np.sqrt(1 - l_r**2)/np.tan(theta1) - l_r
    return x_min, delta, epsilon


def calculate_x_min_3D(d1, d2, d3, d4, theta1):
    d13 = (d1 + d3)/2
    d23 = (d2 + d3)/2
    d14 = (d1 + d4)/2
    d34 = (d3 + d4)/2
    kmax = (d13**2 + d34**2 - d14**2)/(2*d13*d34)
    k = kmax/2 #np.random.uniform(0, kmax)
    d13 = d13*np.sqrt(1.0 - k**2)
    d23 = d23*np.sqrt(1.0 - (k*d13/d23)**2)
    cos0 = (d13**2 + d23**2 - 1)/(2*d13*d23)
    cos1 = (1 + d13**2 - d23**2)/(2*d13)
    cos2 = (1 + d23**2 - d13**2)/(2*d23)
    alpha = np.arccos(cos0)
    gamma = np.arccos(cos1)
    delta = np.arccos(cos2)
    epsilon = np.pi - gamma - delta
    #theta_c = gamma + delta - np.pi/2
    theta_c = np.pi/2 - alpha
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

def calculate_x_min_shallow(d1, d2, d3, theta1):
    d13 = (d1 + d3)/2
    d23 = (d2 + d3)/2
    cos1 = (1 + d13**2 - d23**2)/(2*d13)
    cos2 = (1 + d23**2 - d13**2)/(2*d23)
    gamma = np.arccos(cos1)
    delta = np.arccos(cos2)
    epsilon = np.pi - gamma - delta
    x_min = d13/np.sin(theta1) - d23
    return x_min, delta, epsilon

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

def rebound_angle_bar(x, alpha, beta, theta1):
    e_vx = -alpha*np.cos(theta1) + (alpha + beta)*x**2*np.sin(theta1)**2*np.cos(theta1) + (alpha + beta)*x*np.sin(theta1)**2*np.sqrt(1 - x**2*np.sin(theta1)**2)
    e_vz = alpha*np.sin(theta1) - (alpha + beta)*x**2*np.sin(theta1)**3 + (alpha + beta)*x*np.sin(theta1)*np.cos(theta1)*np.sqrt(1 - x**2*np.sin(theta1)**2)
    theta2 = np.arctan2(e_vz, e_vx)
    return theta2

def rebound_res_bar(x, alpha, beta, theta1):
    e_vx = -alpha*np.cos(theta1) + (alpha + beta)*x**2*np.sin(theta1)**2*np.cos(theta1) + (alpha + beta)*x*np.sin(theta1)**2*np.sqrt(1 - x**2*np.sin(theta1)**2)
    e_vz = alpha*np.sin(theta1) - (alpha + beta)*x**2*np.sin(theta1)**3 + (alpha + beta)*x*np.sin(theta1)*np.cos(theta1)*np.sqrt(1 - x**2*np.sin(theta1)**2)
    res = np.sqrt(e_vx**2 + e_vz**2)
    return res

def rebound_angle_shallow(x, alpha, beta, theta1):
    e_vx = -alpha + (alpha + beta)*(x**2*theta1**2+x*theta1**2*np.sqrt(1.0 -x**2*theta1**2))
    e_vz = alpha*theta1 - (alpha + beta)*(x**2*theta1**3 - x*theta1*np.sqrt(1.0 - x**2*theta1**2))
    theta2 = np.arctan2(e_vz, e_vx)
    return theta2

def rebound_res_shallow(x, alpha, beta, theta1):
    e_vx = -alpha + (alpha + beta)*(x**2*theta1**2+x*theta1**2*np.sqrt(1.0 -x**2*theta1**2))
    e_vz = alpha*theta1 - (alpha + beta)*(x**2*theta1**3 - x*theta1*np.sqrt(1.0 - x**2*theta1**2))
    res = np.sqrt(e_vx**2 + e_vz**2)
    return res

def rebound_angle_shallow_simp(alpha, beta, d1, d2, d3, theta1):
    d13 = (d1 + d3)/2
    d23 = (d2 + d3)/2
    K = 1.0 + alpha/beta
    L = d23 - (d13 - 1.0)/theta1
    S = np.sqrt(2.0*(d23*theta1 - d13 + 1.0))
    theta2 = 1.0/(theta1*L)*(K*S**3/3.0 - S**2*theta1/2.0)
    return theta2

def rebound_res_shallow_simp(alpha, beta, d1, d2, d3, theta1):
    d13 = (d1 + d3)/2
    d23 = (d2 + d3)/2
    e = beta - (beta**2 - alpha**2)/(2.0*beta)*(d23*theta1 - d13 + 1.0)
    return e

def theta2_simp(alpha, beta, d1, d2, d3, theta1, x):
    d13 = (d1 + d3)/2
    d23 = (d2 + d3)/2
    if not x : #<= d13/theta1:
        t = x + d23 - d13/theta1
        epsilon = (d23 - t)*theta1
        sq = 1.0 - d13**2 + 2.0*d13*epsilon
        N = (1.0 - d13**2)*alpha*theta1 - beta*d13**2*theta1 + (alpha + beta)*d13*np.sqrt(sq)
        D = (d13**2 - 1.0)*alpha + d13**2*beta
        tan_theta2 = N/D
    else:
    #    t = x - 1.0/theta1
    #    epsilon = -t*theta1
    #    sq = 2*epsilon
    #    N = -beta*theta1 + (alpha + beta)*np.sqrt(sq)
    #    D = beta
    #    tan_theta2 = N/D
    #    #tan_theta2 = -theta1 + (1.0 + alpha/beta)*np.sqrt(sq)
        t = x*theta1
        B = (alpha + beta)*t*np.sqrt(1 - t**2)
        A = alpha - (alpha + beta)*t**2
        C = -A
        D = -B
        #tan_theta2 = B/C + (B**2 - C**2)/C**2*theta1
        tan_theta2 = (B+A*theta1)/(C + D*theta1)
    #tan_theta2 = np.clip(tan_theta2, -np.pi, np.pi)  # 限制在[-pi, pi]范围内
    #if abs(tan_theta2) > np.pi:
    #    theta2 = 0.0
    #theta2 = tan_theta2
    theta2 = np.arctan(tan_theta2)
    return theta2

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
output_theta2_dist = False # output rebound angle distribution
output_vn = False # output eject velocity
output_Nej = False # output eject number


distribution = 1 # 0:uniform, 1:lognormal, 2:bidisperse, 3:polydisperse, 4:normal
variation_param = 1 # 0: v1, 1: theta1
shallow = False # shallow impact
simplify = False # first order approximation
lognormal_param = True # lognormal distribution parameters
Three_D = True # 3D bed
impactor = 3 #impactor type: 0: bed, 1: coarse, 2: medium, 3: fine

d_min = 1.5e-4
d_max = 6e-4
d_min_impactor_c = 3.55e-4
d_max_impactor_c = d_max
d_min_impactor_m = 2.5e-4
d_max_impactor_m = d_min_impactor_c
d_min_impactor_f = d_min
d_max_impactor_f = d_min_impactor_m
normal_E = 2e-4
normal_D = 4e-4
mu = -8.30271
sigma = 0.25778
mu1 = (d_min + d_max) * 0.3  # 第一个峰靠左
mu2 = (d_min + d_max) * 0.7  # 第二个峰靠右
sigma1 = (d_max - d_min) * 0.1
sigma2 = (d_max - d_min) * 0.1
weight1 = 0.5 # 第一个峰的权重
gamma_ej = 0.049 # the fraction of remaining energy to eject particles
rho = 2650
g = 9.8*(1 - 1.263/rho)
#epsilon = 0.71
#nu = -0.88
#epsilon = 0.73
#nu = -0.90
epsilon = 0.6605
nu = -0.8310
epsilon_La = 0.7538
nu_La = -0.9131
v1_single = 3.38
theta1_single = np.pi/18
v1_start = 1
v1_end = 10
case_num = 45
theta1_start = 0.5*np.pi/180
if shallow:
    theta1_end = 30/180*np.pi
else:
    theta1_end = 89.5/180*np.pi
num_samples = 100
#------------------------------------------------------------------
# impactor diameter
if distribution == 0:
    d1_array = np.random.uniform(d_min, d_max, size=num_samples)
    d_mid = np.percentile(d1_array, 50)
    if impactor == 1:
        d1 = np.percentile(d1_array, 90)
    elif impactor == 2:
        d1 = d_mid
    elif impactor == 3:
        d1 = np.percentile(d1_array, 10)
elif distribution == 1:
    if not lognormal_param:
        mu, sigma = get_normal_params(normal_E, normal_D)
    d1_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, num_samples)
    d_mid = np.percentile(d1_array, 50)
    if impactor == 1:
        d1 = 1.4*d_mid
        #d1 = np.mean(d1_array[(d1_array >= 0.000425) & (d1_array <= 0.0006)])
        #d1 = np.percentile(d1_array, 90)
    elif impactor == 2:
        d1 = d_mid
        #d1 = np.mean(d1_array[(d1_array >= 0.0003) & (d1_array <= 0.000355)])
    elif impactor == 3:
        d1 = 0.73*d_mid
        #d1 = np.mean(d1_array[(d1_array >= 0.00015) & (d1_array <= 0.00025)])
        #d1 = np.percentile(d1_array, 10)
    #test = d1/d_mid
elif distribution == 2:
    d_mid = (d_min + d_max) * 0.5
    if impactor == 1:
        d1 = d_max
    elif impactor == 2:
        d1 = d_mid
    elif impactor == 3:
        d1 = d_min
elif distribution == 3:
    d1_array = generate_bimodal(mu1, sigma1, mu2, sigma2, weight1, d_min, d_max, num_samples)
    d_mid = np.percentile(d1_array, 50)
    if impactor == 1:
        d1 = np.percentile(d1_array, 90)
    elif impactor == 2:
        d1 = d_mid
    elif impactor == 3:
        d1 = np.percentile(d1_array, 10)
elif distribution == 4:
    d1_array = np.random.normal(loc=normal_E, scale=normal_D, size=num_samples)
    d_mid = np.percentile(d1_array, 50)
    if impactor == 1:
        d1 = np.percentile(d1_array, 90)
    elif impactor == 2:
        d1 = d_mid
    elif impactor == 3:
        d1 = np.percentile(d1_array, 10)
# impact velocity
if variation_param == 0:
    x_array = np.linspace(v1_start, v1_end, case_num)
else:
    v1 = v1_single
# impact angle
if variation_param == 1:
    x_array = np.linspace(theta1_start, theta1_end, case_num)
else:
    theta1 = theta1_single

rebound_ratio_list_0 = []
rebound_ratio_list_1 = []
rebound_ratio0_list_0 = []
e_bar_list_0 = []
e_bar_list_1 = []
e0_bar_list_0 = []
ez_bar_list_0 = []
ez_bar_list_1 = []
ex_bar_list_0 = []
ex_bar_list_1 = []
ez0_bar_list_0 = []
ex0_bar_list_0 = []
theta2_bar_list_0 = []
theta2_bar_list_1 = []
theta20_bar_list_0 = []
Nej_bar_list_0 = []
Nej_bar_list_1 = []
vn_bar_list_0 = []
vn_bar_list_1 = []
for n, x in enumerate(tqdm(x_array)):
    rebound_num_array = [0, 0, 0]
    rebound_num_array0 = [0, 0, 0]
    rebound_ratio_array = [0, 0, 0]
    rebound_ratio_array0 = [0, 0, 0]
    d3_list = []
    e_list_0 = []
    e_list_1 = []
    e0_list_0 = []
    ex_list_0 = []
    ex_list_1 = []
    ex0_list_0 = []
    ez_list_0 = []
    ez_list_1 = []
    ez0_list_0 = []
    theta2_list_0 = []
    theta2_list_1 = []
    theta20_list_0 = []
    Nej_list_0 = []
    Nej_list_1 = []
    vn_list_0 = []
    vn_list_1 = []
    while len(d3_list) < num_samples:
        if distribution == 0:
            d_array = np.random.uniform(d_min, d_max, size=3)
            if impactor == 0:
                d1 = d_array[0]
            d2 = d_array[1]
            d3 = d_array[2]
            d3_list.append(d3)
        elif distribution == 1:
            if not lognormal_param:
                mu, sigma = get_normal_params(normal_E, normal_D)
            d_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, 4)
            if impactor == 0:
                d1 = d_array[0]
            d2 = d_array[1]
            d3 = d_array[2]
            d4 = d_array[3]
            d3_list.append(d3)
        elif distribution == 2:
            r = np.random.uniform(0.5,1.5,size=3)
            rr = np.floor(r)
            d_bin = d_max - d_min
            if impactor == 0:
                d1 = d_min + rr[0]*d_bin
            d2 = d_min + rr[1]*d_bin
            d3 = d_min + rr[2]*d_bin
            d3_list.append(d3)
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
            if (d_min <= d1 <= d_max and
                d_min <= d2 <= d_max and
                d2 >= d2_min
                ):
                d3_list.append(d3)
        elif distribution == 4:
            d_array = np.random.normal(loc=normal_E, scale=normal_D, size=3)
            if impactor == 0:
                d1 = d_array[0]
            d2 = d_array[1]
            d3 = d_array[2]
            d3_list.append(d3)
        if variation_param == 0:
            v1 = x
            x_array[n] = v1/np.sqrt(g*d1)
        elif variation_param == 1:
            theta1 = x
        # d1, d2, d3分别为不同的直径
        current_eps = epsilon
        current_nu = nu

        d = (d1 + d2)/2
        d1_hat = d1/d
        d2_hat = d2/d
        d3_hat = d3/d
        d4_hat = d4/d
        # restitution coefficient
        mu_re = current_eps*d1_hat**3/(d1_hat**3 + current_eps*d2_hat**3)
        alpha = (1 + current_eps)/(1 + mu_re) - 1
        beta = 1 - (2/7)*(1 - current_nu)/(1 + mu_re)
        if shallow:
            x_min, psi1, psi2 = calculate_x_min_shallow(d1_hat, d2_hat, d3_hat, theta1)
            x_max = calculate_x_max_shallow(theta1)
        else:
            x_max = calculate_x_max(alpha, beta, theta1)
            if Three_D:
                x_min, psi1, psi2 = calculate_x_min_3D(d1_hat, d2_hat, d3_hat, d4_hat, theta1)
            else:
                x_min, psi1, psi2 = calculate_x_min(d1_hat, d2_hat, d3_hat, theta1)

        if shallow:
            if x_min < x_max:
                dont_change_flag = True
            else:
                dont_change_flag = False
        else:
            if x_min < x_max:
                dont_change_flag = True
            else:
                dont_change_flag = False

        if dont_change_flag:
            x0_hat = np.random.uniform(x_min, x_max)
            if shallow:
                if simplify:
                    e, evx, evz = calculate_e_bar_shallow_simp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1)
                    #theta2 = rebound_angle_shallow_simp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1)
                    theta2 = theta2_simp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1, x0_hat)
                    e = rebound_res_shallow_simp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1)
                else:
                    e, evx, evz = calculate_e_bar_shallow_comp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1, x_min, x_max)
                    theta2, error = quad(rebound_angle_shallow, x_min, x_max, args=(alpha, beta, theta1))
                    theta2 = theta2 / (x_max - x_min)
                    e, error = quad(rebound_res_shallow, x_min, x_max, args=(alpha, beta, theta1))
                    e = e / (x_max - x_min)
            else:
                e, evx, evz = calculate_e_bar(alpha, beta, x_min, x_max, theta1)
                theta2, error = quad(rebound_angle_bar, x_min, x_max, args=(alpha, beta, theta1))
                theta2 = theta2 / (x_max - x_min)
                e, error = quad(rebound_res_bar, x_min, x_max, args=(alpha, beta, theta1))
                e = e / (x_max - x_min)
            e0, evx0, evz0 = calculate_e(alpha, beta, x0_hat, theta1)
        else:
            d_temp = d2
            d2 = d3
            d3 = d_temp
            d = (d1 + d2)/2
            d1_hat = d1/d
            d2_hat = d2/d
            d3_hat = d3/d
            d4_hat = d4/d
            mu_re = current_eps*d1_hat**3/(d1_hat**3 + current_eps*d2_hat**3)
            alpha = (1 + current_eps)/(1 + mu_re) - 1
            beta = 1 - (2/7)*(1 - current_nu)/(1 + mu_re)
            if shallow:
                x_min, psi1, psi2 = calculate_x_min_shallow(d1_hat, d2_hat, d3_hat, theta1)
                x_max = calculate_x_max_shallow(theta1)
            else:
                x_max = calculate_x_max(alpha, beta, theta1)
                if Three_D:
                    x_min, psi1, psi2 = calculate_x_min_3D(d1_hat, d2_hat, d3_hat, d4_hat, theta1)
                else:
                    x_min, psi1, psi2 = calculate_x_min(d1_hat, d2_hat, d3_hat, theta1)

            x0_hat = np.random.uniform(x_min, x_max)
            if shallow:
                if simplify:
                    e, evx, evz = calculate_e_bar_shallow_simp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1)
                    #theta2 = rebound_angle_shallow_simp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1)
                    theta2 = theta2_simp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1, x0_hat)
                    e = rebound_res_shallow_simp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1)
                else:
                    e, evx, evz = calculate_e_bar_shallow_comp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1, x_min, x_max)
                    theta2, error = quad(rebound_angle_shallow, x_min, x_max, args=(alpha, beta, theta1))
                    theta2 = theta2 / (x_max - x_min)
                    e, error = quad(rebound_res_shallow, x_min, x_max, args=(alpha, beta, theta1))
                    e = e / (x_max - x_min)
            else:
                e, evx, evz = calculate_e_bar(alpha, beta, x_min, x_max, theta1)
                theta2, error = quad(rebound_angle_bar, x_min, x_max, args=(alpha, beta, theta1))
                theta2 = theta2 / (x_max - x_min)
                e, error = quad(rebound_res_bar, x_min, x_max, args=(alpha, beta, theta1))
                e = e / (x_max - x_min)
            e0, evx0, evz0 = calculate_e(alpha, beta, x0_hat, theta1)

        v2 = e*v1
        v20 = e0*v1
        v2_z = v1*evz
        v2_z0 = v1*evz0
        v2_x = v1*evx
        v2_x0 = v1*evx0

        E1 = 0.5*rho*(d1**3*np.pi/6)*v1**2
        Ec = rho*(d_mid**3*np.pi/6)*g*d_mid
        k_max = (1.0 - e0**2)*E1/Ec
        lambda_ej = 2.0*np.log(k_max)
        mu_ej = np.log((1.0 - e0**2)*E1) - lambda_ej*np.log(2.0)
        sigma_ej = np.sqrt(lambda_ej)*np.log(2.0)
        En_bar = np.exp(mu_ej + sigma_ej**2/2)
        vn = np.sqrt(2.0*En_bar/(rho*(d2**3*np.pi/6)))
        Nej = gamma_ej*(1.0 - e0**2)*E1/En_bar/2.0*erfc((np.log(Ec) - mu_ej)/(np.sqrt(2.0)*sigma_ej))
        Nej = max(Nej, 0.0)
        vn_list_0.append(vn)
        Nej_list_0.append(Nej)

        is_survive = calculate_survive(d1, d2, d3, psi1, psi2, g, v2_x, v2_z, e)
        is_survive0 = calculate_survive(d1, d2, d3, psi1, psi2, g, v2_x0, v2_z0, e0)
        #is_survive = calculate_survive1(d1, d2, v1, theta1, x0_hat, g, v2_x, v2_z, e0)
        if is_survive:
            rebound_num_array[0] += 1
        ez = evz/np.sin(theta1)
        ex = evx/np.cos(theta1)
        e_list_0.append(e)
        ez_list_0.append(ez)
        ex_list_0.append(ex)
        theta2_list_0.append(theta2)
        if is_survive0:
            rebound_num_array0[0] += 1
        ez0 = evz0/np.sin(theta1)
        ex0 = evx0/np.cos(theta1)
        theta20 = np.arctan2(evz0, evx0)
        e0_list_0.append(e0)
        ez0_list_0.append(ez0)
        ex0_list_0.append(ex0)
        theta20_list_0.append(theta20)

    # d1, d2=d3=d_mid
    d2 = d_mid
    d3 = d_mid
    d4 = d_mid
    current_eps = epsilon_La
    current_nu = nu_La

    d = (d1 + d2)/2
    d1_hat = d1/d
    d2_hat = d2/d
    d3_hat = d3/d
    d4_hat = d4/d
    # restitution coefficient
    mu_re = current_eps*d1_hat**3/(d1_hat**3 + current_eps*d2_hat**3)
    alpha = (1 + current_eps)/(1 + mu_re) - 1
    beta = 1 - (2/7)*(1 - current_nu)/(1 + mu_re)
    if shallow:
        x_min, psi1, psi2 = calculate_x_min_shallow(d1_hat, d2_hat, d3_hat, theta1)
        x_max = calculate_x_max_shallow(theta1)
    else:
        x_max = calculate_x_max(alpha, beta, theta1)
        if Three_D:
            x_min, psi1, psi2 = calculate_x_min_3D(d1_hat, d2_hat, d3_hat, d4_hat, theta1)
        else:
            x_min, psi1, psi2 = calculate_x_min(d1_hat, d2_hat, d3_hat, theta1)

    if shallow:
        if x_min < x_max:
            dont_change_flag = True
        else:
            dont_change_flag = False
    else:
        if x_min < x_max:
            dont_change_flag = True
        else:
            dont_change_flag = False

    e, evx, evz = calculate_e_bar(alpha, beta, x_min, x_max, theta1)
    theta2, error = quad(rebound_angle_bar, x_min, x_max, args=(alpha, beta, theta1))
    theta2 = theta2 / (x_max - x_min)
    e, error = quad(rebound_res_bar, x_min, x_max, args=(alpha, beta, theta1))
    e = e / (x_max - x_min)

    v2 = e*v1
    v2_z = v1*evz
    v2_x = v1*evx

    E1 = 0.5*rho*(d1**3*np.pi/6)*v1**2
    Ec = rho*(d_mid**3*np.pi/6)*g*d_mid
    k_max = (1.0 - e**2)*E1/Ec
    lambda_ej = 2.0*np.log(k_max)
    mu_ej = np.log((1.0 - e**2)*E1) - lambda_ej*np.log(2.0)
    sigma_ej = np.sqrt(lambda_ej)*np.log(2.0)
    En_bar = np.exp(mu_ej + sigma_ej**2/2)
    vn = np.sqrt(2.0*En_bar/(rho*(d2**3*np.pi/6)))
    Nej = gamma_ej*(1.0 - e**2)*E1/En_bar/2.0*erfc((np.log(Ec) - mu_ej)/(np.sqrt(2.0)*sigma_ej))
    Nej = max(Nej, 0.0)
    vn_list_1.append(vn)
    Nej_list_1.append(Nej)

    is_survive = calculate_survive(d1, d2, d3, psi1, psi2, g, v2_x, v2_z, e)
    if is_survive:
        rebound_num_array[1] += 1
    ez = evz/np.sin(theta1)
    ex = evx/np.cos(theta1)
    e_list_1.append(e)
    ez_list_1.append(ez)
    ex_list_1.append(ex)
    theta2_list_1.append(theta2)
    """d1, d2, d3"""
    e_bar_0 = np.mean(e_list_0)
    ez_bar_0 = np.mean(ez_list_0)
    ex_bar_0 = np.mean(ex_list_0)
    theta2_bar_0 = np.mean(theta2_list_0)
    e_bar_list_0.append(e_bar_0)
    ez_bar_list_0.append(ez_bar_0)
    ex_bar_list_0.append(ex_bar_0)
    theta2_bar_list_0.append(theta2_bar_0)
    e0_bar_0 = np.mean(e0_list_0)
    ez0_bar_0 = np.mean(ez0_list_0)
    ex0_bar_0 = np.mean(ex0_list_0)
    #e0_bar_0 = np.sqrt(ex0_bar_0**2 + ez0_bar_0**2)
    theta20_bar_0 = np.mean(theta20_list_0)
    e0_bar_list_0.append(e0_bar_0)
    ez0_bar_list_0.append(ez0_bar_0)
    ex0_bar_list_0.append(ex0_bar_0)
    theta20_bar_list_0.append(theta20_bar_0)
    """d2=d3"""
    e_bar_1 = np.mean(e_list_1)
    ez_bar_1 = np.mean(ez_list_1)
    ex_bar_1 = np.mean(ex_list_1)
    theta2_bar_1 = np.mean(theta2_list_1)
    e_bar_list_1.append(e_bar_1)
    ez_bar_list_1.append(ez_bar_1)
    ex_bar_list_1.append(ex_bar_1)
    theta2_bar_list_1.append(theta2_bar_1)
    """ejection velocity"""
    vn_bar_0 = np.mean(vn_list_0)/np.sqrt(g*d_mid)
    vn_bar_1 = np.mean(vn_list_1)/np.sqrt(g*d_mid)
    vn_bar_list_0.append(vn_bar_0)
    vn_bar_list_1.append(vn_bar_1)
    """ejection number"""
    Nej_bar_0 = np.mean(Nej_list_0)
    Nej_bar_1 = np.mean(Nej_list_1)
    Nej_bar_list_0.append(Nej_bar_0)
    Nej_bar_list_1.append(Nej_bar_1)
    """rebound ratio"""
    rebound_ratio_array = [rebound_num_array[i]/num_samples for i in range(3)]
    rebound_ratio_array0 = [rebound_num_array0[i]/num_samples for i in range(3)]
    rebound_ratio_list_0.append(rebound_ratio_array[0])
    rebound_ratio_list_1.append(rebound_ratio_array[1])
    rebound_ratio0_list_0.append(rebound_ratio_array0[0])

"""other's data"""
Beladjine07_v26= {
    'ang_in': [10, 20, 40, 60, 90],
    'ang_re': [21.21, 26.68, 33.87, 40.94, 90.01],
    'e': [0.77, 0.61, 0.43, 0.26, 0.22]
}
Zhou06_v_many= {
    'ang_in': [8, 11.5],
    'ang_re_mean': [46.97, 47.53],
    'ang_re_std': [0.59, 0.75],
    'e_mean': [0.63, 0.61],
    'e_std': [0.006, 0.009]
}
Rice95_v_many_coarse= {
    'ang_in': np.mean([13.94, 14.75, 14.73, 15.04]),
    'ang_in_std': np.std([13.94, 14.75, 14.73, 15.04]),
    'ang_re': np.mean([20.95, 23.03, 22.55, 25.63]),
    'ang_re_std': np.std([20.95, 23.03, 22.55, 25.63]),
    'e': np.mean([0.58, 0.56, 0.56, 0.58]),
    'e_std': np.std([0.58, 0.56, 0.56, 0.58]),
    'Nej': np.mean([5.7, 5.09, 5.12, 4.64]),
    'Nej_std': np.std([5.7, 5.09, 5.12, 4.64])
}
Rice95_v_many_medium= {
    'ang_in': np.mean([11.82, 11.47, 11.53, 11.36]),
    'ang_in_std': np.std([11.82, 11.47, 11.53, 11.36]),
    'ang_re': np.mean([29.95, 30.31, 31.06, 29.56]),
    'ang_re_std': np.std([29.95, 30.31, 31.06, 29.56]),
    'e': np.mean([0.57, 0.54, 0.59, 0.58]),
    'e_std': np.std([0.57, 0.54, 0.59, 0.58]),
    'Nej': np.mean([2.68, 3.43, 2.67, 2.76]),
    'Nej_std': np.std([2.68, 3.43, 2.67, 2.76])
}
Rice95_v_many_fine= {
    'ang_in': np.mean([10.85, 10.24, 10.46]),
    'ang_in_std': np.std([10.85, 10.24, 10.46]),
    'ang_re': np.mean([44.63, 38.03, 37.85]),
    'ang_re_std': np.std([44.63, 38.03, 37.85]),
    'e': np.mean([0.52, 0.56, 0.58]),
    'e_std': np.std([0.52, 0.56, 0.58]),
    'Nej': np.mean([1.86, 1.71, 1.58]),
    'Nej_std': np.std([1.86, 1.71, 1.58])
}
Chen18_v_many= {
    'ang_in': [23.2, 21.8, 21.9, 30.7, 30.6, 37.6, 47.1, 46.6, 46],
    'ang_re': [32.83, 35.59, 29.90, 43, 33.28, 44.16, 57.62, 50.97, 39.88],
    'e': [0.38, 0.4, 0.45, 0.32, 0.37, 0.27, 0.2, 0.19, 0.22]
}
Willetts89_v_many_coarse= {
    'ang_in': [12.7, 17.8, 23.2, 27.7],
    'ang_re': [19.1, 25.2, 21.4, 27.2],
    'e': [0.63, 0.57, 0.54, 0.46]
}
Willetts89_v_many_medium= {
    'ang_in': [11.7, 18.2, 21.4, 26.3],
    'ang_re': [24.9, 33.4, 33.3, 44.7],
    'e': [0.61, 0.53, 0.50, 0.40]
}
Willetts89_v_many_fine= {
    'ang_in': [9.5, 15.4, 19.7, 24.9],
    'ang_re': [38.8, 42, 42.2, 42.5],
    'e': [0.57, 0.50, 0.48, 0.46]
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

if impactor == 1:
    output_file_e = 'rebound_e_coarse.txt'
    output_file_th = 'rebound_th_coarse.txt'
elif impactor == 2:
    output_file_e = 'rebound_e_medium.txt'
    output_file_th = 'rebound_th_medium.txt'
elif impactor == 3:
    output_file_e = 'rebound_e_fine.txt'
    output_file_th = 'rebound_th_fine.txt'

if variation_param == 1:
    x_array = np.rad2deg(x_array)

if output_e:
    plt.figure(1)
    if variation_param == 0:
        plt.xlabel('$\\hat{v}_1$')
        x_exp = 'v_in'
    elif variation_param == 1:
        plt.xlabel('$\\theta (degree)$')
        x_exp = 'ang_in'
    plt.ylabel('$e$')
    strlabel1 = '$e$'
    strlabel2 = '$\\overline{e}$'
    plt.plot(x_array, e0_bar_list_0, 'r-', label=strlabel1 + ': d1, d2, d3')
    plt.plot(x_array, e_bar_list_0, 'r--', label=strlabel2 + ': d1, d2, d3')
    plt.plot(x_array, e_bar_list_1, 'b--', label=strlabel2 + ': d1, d2=d3')
    plt.xlim(0, x_array[-1])
    with open(output_file_e, 'w') as f:
        for i in range(len(x_array)):
            f.write(f"{x_array[i]:.6f} {e0_bar_list_0[i]:.6f} {e_bar_list_0[i]:.6f} {e_bar_list_1[i]:.6f}\n")

    if impactor == 0 or impactor == 2:
        #plt.plot(Chen18_v_many[x_exp], Chen18_v_many['e'], 'kP', label='Chen18')
        #plt.plot(Beladjine07_v26[x_exp], Beladjine07_v26['e'], 'k*', label='Beladjine07')
        #plt.plot(Zhou06_v_many[x_exp], Zhou06_v_many['e_mean'], 'k^', label='Zhou06')
        #plt.plot(Rice95_v_many_medium[x_exp], Rice95_v_many_medium['e'], 'kD', label='Rice95 medium')
        plt.errorbar(
            Rice95_v_many_medium[x_exp],
            Rice95_v_many_medium['e'],
            yerr=Rice95_v_many_medium['e_std'],
            fmt='kD',
            label='Rice95 medium',
            capsize=5
        )
        plt.plot(Willetts89_v_many_medium[x_exp], Willetts89_v_many_medium['e'], 'kX', label='Willetts89 medium')
        #plt.plot(Rioual20_v_many[x_exp], Rioual20_v_many['e'], 'k.', label='Rioual20')
        #plt.plot(Gordon21_v_many[x_exp], Gordon21_v_many['e'], 'k1', label='Gordon21')
        #plt.plot(Gordon09_v_many[x_exp], Gordon09_v_many['e'], 'k2', label='Gordon09')
    elif impactor == 1:
        #plt.plot(Rice95_v_many_coarse[x_exp], Rice95_v_many_coarse['e'], 'ks', label='Rice95 coarse')
        plt.errorbar(
            Rice95_v_many_coarse[x_exp],
            Rice95_v_many_coarse['e'],
            yerr=Rice95_v_many_coarse['e_std'],
            fmt='ks',
            label='Rice95 coarse',
            capsize=5
        )
        plt.plot(Willetts89_v_many_coarse[x_exp], Willetts89_v_many_coarse['e'], 'kH', label='Willetts89 coarse')
    elif impactor == 3:
        #plt.plot(Rice95_v_many_fine[x_exp], Rice95_v_many_fine['e'], 'ko', label='Rice95 fine')
        plt.errorbar(
            Rice95_v_many_fine[x_exp],
            Rice95_v_many_fine['e'],
            yerr=Rice95_v_many_fine['e_std'],
            fmt='ko',
            label='Rice95 fine',
            capsize=5
        )
        plt.plot(Willetts89_v_many_fine[x_exp], Willetts89_v_many_fine['e'], 'k+', label='Willetts89 fine')

    plt.legend()

if output_ez:
    plt.figure(2)
    if variation_param == 0:
        plt.xlabel('$\\hat{v}_1$')
    elif variation_param == 1:
        plt.xlabel('$\\theta (degree)$')
    plt.ylabel('$e_z$')
    strlabel1 = '$e_z$'
    strlabel2 = '$\\overline{e_z}$'
    plt.plot(x_array, ez0_bar_list_0, 'r-', label=strlabel1 + ': d1, d2, d3')
    plt.plot(x_array, ez_bar_list_0, 'r--', label=strlabel2 + ': d1, d2, d3')
    plt.plot(x_array, ez_bar_list_1, 'b--', label=strlabel2 + ': d1, d2=d3')
    plt.legend()

if output_ex:
    plt.figure(3)
    if variation_param == 0:
        plt.xlabel('$\\hat{v}_1$')
    elif variation_param == 1:
        plt.xlabel('$\\theta (degree)$')
    plt.ylabel('$e_x$')
    strlabel1 = '$e_x$'
    strlabel2 = '$\\overline{e_x}$'
    plt.plot(x_array, ex0_bar_list_0, 'r-', label=strlabel1 + ': d1, d2, d3')
    plt.plot(x_array, ex_bar_list_0, 'r--', label=strlabel2 + ': d1, d2, d3')
    plt.plot(x_array, ex_bar_list_1, 'b--', label=strlabel2 + ': d1, d2=d3')
    plt.legend()

if output_rebound_ratio:
    plt.figure(4)
    if variation_param == 0:
        plt.xlabel('$\\hat{v}_1$')
        plt.semilogx(x_array, rebound_ratio0_list_0, 'r-', label='e: d1, d2, d3')
        plt.semilogx(x_array, rebound_ratio_list_0, 'r--', label='ebar: d1, d2, d3')
        plt.semilogx(x_array, rebound_ratio_list_1, 'b--', label='ebar: d1, d2=d3')
    elif variation_param == 1:
        plt.xlabel('$\\theta (degree)$')
        plt.plot(x_array, rebound_ratio0_list_0, 'r--', label='e: d1, d2, d3')
        plt.plot(x_array, rebound_ratio_list_0, 'r-', label='ebar: d1, d2, d3')
        plt.plot(x_array, rebound_ratio_list_1, 'b-', label='ebar: d1, d2=d3')
    plt.ylabel('Rebound Ratio')
    plt.legend()

if output_theta2_dist:
    plt.figure(5)
    theta20_list_0 = np.rad2deg(theta20_list_0)
    sns.histplot(theta20_list_0, bins=40, kde=True, color='red', stat='density', label='d1, d2, d3', element='step', alpha=0.1)
    plt.xlabel("$\\theta \'$", fontsize=16)
    plt.ylabel('pdf', fontsize=16)
    plt.legend()

if output_theta2:
    plt.figure(6)
    if variation_param == 0:
        plt.xlabel('$\\hat{v}_1$')
    elif variation_param == 1:
        plt.xlabel('$\\theta (degree)$')
    plt.ylabel('$\\theta \'$')
    strlabel1 = '$\\theta \'$'
    strlabel2 = '$\\overline{\\theta \'}$'
    theta20_bar_list_0 = np.rad2deg(theta20_bar_list_0)
    theta2_bar_list_0 = np.rad2deg(theta2_bar_list_0)
    theta2_bar_list_1 = np.rad2deg(theta2_bar_list_1)
    plt.plot(x_array, theta20_bar_list_0, 'r-', label=strlabel1 + ': d1, d2, d3')
    plt.plot(x_array, theta2_bar_list_0, 'r--', label=strlabel2 + ': d1, d2, d3')
    plt.plot(x_array, theta2_bar_list_1, 'b--', label=strlabel2 + ': d1, d2=d3')
    plt.xlim(0, x_array[-1])
    with open(output_file_th, 'w') as f:
        for i in range(len(x_array)):
            f.write(f"{x_array[i]:.6f} {theta20_bar_list_0[i]:.6f} {theta2_bar_list_0[i]:.6f} {theta2_bar_list_1[i]:.6f}\n")

    if impactor == 0 or impactor == 2:
        #plt.plot(Chen18_v_many[x_exp], Chen18_v_many['ang_re'], 'kP', label='Chen et al. 2018')
        #plt.plot(Beladjine07_v26[x_exp], Beladjine07_v26['ang_re'], 'k*', label='Beladjine et al. 2007')
        #plt.plot(Zhou06_v_many[x_exp], Zhou06_v_many['ang_re_mean'], 'k^', label='Zhou et al. 2006')
        #plt.plot(Rice95_v_many_medium[x_exp], Rice95_v_many_medium['ang_re'], 'ko', label='Rice et al. 1995 (medium)')
        plt.errorbar(
            Rice95_v_many_medium[x_exp],
            Rice95_v_many_medium['ang_re'],
            yerr=Rice95_v_many_medium['ang_re_std'],
            fmt='kD',
            label='Rice et al. 1995 (medium)',
            capsize=5
        )
        plt.plot(Willetts89_v_many_medium[x_exp], Willetts89_v_many_medium['ang_re'], 'kX', label='Willetts et al. 1989 (medium)')
    elif impactor == 1:
        #plt.plot(Rice95_v_many_coarse[x_exp], Rice95_v_many_coarse['ang_re'], 'kD', label='Rice et al. 1995 (coarse)')
        plt.errorbar(
            Rice95_v_many_coarse[x_exp],
            Rice95_v_many_coarse['ang_re'],
            yerr=Rice95_v_many_coarse['ang_re_std'],
            fmt='kD',
            label='Rice et al. 1995 (coarse)',
            capsize=5
        )
        plt.plot(Willetts89_v_many_coarse[x_exp], Willetts89_v_many_coarse['ang_re'], 'kH', label='Willetts et al. 1989 (coarse)')
    elif impactor == 3:
        #plt.plot(Rice95_v_many_fine[x_exp], Rice95_v_many_fine['ang_re'], 'ks', label='Rice et al. 1995 (fine)')
        plt.errorbar(
            Rice95_v_many_fine[x_exp],
            Rice95_v_many_fine['ang_re'],
            yerr=Rice95_v_many_fine['ang_re_std'],
            fmt='ks',
            label='Rice et al. 1995 (fine)',
            capsize=5
        )
        plt.plot(Willetts89_v_many_fine[x_exp], Willetts89_v_many_fine['ang_re'], 'k+', label='Willetts et al. 1989 (fine)')
    plt.legend()

if output_Nej:
    plt.figure(7)
    if variation_param == 0:
        plt.xlabel('$\\hat{v}_1$')
    elif variation_param == 1:
        plt.xlabel('$\\theta (degree)$')
    plt.ylabel('$N_{ej}$')
    plt.plot(x_array, Nej_bar_list_0, 'r-', label='d1, d2, d3')
    plt.plot(x_array, Nej_bar_list_1, 'b-', label='d1, d2=d3')

    if impactor == 0 or impactor == 2:
        plt.plot(Rice95_v_many_medium[x_exp], Rice95_v_many_medium['Nej'], 'ko', label='Rice95 medium')
    if impactor == 1:
        plt.plot(Rice95_v_many_coarse[x_exp], Rice95_v_many_coarse['Nej'], 'kD', label='Rice95 coarse')
    if impactor == 3:
        plt.plot(Rice95_v_many_fine[x_exp], Rice95_v_many_fine['Nej'], 'ks', label='Rice95 fine')
    plt.legend()

if output_vn:
    plt.figure(8)
    if variation_param == 0:
        plt.xlabel('$\\hat{v}_1$')
    elif variation_param == 1:
        plt.xlabel('$\\theta (degree)$')
    plt.ylabel('$v_n$')
    plt.plot(x_array, vn_bar_list_0, 'r-', label='d1, d2, d3')
    plt.plot(x_array, vn_bar_list_1, 'b-', label='d1, d2=d3')
    plt.legend()

plt.show()