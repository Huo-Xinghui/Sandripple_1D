import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import truncnorm
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

def calculate_rebound_angle(d2_hat, theta1, alpha, beta):
    gama = 4.0/9.0*beta**2/(alpha + beta)**2/d2_hat
    theta2_min = -theta1
    theta2_max = np.sqrt(theta1/gama)*2.0 - theta1
    theta2_max = min(theta2_max, np.pi)
    theta2_bin = (theta2_max - theta2_min)/theta2_num
    theta2_hist = np.zeros(theta2_num)
    for i in range(theta2_num):
        x = theta2_min + theta2_bin*(i + 0.5)
        theta2_hist[i] = gama*(x + theta1)/theta1*np.log(2.0*theta1/gama/(x + theta1)**2)
        theta2_hist[i] = max(theta2_hist[i], 0.0)
    theta2_hist = theta2_hist/sum(theta2_hist)
    random_number = np.random.rand()
    cumulative_prob = 0.0
    for i in range(theta2_num):
        cumulative_prob += theta2_hist[i]
        if random_number <= cumulative_prob:
            cumulative_prob -= theta2_hist[i]
            theta2 = theta2_bin*(i + 0.5) + (random_number - cumulative_prob)/theta2_hist[i]*theta2_bin
            break
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
output_vn = True # output eject velocity
output_Nej = True # output eject number

just_survive = False # only calculate survival

distribution = 1 # 0:uniform, 1:lognormal, 2:bidisperse, 3:polydisperse, 4:normal
variation_param = 1 # 0: v1, 1: theta1
shallow = False # shallow impact
simplify = False # first order approximation
lognormal_param = True # lognormal distribution parameters
impactor = 2 #impactor type: 0: bed, 1: coarse, 2: medium, 3: fine, 4: dist

d_min = 0.5e-4
d_max = 6e-4
d_min_impactor = 3.55e-4
d_max_impactor = 6e-4
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
epsilon = 0.6
nu = -0.3
#epsilon = 0.73
#nu = -0.18
v1_single = 2.7323
theta1_single = np.pi/18
v1_start = 1
v1_end = 10
case_num = 100
theta1_start = np.pi/180
if shallow:
    theta1_end = np.pi/6
else:
    theta1_end = np.pi/2
theta2_num = 60
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
        #d1 = 1.4*np.percentile(d1_array, 50)
        d1 = np.percentile(d1_array, 90)
    elif impactor == 2:
        d1 = d_mid
    elif impactor == 3:
        #d1 = 0.73*np.percentile(d1_array, 50)
        d1 = np.percentile(d1_array, 10)
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
    x_array = x_array/np.sqrt(g*d1)
else:
    v1_hat_single = v1_single/np.sqrt(g*d1)
    v1_hat = v1_hat_single
# impact angle
if variation_param == 1:
    x_array = np.linspace(theta1_start, theta1_end, case_num)
else:
    theta1 = theta1_single

d2_list = []
rebound_ratio_list_0 = []
rebound_ratio_list_1 = []
rebound_ratio_list_2 = []
rebound_ratio0_list_0 = []
rebound_ratio0_list_1 = []
rebound_ratio0_list_2 = []
e_bar_list_0 = []
e_bar_list_1 = []
e_bar_list_2 = []
e0_bar_list_0 = []
e0_bar_list_1 = []
e0_bar_list_2 = []
ez_bar_list_0 = []
ez_bar_list_1 = []
ez_bar_list_2 = []
ex_bar_list_0 = []
ex_bar_list_1 = []
ex_bar_list_2 = []
ez0_bar_list_0 = []
ez0_bar_list_1 = []
ez0_bar_list_2 = []
ex0_bar_list_0 = []
ex0_bar_list_1 = []
ex0_bar_list_2 = []
theta2_bar_list_0 = []
theta2_bar_list_1 = []
theta2_bar_list_2 = []
theta20_bar_list_0 = []
theta20_bar_list_1 = []
theta20_bar_list_2 = []
Nej_bar_list_0 = []
Nej_bar_list_1 = []
Nej_bar_list_2 = []
vn_bar_list_0 = []
vn_bar_list_1 = []
vn_bar_list_2 = []
x_array0 = []
x_array1 = []
x_array2 = []
x0_array0 = []
x0_array1 = []
x0_array2 = []
for x in tqdm(x_array):
    if variation_param == 0:
        v1_hat = x
    elif variation_param == 1:
        theta1 = x
    rebound_num_array = [0, 0, 0]
    rebound_num_array0 = [0, 0, 0]
    rebound_ratio_array = [0, 0, 0]
    rebound_ratio_array0 = [0, 0, 0]
    d3_list = []
    e_list_0 = []
    e_list_1 = []
    e_list_2 = []
    e0_list_0 = []
    e0_list_1 = []
    e0_list_2 = []
    ex_list_0 = []
    ex_list_1 = []
    ex_list_2 = []
    ex0_list_0 = []
    ex0_list_1 = []
    ex0_list_2 = []
    ez_list_0 = []
    ez_list_1 = []
    ez_list_2 = []
    ez0_list_0 = []
    ez0_list_1 = []
    ez0_list_2 = []
    theta2_list_0 = []
    theta2_list_1 = []
    theta2_list_2 = []
    theta20_list_0 = []
    theta20_list_1 = []
    theta20_list_2 = []
    Nej_list_0 = []
    Nej_list_1 = []
    Nej_list_2 = []
    vn_list_0 = []
    vn_list_1 = []
    vn_list_2 = []
    have_survived = [False, False, False]
    have_survived0 = [False, False, False]
    while len(d3_list) < num_samples:
        if distribution == 0:
            d_array = np.random.uniform(d_min, d_max, size=3)
            if impactor == 0:
                d1 = d_array[0]
            elif impactor == 4:
                d1 = np.random.uniform(d_min_impactor, d_max_impactor, size=1)[0]
            d2 = d_array[1]
            d3 = d_array[2]
            d2_list.append(d2)
            d3_list.append(d3)
        elif distribution == 1:
            mu, sigma = get_normal_params(normal_E, normal_D)
            d_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, 3)
            if impactor == 0:
                d1 = d_array[0]
            elif impactor == 4:
                d1 = generate_truncated_lognormal(mu, sigma, d_min_impactor, d_max_impactor, 1)[0]
            d2 = d_array[1]
            d3 = d_array[2]
            d2_list.append(d2)
            d3_list.append(d3)
        elif distribution == 2:
            r = np.random.uniform(0.5,1.5,size=3)
            rr = np.floor(r)
            d_bin = d_max - d_min
            if impactor == 0 or impactor == 4:
                d1 = d_min + rr[0]*d_bin
            d2 = d_min + rr[1]*d_bin
            d3 = d_min + rr[2]*d_bin
            d2_list.append(d2)
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
                d2_list.append(d2)
                d3_list.append(d3)
        elif distribution == 4:
            d_array = np.random.normal(loc=normal_E, scale=normal_D, size=3)
            if impactor == 0 or impactor == 4:
                d1 = d_array[0]
            d2 = d_array[1]
            d3 = d_array[2]
            d2_list.append(d2)
            d3_list.append(d3)
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
            mu = epsilon*d1_hat**3/(d1_hat**3 + epsilon*d2_hat**3)
            alpha = (1 + epsilon)/(1 + mu) - 1
            beta = 1 - (2/7)*(1 - nu)/(1 + mu)
            x_min, psi1, psi2 = calculate_x_min(d1_hat, d2_hat, d3_hat, theta1)

            if shallow:
                x_max = calculate_x_max_shallow(theta1)
                if x_min < x_max and d2 >= d3:
                    dont_change_flag = True
                else:
                    dont_change_flag = False
            else:
                x_max = calculate_x_max(alpha, beta, theta1)
                if x_min < x_max:
                    dont_change_flag = True
                else:
                    dont_change_flag = False

            if dont_change_flag:
                if shallow:
                    if simplify:
                        e, evx, evz = calculate_e_bar_shallow_simp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1)
                    else:
                        e, evx, evz = calculate_e_bar_shallow_comp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1, x_min, x_max)
                else:
                    e, evx, evz = calculate_e_bar(alpha, beta, x_min, x_max, theta1)
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
                mu = epsilon*d1_hat**3/(d1_hat**3 + epsilon*d2_hat**3)
                alpha = (1 + epsilon)/(1 + mu) - 1
                beta = 1 - (2/7)*(1 - nu)/(1 + mu)
                x_min, psi1, psi2 = calculate_x_min(d1_hat, d2_hat, d3_hat, theta1)

                if shallow:
                    x_max = calculate_x_max_shallow(theta1)
                    if simplify:
                        e, evx, evz = calculate_e_bar_shallow_simp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1)
                    else:
                        e, evx, evz = calculate_e_bar_shallow_comp(alpha, beta, d1_hat, d2_hat, d3_hat, theta1, x_min, x_max)
                else:
                    x_max = calculate_x_max(alpha, beta, theta1)
                    e, evx, evz = calculate_e_bar(alpha, beta, x_min, x_max, theta1)

                x0_hat = np.random.uniform(x_min, x_max)
                e0, evx0, evz0 = calculate_e(alpha, beta, x0_hat, theta1)

            v1 = v1_hat*np.sqrt(g*d1)
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
            if i == 0:
                vn_list_0.append(vn)
                Nej_list_0.append(Nej)
            elif i == 1:
                vn_list_1.append(vn)
                Nej_list_1.append(Nej)
            elif i == 2:
                vn_list_2.append(vn)
                Nej_list_2.append(Nej)

            is_survive = calculate_survive(d1, d2, d3, psi1, psi2, g, v2_x, v2_z, e)
            is_survive0 = calculate_survive(d1, d2, d3, psi1, psi2, g, v2_x0, v2_z0, e0)
            #is_survive = calculate_survive1(d1, d2, v1, theta1, x0_hat, g, v2_x, v2_z, e0)
            if is_survive:
                rebound_num_array[i] += 1
            if (just_survive and is_survive) or not just_survive:
                have_survived[i] = True
                ez = evz #/np.sin(theta1)
                ex = evx #/np.cos(theta1)
                theta2 = np.arctan(evz/evx)
                if i == 0:
                    e_list_0.append(e)
                    ez_list_0.append(ez)
                    ex_list_0.append(ex)
                    theta2_list_0.append(theta2)
                elif i == 1:
                    e_list_1.append(e)
                    ez_list_1.append(ez)
                    ex_list_1.append(ex)
                    theta2_list_1.append(theta2)
                elif i == 2:
                    e_list_2.append(e)
                    ez_list_2.append(ez)
                    ex_list_2.append(ex)
                    theta2_list_2.append(theta2)
                v2 = e*v1
                v2_z = v1*evz
                v2_x = v1*evx
            if is_survive0:
                rebound_num_array0[i] += 1
            if (just_survive and is_survive0) or not just_survive:
                have_survived0[i] = True
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
    if have_survived[0]:
        x_array0.append(x)
        e_bar_0 = np.mean(e_list_0)
        ez_bar_0 = np.mean(ez_list_0)
        ex_bar_0 = np.mean(ex_list_0)
        theta2_bar_0 = np.mean(theta2_list_0)
        e_bar_list_0.append(e_bar_0)
        ez_bar_list_0.append(ez_bar_0)
        ex_bar_list_0.append(ex_bar_0)
        theta2_bar_list_0.append(theta2_bar_0)
    if have_survived0[0]:
        x0_array0.append(x)
        ez0_bar_0 = np.mean(ez0_list_0)
        ex0_bar_0 = np.mean(ex0_list_0)
        e0_bar_0 = np.sqrt(ex0_bar_0**2 + ez0_bar_0**2)
        theta20_bar_0 = np.mean(theta20_list_0)
        e0_bar_list_0.append(e0_bar_0)
        ez0_bar_list_0.append(ez0_bar_0)
        ex0_bar_list_0.append(ex0_bar_0)
        theta20_bar_list_0.append(theta20_bar_0)
    """d2=d3"""
    if have_survived[1]:
        x_array1.append(x)
        e_bar_1 = np.mean(e_list_1)
        ez_bar_1 = np.mean(ez_list_1)
        ex_bar_1 = np.mean(ex_list_1)
        theta2_bar_1 = np.mean(theta2_list_1)
        e_bar_list_1.append(e_bar_1)
        ez_bar_list_1.append(ez_bar_1)
        ex_bar_list_1.append(ex_bar_1)
        theta2_bar_list_1.append(theta2_bar_1)
    if have_survived0[1]:
        x0_array1.append(x)
        ez0_bar_1 = np.mean(ez0_list_1)
        ex0_bar_1 = np.mean(ex0_list_1)
        e0_bar_1 = np.sqrt(ex0_bar_1**2 + ez0_bar_1**2)
        theta20_bar_1 = np.mean(theta20_list_1)
        e0_bar_list_1.append(e0_bar_1)
        ez0_bar_list_1.append(ez0_bar_1)
        ex0_bar_list_1.append(ex0_bar_1)
        theta20_bar_list_1.append(theta20_bar_1)
    """d1=d2=d3"""
    if have_survived[2]:
        x_array2.append(x)
        e_bar_2 = np.mean(e_list_2)
        ez_bar_2 = np.mean(ez_list_2)
        ex_bar_2 = np.mean(ex_list_2)
        theta2_bar_2 = np.mean(theta2_list_2)
        e_bar_list_2.append(e_bar_2)
        ez_bar_list_2.append(ez_bar_2)
        ex_bar_list_2.append(ex_bar_2)
        theta2_bar_list_2.append(theta2_bar_2)
    if have_survived0[2]:
        x0_array2.append(x)
        ez0_bar_2 = np.mean(ez0_list_2)
        ex0_bar_2 = np.mean(ex0_list_2)
        e0_bar_2 = np.sqrt(ex0_bar_2**2 + ez0_bar_2**2)
        theta20_bar_2 = np.mean(theta20_list_2)
        e0_bar_list_2.append(e0_bar_2)
        ez0_bar_list_2.append(ez0_bar_2)
        ex0_bar_list_2.append(ex0_bar_2)
        theta20_bar_list_2.append(theta20_bar_2)
    #e0_bar_0 = np.mean(e0_list_0)
    #e0_bar_1 = np.mean(e0_list_1)
    #e0_bar_2 = np.mean(e0_list_2)
    #e_bar_0 = np.sqrt(ex_bar_0**2 + ez_bar_0**2)
    #e_bar_1 = np.sqrt(ex_bar_1**2 + ez_bar_1**2)
    #e_bar_2 = np.sqrt(ex_bar_2**2 + ez_bar_2**2)
    #theta2_bar_0 = np.arctan(ez_bar_0/ex_bar_0)
    #theta2_bar_1 = np.arctan(ez_bar_1/ex_bar_1)
    #theta2_bar_2 = np.arctan(ez_bar_2/ex_bar_2)
    #theta20_bar_0 = np.arctan(ez0_bar_0/ex0_bar_0)
    #theta20_bar_1 = np.arctan(ez0_bar_1/ex0_bar_1)
    #theta20_bar_2 = np.arctan(ez0_bar_2/ex0_bar_2)
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
    """rebound ratio"""
    rebound_ratio_array = [rebound_num_array[i]/num_samples for i in range(3)]
    rebound_ratio_array0 = [rebound_num_array0[i]/num_samples for i in range(3)]
    rebound_ratio_list_0.append(rebound_ratio_array[0])
    rebound_ratio_list_1.append(rebound_ratio_array[1])
    rebound_ratio_list_2.append(rebound_ratio_array[2])
    rebound_ratio0_list_0.append(rebound_ratio_array0[0])
    rebound_ratio0_list_1.append(rebound_ratio_array0[1])
    rebound_ratio0_list_2.append(rebound_ratio_array0[2])

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
    'ang_in': [13.94, 14.75, 14.73, 15.04],
    'ang_re': [20.95, 23.03, 22.55, 25.63],
    'e': [0.58, 0.56, 0.56, 0.58],
    'Nej': [5.7, 5.09, 5.12, 4.64]
}
Rice95_v_many_medium= {
    'ang_in': [11.82, 11.47, 11.53, 11.36],
    'ang_re': [29.95, 30.31, 31.06, 29.56],
    'e': [0.57, 0.54, 0.59, 0.58],
    'Nej': [2.68, 3.43, 2.67, 2.76]
}
Rice95_v_many_fine= {
    'ang_in': [10.85, 10.24, 10.46],
    'ang_re': [44.63, 38.03, 37.85],
    'e': [0.52, 0.56, 0.58],
    'Nej': [1.86, 1.71, 1.58]
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

if variation_param == 1:
    x_array = np.rad2deg(x_array)
    x_array0 = np.rad2deg(x_array0)
    x_array1 = np.rad2deg(x_array1)
    x_array2 = np.rad2deg(x_array2)
    x0_array0 = np.rad2deg(x0_array0)
    x0_array1 = np.rad2deg(x0_array1)
    x0_array2 = np.rad2deg(x0_array2)

if output_e:
    plt.figure(1)
    if variation_param == 0:
        plt.xlabel(r'$\hat{v}_1$')
        plt.ylabel(r'$\langle e \rangle_{\theta}, \langle E \rangle_{\theta}$')
        strlabel1 = r'$\langle e \rangle_{\theta}$'
        strlabel2 = r'$\langle E \rangle_{\theta}$'
        x_exp = 'v_in'
    elif variation_param == 1:
        plt.xlabel(r'$\theta (degree)$')
        plt.ylabel(r'$\langle e \rangle_{\hat{v}}, \langle E \rangle_{\hat{v}}$')
        strlabel1 = r'$\langle e \rangle_{\hat{v}}$'
        strlabel2 = r'$\langle E \rangle_{\hat{v}}$'
        x_exp = 'ang_in'
    plt.plot(x0_array0, e0_bar_list_0, 'r-', label=strlabel1 + ': d1, d2, d3')
    plt.plot(x0_array1, e0_bar_list_1, 'b-', label=strlabel1 + ': d1, d2=d3')
    plt.plot(x0_array2, e0_bar_list_2, 'g-', label=strlabel1 + ': d1=d2=d3')
    #plt.plot(x_array0, e_bar_list_0, 'r--', label=strlabel2 + ': d1, d2, d3')
    #plt.plot(x_array1, e_bar_list_1, 'b--', label=strlabel2 + ': d1, d2=d3')
    #plt.plot(x_array2, e_bar_list_2, 'g--', label=strlabel2 + ': d1=d2=d3')

    if impactor == 0 or impactor == 2:
        plt.plot(Chen18_v_many[x_exp], Chen18_v_many['e'], 'kP', label='Chen18')
        plt.plot(Beladjine07_v26[x_exp], Beladjine07_v26['e'], 'k*', label='Beladjine07')
        #plt.plot(Zhou06_v_many[x_exp], Zhou06_v_many['e_mean'], 'k^', label='Zhou06')
        plt.plot(Rice95_v_many_medium[x_exp], Rice95_v_many_medium['e'], 'kD', label='Rice95 medium')
        plt.plot(Willetts89_v_many_medium[x_exp], Willetts89_v_many_medium['e'], 'kX', label='Willetts89 medium')
        #plt.plot(Rioual20_v_many[x_exp], Rioual20_v_many['e'], 'k.', label='Rioual20')
        #plt.plot(Gordon21_v_many[x_exp], Gordon21_v_many['e'], 'k1', label='Gordon21')
        #plt.plot(Gordon09_v_many[x_exp], Gordon09_v_many['e'], 'k2', label='Gordon09')
    elif impactor == 1 or impactor == 4:
        plt.plot(Rice95_v_many_coarse[x_exp], Rice95_v_many_coarse['e'], 'ks', label='Rice95 coarse')
        plt.plot(Willetts89_v_many_coarse[x_exp], Willetts89_v_many_coarse['e'], 'kH', label='Willetts89 coarse')
    elif impactor == 3:
        plt.plot(Rice95_v_many_fine[x_exp], Rice95_v_many_fine['e'], 'ko', label='Rice95 fine')
        plt.plot(Willetts89_v_many_fine[x_exp], Willetts89_v_many_fine['e'], 'k+', label='Willetts89 fine')

    plt.legend()

if output_ez:
    plt.figure(2)
    if variation_param == 0:
        plt.xlabel(r'$\hat{v}_1$')
        plt.ylabel(r'$\langle e_z \rangle_{\theta}, \langle E_z \rangle_{\theta}$')
        strlabel1 = r'$\langle e_z \rangle_{\theta}$'
        strlabel2 = r'$\langle E_z \rangle_{\theta}$'
    elif variation_param == 1:
        plt.xlabel(r'$\theta (degree)$')
        plt.ylabel(r'$\langle e_z \rangle_{\hat{v}}, \langle E_z \rangle_{\hat{v}}$')
        strlabel1 = r'$\langle e_z \rangle_{\hat{v}}$'
        strlabel2 = r'$\langle E_z \rangle_{\hat{v}}$'
    plt.plot(x0_array0, ez0_bar_list_0, 'r-', label=strlabel1 + ': d1, d2, d3')
    plt.plot(x0_array1, ez0_bar_list_1, 'b-', label=strlabel1 + ': d1, d2=d3')
    plt.plot(x0_array2, ez0_bar_list_2, 'g-', label=strlabel1 + ': d1=d2=d3')
    plt.plot(x_array0, ez_bar_list_0, 'r--', label=strlabel2 + ': d1, d2, d3')
    plt.plot(x_array1, ez_bar_list_1, 'b--', label=strlabel2 + ': d1, d2=d3')
    plt.plot(x_array2, ez_bar_list_2, 'g--', label=strlabel2 + ': d1=d2=d3')
    plt.legend()

if output_ex:
    plt.figure(3)
    if variation_param == 0:
        plt.xlabel(r'$\hat{v}_1$')
        plt.ylabel(r'$\langle e_x \rangle_{\theta}, \langle E_x \rangle_{\theta}$')
        strlabel1 = r'$\langle e_x \rangle_{\theta}$'
        strlabel2 = r'$\langle E_x \rangle_{\theta}$'
    elif variation_param == 1:
        plt.xlabel(r'$\theta (degree)$')
        plt.ylabel(r'$\langle e_x \rangle_{\hat{v}}, \langle E_x \rangle_{\hat{v}}$')
        strlabel1 = r'$\langle e_x \rangle_{\hat{v}}$'
        strlabel2 = r'$\langle E_x \rangle_{\hat{v}}$'
    plt.plot(x0_array0, ex0_bar_list_0, 'r-', label=strlabel1 + ': d1, d2, d3')
    plt.plot(x0_array1, ex0_bar_list_1, 'b-', label=strlabel1 + ': d1, d2=d3')
    plt.plot(x0_array2, ex0_bar_list_2, 'g-', label=strlabel1 + ': d1=d2=d3')
    plt.plot(x_array0, ex_bar_list_0, 'r--', label=strlabel2 + ': d1, d2, d3')
    plt.plot(x_array1, ex_bar_list_1, 'b--', label=strlabel2 + ': d1, d2=d3')
    plt.plot(x_array2, ex_bar_list_2, 'g--', label=strlabel2 + ': d1=d2=d3')
    plt.legend()

if output_rebound_ratio:
    plt.figure(4)
    if variation_param == 0:
        plt.xlabel(r'$\hat{v}_1$')
        plt.semilogx(x_array, rebound_ratio0_list_0, 'r-', label='e: d1, d2, d3')
        plt.semilogx(x_array, rebound_ratio0_list_1, 'b-', label='e: d1, d2=d3')
        plt.semilogx(x_array, rebound_ratio0_list_2, 'g-', label='e: d1=d2=d3')
        plt.semilogx(x_array, rebound_ratio_list_0, 'r--', label='ebar: d1, d2, d3')
        plt.semilogx(x_array, rebound_ratio_list_1, 'b--', label='ebar: d1, d2=d3')
        plt.semilogx(x_array, rebound_ratio_list_2, 'g--', label='ebar: d1=d2=d3')
    elif variation_param == 1:
        plt.xlabel(r'$\theta (degree)$')
        plt.plot(x_array, rebound_ratio0_list_0, 'r--', label='e: d1, d2, d3')
        plt.plot(x_array, rebound_ratio0_list_1, 'b--', label='e: d1, d2=d3')
        plt.plot(x_array, rebound_ratio0_list_2, 'g--', label='e: d1=d2=d3')
        plt.plot(x_array, rebound_ratio_list_0, 'r-', label='ebar: d1, d2, d3')
        plt.plot(x_array, rebound_ratio_list_1, 'b-', label='ebar: d1, d2=d3')
        plt.plot(x_array, rebound_ratio_list_2, 'g-', label='ebar: d1=d2=d3')
    plt.ylabel('Rebound Ratio')
    plt.legend()

if output_theta2_dist:
    plt.figure(5)
    theta20_list_0 = np.rad2deg(theta20_list_0)
    theta20_list_1 = np.rad2deg(theta20_list_1)
    theta20_list_2 = np.rad2deg(theta20_list_2)
    sns.histplot(theta20_list_0, bins=40, kde=True, color='red', stat='density', label='d1, d2, d3', element='step', alpha=0.1)
    sns.histplot(theta20_list_1, bins=40, kde=True, color='blue', stat='density', label='d1, d2=d3', element='step', alpha=0.1)
    sns.histplot(theta20_list_2, bins=40, kde=True, color='green', stat='density', label='d1=d2=d3', element='step', alpha=0.1)
    plt.xlabel(r"$\theta\'$", fontsize=16)
    plt.ylabel('pdf', fontsize=16)
    plt.legend()

if output_theta2:
    plt.figure(6)
    if variation_param == 0:
        plt.xlabel(r'$\hat{v}_1$')
        plt.ylabel(r'$\langle \theta\' \rangle_{\theta}, \langle \Theta\' \rangle_{\theta}$')
        strlabel1 = r'$\langle \theta\' \rangle_{\theta}$'
        strlabel2 = r'$\langle \Theta\' \rangle_{\theta}$'
    elif variation_param == 1:
        plt.xlabel(r'$\theta (degree)$')
        plt.ylabel(r'$\langle \theta\' \rangle_{\hat{v}}, \langle \Theta\' \rangle_{\hat{v}}$')
        strlabel1 = r'$\langle \theta\' \rangle_{\hat{v}}$'
        strlabel2 = r'$\langle \Theta\' \rangle_{\hat{v}}$'
    theta20_bar_list_0 = np.rad2deg(theta20_bar_list_0)
    theta20_bar_list_1 = np.rad2deg(theta20_bar_list_1)
    theta20_bar_list_2 = np.rad2deg(theta20_bar_list_2)
    theta2_bar_list_0 = np.rad2deg(theta2_bar_list_0)
    theta2_bar_list_1 = np.rad2deg(theta2_bar_list_1)
    theta2_bar_list_2 = np.rad2deg(theta2_bar_list_2)
    plt.plot(x0_array0, theta20_bar_list_0, 'r-', label=strlabel1 + ': d1, d2, d3')
    plt.plot(x0_array1, theta20_bar_list_1, 'b-', label=strlabel1 + ': d1, d2=d3')
    plt.plot(x0_array2, theta20_bar_list_2, 'g-', label=strlabel1 + ': d1=d2=d3')
    #plt.plot(x_array0, theta2_bar_list_0, 'r--', label=strlabel2 + ': d1, d2, d3')
    #plt.plot(x_array1, theta2_bar_list_1, 'b--', label=strlabel2 + ': d1, d2=d3')
    #plt.plot(x_array2, theta2_bar_list_2, 'g--', label=strlabel2 + ': d1=d2=d3')

    if impactor == 0 or impactor == 2:
        plt.plot(Chen18_v_many[x_exp], Chen18_v_many['ang_re'], 'kP', label='Chen et al. 2018')
        plt.plot(Beladjine07_v26[x_exp], Beladjine07_v26['ang_re'], 'k*', label='Beladjine et al. 2007')
        #plt.plot(Zhou06_v_many[x_exp], Zhou06_v_many['ang_re_mean'], 'k^', label='Zhou et al. 2006')
        plt.plot(Rice95_v_many_medium[x_exp], Rice95_v_many_medium['ang_re'], 'ko', label='Rice et al. 1995 (medium)')
        plt.plot(Willetts89_v_many_medium[x_exp], Willetts89_v_many_medium['ang_re'], 'kX', label='Willetts et al. 1989 (medium)')
    elif impactor == 1 or impactor == 4:
        plt.plot(Rice95_v_many_coarse[x_exp], Rice95_v_many_coarse['ang_re'], 'kD', label='Rice et al. 1995 (coarse)')
        plt.plot(Willetts89_v_many_coarse[x_exp], Willetts89_v_many_coarse['ang_re'], 'kH', label='Willetts et al. 1989 (coarse)')
    elif impactor == 3:
        plt.plot(Rice95_v_many_fine[x_exp], Rice95_v_many_fine['ang_re'], 'ks', label='Rice et al. 1995 (fine)')
        plt.plot(Willetts89_v_many_fine[x_exp], Willetts89_v_many_fine['ang_re'], 'k+', label='Willetts et al. 1989 (fine)')
    plt.legend()

if output_Nej:
    plt.figure(7)
    if variation_param == 0:
        plt.xlabel(r'$\hat{v}_1$')
        plt.ylabel(r'$\langle N_{ej} \rangle_{\theta}$')
    elif variation_param == 1:
        plt.xlabel(r'$\theta (degree)$')
        plt.ylabel(r'$\langle N_{ej} \rangle_{\hat{v}}$')
    plt.plot(x_array, Nej_bar_list_0, 'r-', label='d1, d2, d3')
    plt.plot(x_array, Nej_bar_list_1, 'b-', label='d1, d2=d3')
    plt.plot(x_array, Nej_bar_list_2, 'g-', label='d1=d2=d3')

    if impactor == 0 or impactor == 2:
        plt.plot(Rice95_v_many_medium[x_exp], Rice95_v_many_medium['Nej'], 'ko', label='Rice95 medium')
    if impactor == 1 or impactor == 4:
        plt.plot(Rice95_v_many_coarse[x_exp], Rice95_v_many_coarse['Nej'], 'kD', label='Rice95 coarse')
    if impactor == 3:
        plt.plot(Rice95_v_many_fine[x_exp], Rice95_v_many_fine['Nej'], 'ks', label='Rice95 fine')
    plt.legend()

if output_vn:
    plt.figure(8)
    if variation_param == 0:
        plt.xlabel(r'$\hat{v}_1$')
        plt.ylabel(r'$\langle v_n \rangle_{\theta}$')
    elif variation_param == 1:
        plt.xlabel(r'$\theta (degree)$')
        plt.ylabel(r'$\langle v_n \rangle_{\hat{v}}$')
    plt.plot(x_array, vn_bar_list_0, 'r-', label='d1, d2, d3')
    plt.plot(x_array, vn_bar_list_1, 'b-', label='d1, d2=d3')
    plt.plot(x_array, vn_bar_list_2, 'g-', label='d1=d2=d3')
    plt.legend()

plt.show()