import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import truncnorm
from tqdm import tqdm

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
    theta_c = np.arccos(cos1) + np.arccos(cos2) - np.pi/2
    l_r = (1 + d23**2 - d13**2)/(2*d23)
    if theta1 <= theta_c:
        x_min = d13/np.sin(theta1) - d23
    else:
        x_min = np.sqrt(1 - l_r**2)/np.tan(theta1) - l_r
    return x_min

def calculate_x_max(alpha, beta, theta1):
    x_max_try = 1/np.sin(theta1)
    x_min_try = 1/np.tan(theta1)
    attempt_num = 100
    delta_x = (x_max_try - x_min_try)/attempt_num
    for n in range(1, attempt_num + 1):
        x = x_min_try + n*delta_x
        x_sin_square = x**2*np.sin(theta1)**2
        x_cos = x*np.cos(theta1)
        neq_LHS = 1/(1 + beta/alpha)
        neq_RHS = x_sin_square - x_cos*np.sqrt(1 - x_sin_square)
        if neq_LHS < neq_RHS:
            break
    x_max = x
    return x_max


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
    e_bar = np.sqrt(e_vx**2 + e_vz**2)
    return e_bar, e_vz

def calculate_e_bar_simple(d2_hat, theta1, alpha, beta):
    e_bar = beta - (beta**2 - alpha**2)*d2_hat*theta1/(2*beta)
    ez_bar = -beta + (2/3)*(alpha + beta)*np.sqrt(2*d2_hat/theta1)
    e_z_bar = ez_bar*np.sin(theta1)
    return e_bar, e_z_bar

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

def calculate_survive(d1, d2, d3, theta1, g, v2_x, v2_z, e):
    x = d2*(d1+d2+d3) - d1*d3
    y = (d1 + d2)*(d2 + d3)
    psi = np.arccos(x/y)
    psi = min(psi, theta1)
    zb = (1.0 - np.sin(psi))*0.5*(d1 + d2)
    z_final = v2_z**2/(2.0*g)
    if z_final > zb:
        res_x = (np.sqrt((v2_z/g)**2 - (2.0*zb)/g) + v2_z/g)*v2_x - 0.5*d2
    else:
        res_x = -1
    if res_x <= 0.0 or e <= 0.0:
        is_survive = False
    else:
        is_survive = True
    return is_survive

#-------------------------------------------------------------------
distribution = 1 # 0:uniform, 1:lognormal, 2:bidisperse, 3:polydisperse
variation_param = 0 # 0: v1, 1: theta1
d_min = 1e-4
d_max = 10e-4
normal_E = 4e-4
normal_D = 2e-4
mu1 = (d_min + d_max) * 0.3  # 第一个峰靠左
mu2 = (d_min + d_max) * 0.7  # 第二个峰靠右
sigma1 = (d_max - d_min) * 0.1
sigma2 = (d_max - d_min) * 0.1
weight1 = 0.5 # 第一个峰的权重
rho = 2650
g = 9.8*(1 - 1.263/rho)
epsilon = 0.78
nu = -0.13
v1_hat_single = 5
theta1_single = np.pi/12
v1_hat_start = 0
v1_hat_end = 100
case_num = 100
theta1_start = np.pi/180
theta1_end = np.pi/1.9
theta2_num = 60
num_samples = 200
#------------------------------------------------------------------
# impact velocity
if variation_param == 0:
    x_array = np.linspace(v1_hat_start, v1_hat_end, case_num)
else:
    v1_hat = v1_hat_single
# impact angle
if variation_param == 1:
    x_array = np.linspace(theta1_start, theta1_end, case_num)
else:
    theta1 = theta1_single
if distribution == 0:
    d1_array = np.random.uniform(d_min, d_max, size=num_samples)
    d_mid = np.percentile(d1_array, 50)
    #d_mid = np.percentile(d1_array, 90)
elif distribution == 1:
    mu, sigma = get_normal_params(normal_E, normal_D)
    d1_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, num_samples)
    d_mid = np.percentile(d1_array, 50)
    #d_mid = np.percentile(d1_array, 90)
elif distribution == 2:
    d1_mid = d_min
    #d_mid = (d_min + d_max) * 0.5
d1 = d_mid
m_in = np.pi * d1**3 / 6 * rho
E_in = m_in * g * d1 * 10

d2_list = []
rebound_ratio_list = []
rebound_ratio_list_1 = []
rebound_ratio_list_2 = []

for x in tqdm(x_array):
    if variation_param == 0:
        v1_hat = x
    elif variation_param == 1:
        theta1 = x
    rebound_num_array = [0, 0, 0]
    d3_list = []
    while len(d3_list) < num_samples:
        if distribution == 0:
            d_array = np.random.uniform(d_min, d_max, size=3)
            #d1 = d_array[0]
            d2 = d_array[1]
            d3 = d_array[2]
            d2_list.append(d2)
            d3_list.append(d3)
        elif distribution == 1:
            mu, sigma = get_normal_params(normal_E, normal_D)
            d_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, 3)
            #d1 = d_array[0]
            d2 = d_array[1]
            d3 = d_array[1]
            d2_list.append(d2)
            d3_list.append(d3)
        elif distribution == 2:
            r = np.random.uniform(0.5,1.5,size=3)
            rr = np.floor(r)
            d_bin = d_max - d_min
            #d1 = d_min + rr[0]*d_bin
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
        for i in range(3):
            if i == 1:
                d2 = d2
                d3 = d2
            elif i == 2:
                d2 = d1
                d3 = d1

            v1 = v1_hat * np.sqrt(g * d1)
            d = (d1 + d2)/2
            d1_hat = d1/d
            d2_hat = d2/d
            d3_hat = d3/d
            # restitution coefficient
            mu = epsilon*d1_hat**3/(d1_hat**3 + epsilon*d2_hat**3)
            alpha = (1 + epsilon)/(1 + mu) - 1
            beta = 1 - (2/7)*(1 - nu)/(1 + mu)
            x_min = calculate_x_min(d1_hat, d2_hat, d3_hat, theta1)
            x_max = calculate_x_max(alpha, beta, theta1)
            e, evz = calculate_e_bar(alpha, beta, x_min, x_max, theta1)
            v2 = v1*e
            v2_z = v1*evz
            v2_x = np.sqrt(v2**2 - v2_z**2)
            # TODO: check function calculate_survive
            # survive or not
            is_survive = calculate_survive(d1, d2, d3, theta1, g, v2_x, v2_z, e)
            if is_survive:
                rebound_num_array[i] += 1
    rebound_ratio = rebound_num/num_samples
    rebound_ratio_1 = rebound_num_1/num_samples
    rebound_ratio_2 = rebound_num_2/num_samples
    rebound_ratio_list.append(rebound_ratio)
    rebound_ratio_list_1.append(rebound_ratio_1)
    rebound_ratio_list_2.append(rebound_ratio_2)
plt.plot(x_array, rebound_ratio_list, 'ro-', label='d1, d2, d3')
plt.plot(x_array, rebound_ratio_list_1, 'bo-', label='d1, d2=d3')
plt.plot(x_array, rebound_ratio_list_2, 'go-', label='d1=d2=d3')
if variation_param == 0:
    plt.xlabel('v1_hat')
elif variation_param == 1:
    plt.xlabel('Impact Angle (rad)')
plt.ylabel('Rebound Ratio')
plt.legend()
plt.show()