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

def calculate_e(alpha, beta, x, theta1):
    # 任意x位置上的e_vx可能出现小于0的情况，但其在xmin和xmax范围内的积分是正的
    e_vx = -alpha*np.cos(theta1) + (alpha + beta)*x**2*np.sin(theta1)**2*np.cos(theta1) + (alpha + beta)*x*np.sin(theta1)**2*np.sqrt(1 - x**2*np.sin(theta1)**2)
    e_vz = alpha*np.sin(theta1) - (alpha + beta)*x**2*np.sin(theta1)**3 + (alpha + beta)*x*np.sin(theta1)*np.cos(theta1)*np.sqrt(1 - x**2*np.sin(theta1)**2)
    e = np.sqrt(e_vx**2 + e_vz**2)
    return e, e_vx, e_vz

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
distribution = 1 # 0:uniform, 1:lognormal, 2:bidisperse, 3:polydisperse
variation_param = 1 # 0: v1, 1: theta1
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
v1_hat_start = 1
v1_hat_end = 10
case_num = 100
theta1_start = np.pi/180
theta1_end = np.pi/2.1
theta2_num = 60
num_samples = 1000
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
    #d_mid = np.percentile(d1_array, 50)
    d_mid = np.percentile(d1_array, 90)
elif distribution == 2:
    #d_mid = d_min
    d_mid = (d_min + d_max) * 0.5
d1 = d_mid
m_in = np.pi * d1**3 / 6 * rho
E_in = m_in * g * d1 * 10

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
            d3 = d_array[2]
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

            v1 = v1_hat * np.sqrt(g * d1)
            d = (d1 + d2)/2
            d1_hat = d1/d
            d2_hat = d2/d
            d3_hat = d3/d
            # restitution coefficient
            mu = epsilon*d1_hat**3/(d1_hat**3 + epsilon*d2_hat**3)
            alpha = (1 + epsilon)/(1 + mu) - 1
            beta = 1 - (2/7)*(1 - nu)/(1 + mu)
            x_min, psi1, psi2 = calculate_x_min(d1_hat, d2_hat, d3_hat, theta1)
            x_max = calculate_x_max(alpha, beta, theta1)
            if x_min < x_max:
                e, evx, evz = calculate_e_bar(alpha, beta, x_min, x_max, theta1)
                x0_hat = np.random.uniform(x_min, x_max)
                e0, evx0, evz0 = calculate_e(alpha, beta, x0_hat, theta1)
            else:
                d_temp = d2
                d2 = d3
                d3 = d_temp
                #d3 = d2
                d = (d1 + d2)/2
                d1_hat = d1/d
                d2_hat = d2/d
                d3_hat = d3/d
                mu = epsilon*d1_hat**3/(d1_hat**3 + epsilon*d2_hat**3)
                alpha = (1 + epsilon)/(1 + mu) - 1
                beta = 1 - (2/7)*(1 - nu)/(1 + mu)
                x_min, psi1, psi2 = calculate_x_min(d1_hat, d2_hat, d3_hat, theta1)
                x_max = calculate_x_max(alpha, beta, theta1)
                e, evx, evz = calculate_e_bar(alpha, beta, x_min, x_max, theta1)
                x0_hat = np.random.uniform(x_min, x_max)
                e0, evx0, evz0 = calculate_e(alpha, beta, x0_hat, theta1)
            ez = evz/np.sin(theta1)
            ex = evx/np.cos(theta1)
            ez0 = evz0/np.sin(theta1)
            ex0 = evx0/np.cos(theta1)
            theta2 = np.arcsin(evz/e)
            theta20 = np.arcsin(evz0/e0)
            if i == 0:
                e_list_0.append(e)
                e0_list_0.append(e0)
                ez_list_0.append(ez)
                ex_list_0.append(ex)
                ez0_list_0.append(ez0)
                ex0_list_0.append(ex0)
                theta2_list_0.append(theta2)
                theta20_list_0.append(theta20)
            elif i == 1:
                e_list_1.append(e)
                e0_list_1.append(e0)
                ez_list_1.append(ez)
                ex_list_1.append(ex)
                ez0_list_1.append(ez0)
                ex0_list_1.append(ex0)
                theta2_list_1.append(theta2)
                theta20_list_1.append(theta20)
            elif i == 2:
                e_list_2.append(e)
                e0_list_2.append(e0)
                ez_list_2.append(ez)
                ex_list_2.append(ex)
                ez0_list_2.append(ez0)
                ex0_list_2.append(ex0)
                theta2_list_2.append(theta2)
                theta20_list_2.append(theta20)
            v2 = e*v1
            v20 = e0*v1
            v2_z = v1*evz
            v2_z0 = v1*evz0
            v2_x = v1*evx
            v2_x0 = v1*evx0
            is_survive = calculate_survive(d1, d2, d3, psi1, psi2, g, v2_x, v2_z, e)
            is_survive0 = calculate_survive(d1, d2, d3, psi1, psi2, g, v2_x0, v2_z0, e0)
            #is_survive = calculate_survive1(d1, d2, v1, theta1, x0_hat, g, v2_x, v2_z, e0)
            if is_survive:
                rebound_num_array[i] += 1
            if is_survive0:
                rebound_num_array0[i] += 1
    e_bar_0 = np.mean(e_list_0)
    e_bar_1 = np.mean(e_list_1)
    e_bar_2 = np.mean(e_list_2)
    e0_bar_0 = np.mean(e0_list_0)
    e0_bar_1 = np.mean(e0_list_1)
    e0_bar_2 = np.mean(e0_list_2)
    ez_bar_0 = np.mean(ez_list_0)
    ez_bar_1 = np.mean(ez_list_1)
    ez_bar_2 = np.mean(ez_list_2)
    ez0_bar_0 = np.mean(ez0_list_0)
    ez0_bar_1 = np.mean(ez0_list_1)
    ez0_bar_2 = np.mean(ez0_list_2)
    ex_bar_0 = np.mean(ex_list_0)
    ex_bar_1 = np.mean(ex_list_1)
    ex_bar_2 = np.mean(ex_list_2)
    ex0_bar_0 = np.mean(ex0_list_0)
    ex0_bar_1 = np.mean(ex0_list_1)
    ex0_bar_2 = np.mean(ex0_list_2)
    e_bar_list_0.append(e_bar_0)
    e_bar_list_1.append(e_bar_1)
    e_bar_list_2.append(e_bar_2)
    e0_bar_list_0.append(e0_bar_0)
    e0_bar_list_1.append(e0_bar_1)
    e0_bar_list_2.append(e0_bar_2)
    ez_bar_list_0.append(ez_bar_0)
    ez_bar_list_1.append(ez_bar_1)
    ez_bar_list_2.append(ez_bar_2)
    ex_bar_list_0.append(ex_bar_0)
    ex_bar_list_1.append(ex_bar_1)
    ex_bar_list_2.append(ex_bar_2)
    ez0_bar_list_0.append(ez0_bar_0)
    ez0_bar_list_1.append(ez0_bar_1)
    ez0_bar_list_2.append(ez0_bar_2)
    ex0_bar_list_0.append(ex0_bar_0)
    ex0_bar_list_1.append(ex0_bar_1)
    ex0_bar_list_2.append(ex0_bar_2)
    rebound_ratio_array = [rebound_num_array[i]/num_samples for i in range(3)]
    rebound_ratio_array0 = [rebound_num_array0[i]/num_samples for i in range(3)]
    rebound_ratio_list_0.append(rebound_ratio_array[0])
    rebound_ratio_list_1.append(rebound_ratio_array[1])
    rebound_ratio_list_2.append(rebound_ratio_array[2])
    rebound_ratio0_list_0.append(rebound_ratio_array0[0])
    rebound_ratio0_list_1.append(rebound_ratio_array0[1])
    rebound_ratio0_list_2.append(rebound_ratio_array0[2])

plt.figure(1)
if variation_param == 0:
    plt.xlabel(r'$\hat{v}_1$')
    plt.ylabel(r'$\langle e \rangle_{\theta}, \langle E \rangle_{\theta}$')
    strlabel1 = r'$\langle e \rangle_{\theta}$'
    strlabel2 = r'$\langle E \rangle_{\theta}$'
elif variation_param == 1:
    x_array = np.rad2deg(x_array)
    plt.xlabel(r'$\theta (degree)$')
    plt.ylabel(r'$\langle e \rangle_{\hat{v}}, \langle E \rangle_{\hat{v}}$')
    strlabel1 = r'$\langle e \rangle_{\hat{v}}$'
    strlabel2 = r'$\langle E \rangle_{\hat{v}}$'
plt.plot(x_array, e0_bar_list_0, 'r-', label=strlabel1 + ': d1, d2, d3')
plt.plot(x_array, e_bar_list_0, 'r--', label=strlabel2 + ': d1, d2, d3')
plt.plot(x_array, e0_bar_list_1, 'b-', label=strlabel1 + ': d1, d2=d3')
plt.plot(x_array, e_bar_list_1, 'b--', label=strlabel2 + ': d1, d2=d3')
plt.plot(x_array, e0_bar_list_2, 'g-', label=strlabel1 + ': d1=d2=d3')
plt.plot(x_array, e_bar_list_2, 'g--', label=strlabel2 + ': d1=d2=d3')
plt.legend()

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
plt.plot(x_array, ez0_bar_list_0, 'r-', label=strlabel1 + ': d1, d2, d3')
plt.plot(x_array, ez_bar_list_0, 'r--', label=strlabel2 + ': d1, d2, d3')
plt.plot(x_array, ez0_bar_list_1, 'b-', label=strlabel1 + ': d1, d2=d3')
plt.plot(x_array, ez_bar_list_1, 'b--', label=strlabel2 + ': d1, d2=d3')
plt.plot(x_array, ez0_bar_list_2, 'g-', label=strlabel1 + ': d1=d2=d3')
plt.plot(x_array, ez_bar_list_2, 'g--', label=strlabel2 + ': d1=d2=d3')
plt.legend()

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
plt.plot(x_array, ex0_bar_list_0, 'r-', label=strlabel1 + ': d1, d2, d3')
plt.plot(x_array, ex_bar_list_0, 'r--', label=strlabel2 + ': d1, d2, d3')
plt.plot(x_array, ex0_bar_list_1, 'b-', label=strlabel1 + ': d1, d2=d3')
plt.plot(x_array, ex_bar_list_1, 'b--', label=strlabel2 + ': d1, d2=d3')
plt.plot(x_array, ex0_bar_list_2, 'g-', label=strlabel1 + ': d1=d2=d3')
plt.plot(x_array, ex_bar_list_2, 'g--', label=strlabel2 + ': d1=d2=d3')
plt.legend()

plt.figure(4)
if variation_param == 0:
    plt.xlabel(r'$\hat{v}_1$')
    plt.semilogx(x_array, rebound_ratio0_list_0, 'r-', label='e: d1, d2, d3')
    plt.semilogx(x_array, rebound_ratio_list_0, 'r--', label='ebar: d1, d2, d3')
    plt.semilogx(x_array, rebound_ratio0_list_1, 'b-', label='e: d1, d2=d3')
    plt.semilogx(x_array, rebound_ratio_list_1, 'b--', label='ebar: d1, d2=d3')
    plt.semilogx(x_array, rebound_ratio0_list_2, 'g-', label='e: d1=d2=d3')
    plt.semilogx(x_array, rebound_ratio_list_2, 'g--', label='ebar: d1=d2=d3')
elif variation_param == 1:
    plt.xlabel(r'$\theta (degree)$')
    plt.plot(x_array, rebound_ratio_list_0, 'r-', label='ebar: d1, d2, d3')
    plt.plot(x_array, rebound_ratio0_list_0, 'r--', label='e: d1, d2, d3')
    plt.plot(x_array, rebound_ratio_list_1, 'b-', label='ebar: d1, d2=d3')
    plt.plot(x_array, rebound_ratio0_list_1, 'b--', label='e: d1, d2=d3')
    plt.plot(x_array, rebound_ratio_list_2, 'g-', label='ebar: d1=d2=d3')
    plt.plot(x_array, rebound_ratio0_list_2, 'g--', label='e: d1=d2=d3')
plt.ylabel('Rebound Ratio')
plt.legend()

plt.show()