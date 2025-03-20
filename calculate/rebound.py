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

def calculate_e_bar(d2_hat, theta1, alpha, beta):
    csc_theta1 = 1/np.sin(theta1)
    cot_theta1 = 1/np.tan(theta1)
    if 2*np.sin(theta1) < d2_hat:
        x_hat_min = csc_theta1 - d2_hat
    else:
        x_hat_min = cot_theta1 * np.sqrt(1 - (d2_hat/2)**2) - d2_hat/2
    x0 = x_hat_min
    x_hat_max = 1/np.sin(theta1)
    # no secondary collision
    attempt_num = 100
    delta_x_hat = (x_hat_max - x_hat_min)/attempt_num
    for n in range(1, attempt_num + 1):
        x1 = x_hat_min + n*delta_x_hat
        neq_LHS = 1 / (1 + beta/alpha)
        x1_sin_square = x1**2*np.sin(theta1)**2
        x1_sin_square = min(x1_sin_square, 1)
        x1_cos = x1*np.cos(theta1)
        neq_RHS = x1_sin_square - x1_cos*np.sqrt(1 - x1_sin_square)
        if neq_LHS < neq_RHS:
            break
    x0_sin_square = x0**2*np.sin(theta1)**2
    x0_sin_square = min(x0_sin_square, 1)
    x0_cos = x0*np.cos(theta1)
    eqx_x0 = -(alpha + beta)*(1 - x0_sin_square)**(3/2) + x0_cos*(-3*alpha + (alpha + beta)*x0_sin_square)
    eqx_x1 = -(alpha + beta)*(1 - x1_sin_square)**(3/2) + x1_cos*(-3*alpha + (alpha + beta)*x1_sin_square)
    e_x_bar = (eqx_x1 - eqx_x0)/(x1 - x0)/3
    x0_sin_32 = x0**3*np.sin(theta1)**2
    eqz_x0 = -alpha*x0 - 1/3*(alpha + beta)*(x0_sin_32 + cot_theta1*csc_theta1*(1 - x0_sin_square)**(3/2))
    x1_sin_32 = x1**3*np.sin(theta1)**2
    eqz_x1 = -alpha*x1 - 1/3*(alpha + beta)*(x1_sin_32 + cot_theta1*csc_theta1*(1 - x1_sin_square)**(3/2))
    e_z_bar = (eqz_x1 - eqz_x0)/(x1 - x0)*np.sin(theta1)
    e_bar = np.sqrt(e_x_bar**2 + e_z_bar**2)
    return e_bar

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

def calculate_survive(d1, d2, d3, g, v2_x, v2_z, e):
    x = d2*(d1+d2+d3) - d1*d3
    y = (d1 + d2)*(d2 + d3)
    psi = np.arccos(x/y)
    zb = (1.0 - np.sin(psi))*d2
    if v2_z**2/(2.0*g) > zb:
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
plot_type = 0 # 0: Ec
d1 = 0.0001
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
restitution_N = 0.78
restitution_T = -0.13
v1_hat_single = 100
theta1_single = np.pi/12
v1_hat_start = 0
v1_hat_end = 10
case_num = 100
theta1_start = np.pi/180
theta1_end = np.pi/2
theta2_num = 60
num_samples = 100
d_mid = (d_min + d_max) / 2
m_in = np.pi * d_mid**3 / 6 * rho
E_in = m_in * g * d_mid * 10
max_attempts = num_samples + 1
#------------------------------------------------------------------
# impact velocity
#x_array = np.linspace(v1_hat_start, v1_hat_end, case_num)
v1_hat = v1_hat_single
# impact angle
x_array = np.linspace(theta1_start, theta1_end, case_num)
#theta1 = theta1_single

mu, sigma = get_normal_params(normal_E, normal_D)
d1_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, num_samples)
d1 = np.percentile(d1_array, 50)
#d1 = np.percentile(d1_array, 90)

d2_list = []
rebound_ratio_list = []
rebound_ratio_list_1 = []
rebound_ratio_list_2 = []

for x in tqdm(x_array):
    #v1_hat = x
    theta1 = x
    attempts = 0
    rebound_num = 0
    rebound_num_1 = 0
    rebound_num_2 = 0
    d3_list = []
    while len(d3_list) < num_samples and attempts < max_attempts:
        attempts += 1
        if distribution == 0:
            d3 = np.random.uniform(d_min, d_max)
            d2_min = (-(d1+d3)+np.sqrt((d1+d3)**2+4*d1*d3))/2
            d2_min = max(d_min, d2_min)
            d2 = np.random.uniform(d2_min, d_max)
            d2_list.append(d2)
            d3_list.append(d3)
        elif distribution == 1:
            mu, sigma = get_normal_params(normal_E, normal_D)
            #d1 = generate_truncated_lognormal(mu, sigma, d_min, d_max, 1)[0]
            d3 = generate_truncated_lognormal(mu, sigma, d_min, d_max, 1)[0]
            d2_1 = d3
            d2_2 = d1
            d3_2 = d1
            d2_min = (-(d1+d3)+np.sqrt((d1+d3)**2+4*d1*d3))/2
            d2_min = max(d_min, d2_min)
            d2 = generate_truncated_lognormal(mu, sigma, d2_min, d_max, 1)[0]
            d2_list.append(d2)
            d3_list.append(d3)
        elif distribution == 2:
            r = np.random.uniform(0.5,1.5,size=2)
            rr = np.floor(r)
            d_bin = d_max - d_min
            d2 = d_min + rr[0]*d_bin
            d3 = d_min + rr[1]*d_bin
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

        # normalize d1, d2
        d = (d1 + d2) / 2
        d_1 = (d1 + d2_1) / 2
        d_2 = (d1 + d2_2) / 2
        d1_hat = d1 / d
        d1_hat_1 = d1 / d_1
        d1_hat_2 = d1 / d_2
        d2_hat = d2 / d
        d2_hat_1 = d2 / d_1
        d2_hat_2 = d2 / d_2
        v1 = v1_hat * np.sqrt(g * d1)
        # restitution coefficient
        mu = restitution_N*d1_hat**3/(d1_hat**3 + restitution_N*d2_hat**3)
        mu_1 = restitution_N*d1_hat_1**3/(d1_hat_1**3 + restitution_N*d2_hat_1**3)
        mu_2 = restitution_N*d1_hat_2**3/(d1_hat_2**3 + restitution_N*d2_hat_2**3)
        alpha = (1 + restitution_N)/(1 + mu) - 1
        alpha_1 = (1 + restitution_N)/(1 + mu_1) - 1
        alpha_2 = (1 + restitution_N)/(1 + mu_2) - 1
        beta = 1 - (2/7)*(1 - restitution_T)/(1 + mu)
        beta_1 = 1 - (2/7)*(1 - restitution_T)/(1 + mu_1)
        beta_2 = 1 - (2/7)*(1 - restitution_T)/(1 + mu_2)
        e = calculate_e_bar(d2_hat, theta1, alpha, beta)
        e_1 = calculate_e_bar(d2_hat_1, theta1, alpha_1, beta_1)
        e_2 = calculate_e_bar(d2_hat_2, theta1, alpha_2, beta_2)
        v2 = v1*e
        v2_1 = v1*e_1
        v2_2 = v1*e_2
        # rebound angle
        theta2 = calculate_rebound_angle(d2_hat, theta1, alpha, beta)
        theta2_1 = calculate_rebound_angle(d2_hat_1, theta1, alpha_1, beta_1)
        theta2_2 = calculate_rebound_angle(d2_hat_2, theta1, alpha_2, beta_2)
        v2_z = v2*np.sin(theta2)
        v2_z_1 = v2_1*np.sin(theta2_1)
        v2_z_2 = v2_2*np.sin(theta2_2)
        v2_x = np.sqrt(v2**2 - v2_z**2)
        v2_x_1 = np.sqrt(v2_1**2 - v2_z_1**2)
        v2_x_2 = np.sqrt(v2_2**2 - v2_z_2**2)
        # survive or not
        is_survive = calculate_survive(d1, d2, d3, g, v2_x, v2_z, e)
        is_survive_1 = calculate_survive(d1, d2_1, d3, g, v2_x_1, v2_z_1, e_1)
        is_survive_2 = calculate_survive(d1, d2_2, d3, g, v2_x_2, v2_z_2, e_2)
        if is_survive:
            rebound_num += 1
        if is_survive_1:
            rebound_num_1 += 1
        if is_survive_2:
            rebound_num_2 += 1
    rebound_ratio = rebound_num/num_samples
    rebound_ratio_1 = rebound_num_1/num_samples
    rebound_ratio_2 = rebound_num_2/num_samples
    rebound_ratio_list.append(rebound_ratio)
    rebound_ratio_list_1.append(rebound_ratio_1)
    rebound_ratio_list_2.append(rebound_ratio_2)
plt.plot(x_array, rebound_ratio_list, 'ro-', label='d1, d2, d3')
plt.plot(x_array, rebound_ratio_list_1, 'bo-', label='d1, d2=d3')
plt.plot(x_array, rebound_ratio_list_2, 'go-', label='d1=d2=d3')
plt.xlabel('Impact Angle (rad)')
plt.ylabel('Rebound Ratio')
plt.legend()
plt.show()