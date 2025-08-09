import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.stats import truncnorm
from scipy.stats import norm
from tqdm import tqdm
from scipy.special import erfc
from scipy.ndimage import gaussian_filter
from scipy.integrate import quad
from sklearn.cross_decomposition import CCA

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
    alpha = np.arccos(cos0)
    #theta_c = gamma + delta - np.pi/2
    theta_c = np.pi/2 - alpha
    l_r = (1 + d23**2 - d13**2)/(2*d23)
    if theta1 <= theta_c:
        x_min = d13/np.sin(theta1) - d23
    else:
        x_min = np.sqrt(1 - l_r**2)/np.tan(theta1) - l_r
    return x_min

def calculate_x_min_3D(d1, d2, d3, d4, theta1):
    d13 = (d1 + d3)/2
    d23 = (d2 + d3)/2
    d14 = (d1 + d4)/2
    d34 = (d3 + d4)/2
    kmax = (d13**2 + d34**2 - d14**2)/(2*d13*d34)
    k = np.random.uniform(0, kmax)
    d13 = d13*np.sqrt(1.0 - k**2)
    d23 = d23*np.sqrt(1.0 - (k*d13/d23)**2)
    cos0 = (d13**2 + d23**2 - 1)/(2*d13*d23)
    alpha = np.arccos(cos0)
    theta_c = np.pi/2 - alpha
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

#-------------------------------------------------------------------
output_e = True # output e
output_theta2 = True # output rebound angle
output_vn = False # output eject velocity
output_Nej = False # output eject number

distribution = 1 # 0:uniform, 1:lognormal, 2:bidisperse, 3:polydisperse, 4:normal
shallow = False # shallow impact
simplify = False # first order approximation
lognormal_param = True # lognormal distribution parameters
Three_D = False # 3D bed

no_calculate = False # do not calculate, just output figure

d_min = 1.5e-4
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
gamma_ej = 0.049 # the fraction of remaining energy to eject particles
rho = 2650
g = 9.8*(1 - 1.263/rho)
num_samples = 200
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
    'd_in': np.mean([coarse_d1, coarse_d1, coarse_d1, coarse_d1]),
    'v_in': np.mean([2.7323, 2.7013, 2.7070, 2.7979]),
    'ang_in': np.mean([13.94, 14.75, 14.73, 15.04]),
    'ang_re': np.mean([20.95, 23.03, 22.55, 25.63]),
    'e': np.mean([0.58, 0.56, 0.56, 0.58]),
    'Nej': np.mean([5.7, 5.09, 5.12, 4.64]),
    'v_ej': np.mean([0.243, 0.2422, 0.2522, 0.2637])
}
Rice95_v_many_medium= {
    'd_in': np.mean([medium_d1, medium_d1, medium_d1, medium_d1]),
    'v_in': np.mean([3.1649, 3.3638, 3.3604, 3.2963]),
    'ang_in': np.mean([11.82, 11.47, 11.53, 11.36]),
    'ang_re': np.mean([29.95, 30.31, 31.06, 29.56]),
    'e': np.mean([0.57, 0.54, 0.59, 0.58]),
    'Nej': np.mean([2.68, 3.43, 2.67, 2.76]),
    'v_ej': np.mean([0.2544, 0.252, 0.2607, 0.295])
}
Rice95_v_many_fine= {
    'd_in': np.mean([fine_d1, fine_d1, fine_d1]),
    'v_in': np.mean([3.8493, 3.7801, 3.7411]),
    'ang_in': np.mean([10.85, 10.24, 10.46]),
    'ang_re': np.mean([44.63, 38.03, 37.85]),
    'e': np.mean([0.52, 0.56, 0.58]),
    'Nej': np.mean([1.86, 1.71, 1.58]),
    'v_ej': np.mean([0.2757, 0.2566, 0.2773])
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

d1_array = [coarse_d1, medium_d1, fine_d1]
d1_array = d1_array + Willetts89_v_many_coarse['d_in'] + Willetts89_v_many_medium['d_in'] + Willetts89_v_many_fine['d_in']
v1_array = [Rice95_v_many_coarse['v_in'], Rice95_v_many_medium['v_in'], Rice95_v_many_fine['v_in']]
v1_array = v1_array + Willetts89_v_many_coarse['v_in'] + Willetts89_v_many_medium['v_in'] + Willetts89_v_many_fine['v_in']
theta1_array = [Rice95_v_many_coarse['ang_in'], Rice95_v_many_medium['ang_in'], Rice95_v_many_fine['ang_in']]
theta1_array = theta1_array + Willetts89_v_many_coarse['ang_in'] + Willetts89_v_many_medium['ang_in'] + Willetts89_v_many_fine['ang_in']
iteration_num = len(d1_array)

d_epsilon = 0.1
d_nu = 0.1
epsilon_list = np.arange(0.3, 1.2 + 1.1*d_epsilon, d_epsilon).tolist()
nu_list = np.arange(-1.2, -0.3 + 1.1*d_nu, d_nu).tolist()
rmse_list = []

if no_calculate:
    data = np.loadtxt('eject_compair_rmse.txt', delimiter=',')
    read_x = data[:, 0]
    read_y = data[:, 1]
    read_z = data[:, 2]
    mesh_x = np.unique(read_x)
    mesh_y = np.unique(read_y)
    mesh_z = read_z.reshape(len(mesh_y), len(mesh_x))
    mesh_z_smooth = gaussian_filter(mesh_z, sigma=3.0)
    read_z_smooth = mesh_z_smooth.flatten()
    plt.contourf(mesh_x, mesh_y, mesh_z_smooth, cmap='viridis', levels=100)
    cbar = plt.colorbar(ticks=np.linspace(0, 5, 11))
    cbar.set_label(r'$\Delta^{\mathrm{norm}}_{\mathrm{ave}}$', rotation=0, labelpad=20, fontsize=14)
    cbar.ax.tick_params(labelsize=10)
    plt.xlabel(r'$\epsilon$', fontsize=14)
    plt.ylabel(r'$\nu$', fontsize=14)
    plt.tick_params(labelsize=10)
    plt.show()
    min_indices = np.argpartition(read_z_smooth, 10)[:10]  # Get indices of the 10 smallest RMSE values
    for min_index in min_indices:
        min_x = read_x[min_index]
        min_y = read_y[min_index]
        min_z = read_z[min_index]
        print(f'Minimum RMSE: {min_z:.4f} at Epsilon: {min_x:.4f}, Nu: {min_y:.4f}')
    sys.exit(0)

os.remove('epsilon_nu_CCA_weights.txt') if os.path.exists('epsilon_nu_CCA_weights.txt') else None
os.remove('epsilon_nu_corr.txt') if os.path.exists('epsilon_nu_corr.txt') else None
for epsilon in tqdm(epsilon_list):
    for nu in tqdm(nu_list):
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
        for n in range(iteration_num):
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
                    d_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, 4)
                    d2 = d_array[1]
                    d3 = d_array[2]
                    d4 = d_array[3]
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

                d = (d1 + d2)/2
                d1_hat = d1/d
                d2_hat = d2/d
                d3_hat = d3/d
                d4_hat = d4/d
                # restitution coefficient
                mu_re = epsilon*d1_hat**3/(d1_hat**3 + epsilon*d2_hat**3)
                alpha = (1 + epsilon)/(1 + mu_re) - 1
                beta = 1 - (2/7)*(1 - nu)/(1 + mu_re)
                if Three_D:
                    x_min = calculate_x_min_3D(d1_hat, d2_hat, d3_hat, d4_hat, theta1)
                else:
                    x_min = calculate_x_min(d1_hat, d2_hat, d3_hat, theta1)

                x_max = calculate_x_max(alpha, beta, theta1)
                if x_min < x_max:
                    dont_change_flag = True
                else:
                    dont_change_flag = False

                if dont_change_flag:
                    #x0_hat = np.random.uniform(x_min, x_max)
                    #e0, evx0, evz0 = calculate_e(alpha, beta, x0_hat, theta1)
                    theta20, error = quad(rebound_angle_bar, x_min, x_max, args=(alpha, beta, theta1))
                    theta20 = theta20 / (x_max - x_min)
                    e0, error = quad(rebound_res_bar, x_min, x_max, args=(alpha, beta, theta1))
                    e0 = e0 / (x_max - x_min)
                    #e0, evx0, evz0 = calculate_e_bar(alpha, beta, x_min, x_max, theta1)
                else:
                    d_temp = d2
                    d2 = d3
                    d3 = d_temp
                    d = (d1 + d2)/2
                    d1_hat = d1/d
                    d2_hat = d2/d
                    d3_hat = d3/d
                    d4_hat = d4/d
                    mu_re = epsilon*d1_hat**3/(d1_hat**3 + epsilon*d2_hat**3)
                    alpha = (1 + epsilon)/(1 + mu_re) - 1
                    beta = 1 - (2/7)*(1 - nu)/(1 + mu_re)
                    if Three_D:
                        x_min = calculate_x_min_3D(d1_hat, d2_hat, d3_hat, d4_hat, theta1)
                    else:
                        x_min = calculate_x_min(d1_hat, d2_hat, d3_hat, theta1)
                    x_max = calculate_x_max(alpha, beta, theta1)
                    #x0_hat = np.random.uniform(x_min, x_max)
                    #e0, evx0, evz0 = calculate_e(alpha, beta, x0_hat, theta1)
                    theta20, error = quad(rebound_angle_bar, x_min, x_max, args=(alpha, beta, theta1))
                    theta20 = theta20 / (x_max - x_min)
                    e0, error = quad(rebound_res_bar, x_min, x_max, args=(alpha, beta, theta1))
                    e0 = e0 / (x_max - x_min)
                    #e0, evx0, evz0 = calculate_e_bar(alpha, beta, x_min, x_max, theta1)

                if not lognormal_param:
                    mu, sigma = get_normal_params(normal_E, normal_D)
                #dn = generate_truncated_lognormal(mu, sigma, d_min, d_max, 1)[0]
                m1 = rho*(d1**3*np.pi/6)
                #mn = rho*(dn**3*np.pi/6)
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
                vn_list_0.append(vn)
                Nej_list_0.append(Nej)

                #ez0 = evz0 #/np.sin(theta1)
                #ex0 = evx0 #/np.cos(theta1)
                #theta20 = np.arctan(evz0/evx0)
                #if theta20 < 0:
                #    theta20 += np.pi
                e0_list_0.append(e0)
                #ez0_list_0.append(ez0)
                #ex0_list_0.append(ex0)
                theta20_list_0.append(theta20)
                v20 = e0*v1
                #v2_z0 = v1*evz0
                #v2_x0 = v1*evx0
            """d1, d2, d3"""
            e0_bar_0 = np.mean(e0_list_0)
            #ez0_bar_0 = np.mean(ez0_list_0)
            #ex0_bar_0 = np.mean(ex0_list_0)
            theta20_bar_0 = np.mean(theta20_list_0)
            e0_bar_list_0.append(e0_bar_0)
            #ez0_bar_list_0.append(ez0_bar_0)
            #ex0_bar_list_0.append(ex0_bar_0)
            theta20_bar_list_0.append(theta20_bar_0)
            """ejection velocity"""
            vn_bar_0 = np.mean(vn_list_0)
            vn_bar_list_0.append(vn_bar_0)
            """ejection number"""
            Nej_bar_0 = np.mean(Nej_list_0)
            Nej_bar_list_0.append(Nej_bar_0)

        len_W89_c = len(Willetts89_v_many_coarse['e'])
        len_W89_m = len(Willetts89_v_many_medium['e'])
        len_W89_f = len(Willetts89_v_many_fine['e'])
        color_map = 'jet'
        size_c = 80
        size_m = 50
        size_f = 20
        e_data_list_c = []
        e_rslt_list_c = []
        e_data_list_m = []
        e_rslt_list_m = []
        e_data_list_f = []
        e_rslt_list_f = []
        th_data_list_c = []
        th_rslt_list_c = []
        th_data_list_m = []
        th_rslt_list_m = []
        th_data_list_f = []
        th_rslt_list_f = []
        if output_e:
            e_data_list_c += [Rice95_v_many_coarse['e']]
            e_rslt_list_c += [e0_bar_list_0[0]]

            e_data_list_m += [Rice95_v_many_medium['e']]
            e_rslt_list_m += [e0_bar_list_0[1]]

            e_data_list_f += [Rice95_v_many_fine['e']]
            e_rslt_list_f += [e0_bar_list_0[2]]
            len0 = 3

            e_data_list_c += Willetts89_v_many_coarse['e']
            e_rslt_list_c += e0_bar_list_0[len0:len0 + len_W89_c]
            len0 += len_W89_c

            e_data_list_m += Willetts89_v_many_medium['e']
            e_rslt_list_m += e0_bar_list_0[len0:len0 + len_W89_m]
            len0 += len_W89_m

            e_data_list_f += Willetts89_v_many_fine['e']
            e_rslt_list_f += e0_bar_list_0[len0:len0 + len_W89_f]
            len0 += len_W89_f

            e_data_list = e_data_list_c + e_data_list_m + e_data_list_f
            e_rslt_list = e_rslt_list_c + e_rslt_list_m + e_rslt_list_f

        if output_theta2:
            th_data_list_c += [Rice95_v_many_coarse['ang_re']]
            th_rslt_list_c += [np.degrees(theta20_bar_list_0[0])]

            th_data_list_m += [Rice95_v_many_medium['ang_re']]
            th_rslt_list_m += [np.degrees(theta20_bar_list_0[1])]

            th_data_list_f += [Rice95_v_many_fine['ang_re']]
            th_rslt_list_f += [np.degrees(theta20_bar_list_0[2])]
            len0 = 3

            th_data_list_c += Willetts89_v_many_coarse['ang_re']
            th_rslt_list_c += np.degrees(theta20_bar_list_0[len0:len0 + len_W89_c]).tolist()
            len0 += len_W89_c

            th_data_list_m += Willetts89_v_many_medium['ang_re']
            th_rslt_list_m += np.degrees(theta20_bar_list_0[len0:len0 + len_W89_m]).tolist()
            len0 += len_W89_m

            th_data_list_f += Willetts89_v_many_fine['ang_re']
            th_rslt_list_f += np.degrees(theta20_bar_list_0[len0:len0 + len_W89_f]).tolist()
            len0 += len_W89_f

            th_data_list = th_data_list_c + th_data_list_m + th_data_list_f
            th_rslt_list = th_rslt_list_c + th_rslt_list_m + th_rslt_list_f

        e_data_array = np.array(e_data_list)
        e_rslt_array = np.array(e_rslt_list)
        e_corr = np.corrcoef(e_data_array, e_rslt_array)[0, 1]
        th_data_array = np.array(th_data_list)
        th_rslt_array = np.array(th_rslt_list)
        th_corr = np.corrcoef(th_data_array, th_rslt_array)[0, 1]
        total_corr = np.sqrt(e_corr* th_corr)
        with open('epsilon_nu_corr.txt', 'a') as f:
            f.write(f'{epsilon:.4f}, {nu:.4f}, {total_corr:.4f}\n')
        data_matrix = np.vstack((e_data_array, th_data_array)).T
        rslt_matrix = np.vstack((e_rslt_array, th_rslt_array)).T
        cca = CCA(n_components=1)
        cca.fit(data_matrix, rslt_matrix)
        data_c, result_c = cca.transform(data_matrix, rslt_matrix)
        M = data_c[:, 0]
        N = result_c[:, 0]
        max_corr = np.corrcoef(M, N)[0, 1]
        data_weights = cca.x_weights_
        rslt_weights = cca.y_weights_
        w_data = data_weights.flatten()
        w_rslt = rslt_weights.flatten()
        w_data_e = w_data[0]
        w_data_th = w_data[1]
        w_rslt_e = w_rslt[0]
        w_rslt_th = w_rslt[1]
        with open('epsilon_nu_CCA_weights.txt', 'a') as f:
            f.write(f'{epsilon:.4f}, {nu:.4f}, {max_corr:.4f}, {w_data_e:.4f}, {w_data_th:.4f}, {w_rslt_e:.4f}, {w_rslt_th:.4f}\n')

plt.figure(1, figsize=(10, 8))
corr_data = np.loadtxt('epsilon_nu_corr.txt', delimiter=',')
read_epsilon = corr_data[:, 0]
read_nu = corr_data[:, 1]
read_corr = corr_data[:, 2]
unique_epsilon = np.unique(read_epsilon)
unique_nu = np.unique(read_nu)
corr = read_corr.reshape(len(unique_nu), len(unique_epsilon))
plt.contourf(unique_epsilon, unique_nu, corr, cmap='viridis', levels=100)
plt.colorbar(ticks=np.linspace(0, 1, 11), label='Correlation Coefficient')
plt.xlabel('Epsilon')
plt.ylabel('Nu')

plt.figure(2, figsize=(10, 8))
CCA_data = np.loadtxt('epsilon_nu_CCA_weights.txt', delimiter=',')
read_CCA_corr = CCA_data[:, 2]
CCA_corr = read_CCA_corr.reshape(len(unique_nu), len(unique_epsilon))
plt.contourf(unique_epsilon, unique_nu, CCA_corr, cmap='jet', levels=100)
plt.colorbar(ticks=np.linspace(0, 1, 11), label='CCA Correlation Coefficient')
plt.xlabel('Epsilon')
plt.ylabel('Nu')

plt.figure(3, figsize=(10, 8))
read_w_data_e = CCA_data[:, 3]
read_w_data_th = CCA_data[:, 4]
w_data_e = read_w_data_e.reshape(len(unique_nu), len(unique_epsilon))
w_data_th = read_w_data_th.reshape(len(unique_nu), len(unique_epsilon))
plt.quiver(unique_epsilon, unique_nu, w_data_e, w_data_th, color='k')
plt.xlabel('Epsilon')
plt.ylabel('Nu')

plt.figure(4, figsize=(10, 8))
read_w_rslt_e = CCA_data[:, 5]
read_w_rslt_th = CCA_data[:, 6]
w_rslt_e = read_w_rslt_e.reshape(len(unique_nu), len(unique_epsilon))
w_rslt_th = read_w_rslt_th.reshape(len(unique_nu), len(unique_epsilon))
plt.quiver(unique_epsilon, unique_nu, w_rslt_e, w_rslt_th, color='r')
plt.xlabel('Epsilon')
plt.ylabel('Nu')

plt.show()

max_indices = np.argpartition(read_corr, -10)[-10:]  # Get indices of the 10 largest corr values
for max_index in max_indices:
    print(f'Max corr: {read_corr[max_index]:.4f} at Epsilon: {read_epsilon[max_index]:.4f}, Nu: {read_nu[max_index]:.4f}')
max_indices = np.argpartition(read_CCA_corr, -10)[-10:]  # Get indices of the 10 largest CCA corr values
for max_index in max_indices:
    print(f'Max CCA corr: {read_CCA_corr[max_index]:.4f} at Epsilon: {read_epsilon[max_index]:.4f}, Nu: {read_nu[max_index]:.4f}, Weights: Exp_e:{read_w_data_e[max_index]:.4f}, Exp_th:{read_w_data_th[max_index]:.4f}, Mod_e:{read_w_rslt_e[max_index]:.4f}, Mod_th:{read_w_rslt_th[max_index]:.4f}')