import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
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

d1 = 0.006
v1 = 26.0
gamma_ej = 0.049 # the fraction of remaining energy to eject particles
rho = 1768.39
g = 9.8*(1 - 1.263/rho)
epsilon = 0.73
nu = -0.18
num_samples = 100
#------------------------------------------------------------------
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

theta1_array = Beladjine07_v26['ang_in']
iteration_num = len(theta1_array)

e_bar_list = []
ez_bar_list = []
ex_bar_list = []
theta2_bar_list = []
Nej_bar_list = []
vn_bar_list = []
for n in tqdm(range(iteration_num)):
    theta1 = np.radians(theta1_array[n])
    e_list = []
    ex_list = []
    ez_list = []
    theta2_list = []
    Nej_list = []
    vn_list = []
    current_n = 0
    while current_n < num_samples:
        current_n += 1
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
        x_hat = np.random.uniform(x_min, x_max)
        e, evx, evz = calculate_e(alpha, beta, x_hat, theta1)
        m1 = rho*(d1**3*np.pi/6)
        m2 = rho*(d2**3*np.pi/6)
        E1 = 0.5*m1*v1**2
        Ec = m2*g*d2
        k_max = (1.0 - e**2)*E1/Ec
        lambda_ej = 2.0*np.log(k_max)
        mu_ej = np.log((1.0 - e**2)*E1) - lambda_ej*np.log(2.0)
        sigma_ej = np.sqrt(lambda_ej)*np.log(2.0)
        En_bar = np.exp(mu_ej + sigma_ej**2/2)
        En = 0.0
        while En < Ec:
            En = np.random.lognormal(mu_ej, sigma_ej, size=1)[0]
        vn = np.sqrt(2.0*En/m2)
        Nej = gamma_ej*(1.0 - e**2)*E1/En_bar/2.0*erfc((np.log(Ec) - mu_ej)/(np.sqrt(2.0)*sigma_ej))
        Nej = max(Nej, 0.0)
        vn_list.append(vn)
        Nej_list.append(Nej)

        ez = evz #/np.sin(theta1)
        ex = evx #/np.cos(theta1)
        theta2 = np.arctan(evz/evx)
        if theta2 < 0:
            theta2 += np.pi
        e_list.append(e)
        ez_list.append(ez)
        ex_list.append(ex)
        theta2_list.append(theta2)
        v2 = e*v1
        v2_z = v1*evz
        v2_x = v1*evx
    ez_bar = np.mean(ez_list)
    ex_bar = np.mean(ex_list)
    e_bar = np.sqrt(ex_bar**2 + ez_bar**2)
    theta2_bar = np.mean(theta2_list)
    e_bar_list.append(e_bar)
    ez_bar_list.append(ez_bar)
    ex_bar_list.append(ex_bar)
    theta2_bar_list.append(theta2_bar)
    vn_bar = np.mean(vn_list)
    vn_bar_list.append(vn_bar)
    Nej_bar = np.mean(Nej_list)
    Nej_bar_list.append(Nej_bar)

color_map = 'jet'
if output_e:
    plt.figure(1, figsize=(8, 6))
    colors = Beladjine07_v26['ang_in']
    x_array = Beladjine07_v26['e']
    y_array = e_bar_list
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, marker='o', label='Beladjine07')
    min_xy = min(x_array + y_array)
    max_xy = max(x_array + y_array)
    line_points = [min_xy, max_xy]
    plt.plot(line_points, line_points, 'k--', label='y=x')
    plt.xlabel('$e_{exp}$')
    plt.ylabel('$e_{sim}$')
    plt.colorbar(label='Impact angle (degrees)')
    plt.legend()

if output_theta2:
    plt.figure(2, figsize=(8, 6))
    colors = Beladjine07_v26['ang_in']
    x_array = Beladjine07_v26['ang_re']
    y_array = np.degrees(theta2_bar_list)
    y_array = y_array.tolist()
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, marker='o', label='Beladjine07')
    min_xy = min(x_array + y_array)
    max_xy = max(x_array + y_array)
    line_points = [min_xy, max_xy]
    plt.plot(line_points, line_points, 'k--', label='y=x')
    plt.xlabel('$\\theta_{re,exp}$ (degrees)')
    plt.ylabel('$\\theta_{re,sim}$ (degrees)')
    plt.colorbar(label='Impact angle (degrees)')
    plt.legend()

if output_Nej:
    plt.figure(3, figsize=(8, 6))
    colors = Beladjine07_v26['ang_in']
    x_array = Beladjine07_v26['Nej']
    y_array = Nej_bar_list
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, marker='o', label='Beladjine07')
    min_xy = min(x_array + y_array)
    max_xy = max(x_array + y_array)
    line_points = [min_xy, max_xy]
    plt.plot(line_points, line_points, 'k--', label='y=x')
    plt.xlabel('$N_{ej,exp}$')
    plt.ylabel('$N_{ej,sim}$')
    plt.colorbar(label='Impact angle (degrees)')
    plt.legend()

if output_vn:
    plt.figure(4, figsize=(8, 6))
    colors = Beladjine07_v26['ang_in']
    x_array = Beladjine07_v26['v_ej']
    y_array = vn_bar_list
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, marker='o', label='Beladjine07')
    min_xy = min(x_array + y_array)
    max_xy = max(x_array + y_array)
    line_points = [min_xy, max_xy]
    plt.plot(line_points, line_points, 'k--', label='y=x')
    plt.xlabel('$v_{n,exp}$ (m/s)')
    plt.ylabel('$v_{n,sim}$ (m/s)')
    plt.colorbar()

plt.show()