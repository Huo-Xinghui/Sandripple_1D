import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import truncnorm

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

def E_out_min(d1, d2, d3, rho, g):
    v1 = np.pi * d1**3 / 6
    m1 = rho * v1
    x = d2*(d1+d2+d3) - d1*d3
    y = (d1 + d2)*(d2 + d3)
    theta = np.arccos(x/y)
    E_min_norm = (d1 + d2)*(np.sqrt(2*(1 - np.sin(theta))) + (1 - np.sin(theta)))/4/d1
    E_min = E_min_norm*m1*g*d1
    return E_min_norm, E_min

def E_out(E_in, E_min_mean, e):
    lamb = 2*np.log((1 - e**2)*E_in/E_min_mean)
    sigma = np.sqrt(lamb)*np.log(2)
    mu = np.log((1 - e**2)*E_in) - lamb*np.log(2)
    E_out_mean = np.exp(mu + sigma**2/2)
    E_out_sigma = np.exp(2*mu + sigma**2)*(np.exp(sigma**2) - 1)
    return E_out_mean, E_out_sigma

def generate_truncated_lognormal(mu, sigma, dmin, dmax, size = 1):
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

#-------------------------------------------------------------------
distribution = 3 # 0:uniform, 1:lognormal, 2:bidisperse, 3:polydisperse
plot_type = 0 # 0: Ec, 1: d
d_min = 1e-4
d_max = 7e-4
normal_E = 4e-4
normal_D = 1e-4
mu1 = (d_min + d_max) * 0.3  # 第一个峰靠左
mu2 = (d_min + d_max) * 0.7  # 第二个峰靠右
sigma1 = (d_max - d_min) * 0.1
sigma2 = (d_max - d_min) * 0.1
weight1 = 0.5 # 第一个峰的权重
rho = 2650
g = 9.8*(1 - 1.263/rho)
e = 0.5
num_samples = 10000
d_mid = (d_min + d_max) / 2
m_in = np.pi * d_mid**3 / 6 * rho
E_in = m_in * g * d_mid * 10
max_attempts = num_samples + 1
#------------------------------------------------------------------

Ec_list1 = []
Ec_list2 = []
Ec_list3 = []
Ec_norm_list1 = []
Ec_norm_list2 = []
Ec_norm_list3 = []
d1_list = []
d2_list = []
d3_list = []

attempts = 0

while len(d3_list) < num_samples and attempts < max_attempts:
    attempts += 1
    if distribution == 0:
        d1 = np.random.uniform(d_min, d_max)
        d3 = np.random.uniform(d_min, d_max)
        d2_min = (-(d1+d3)+np.sqrt((d1+d3)**2+4*d1*d3))/2
        d2_min = max(d_min, d2_min)
        d2 = np.random.uniform(d2_min, d_max)
        d1_list.append(d1)
        d2_list.append(d2)
        d3_list.append(d3)
    elif distribution == 1:
        mu, sigma = get_normal_params(normal_E, normal_D)
        d1 = generate_truncated_lognormal(mu, sigma, d_min, d_max, 1)[0]
        d3 = generate_truncated_lognormal(mu, sigma, d_min, d_max, 1)[0]
        d2_min = (-(d1+d3)+np.sqrt((d1+d3)**2+4*d1*d3))/2
        d2_min = max(d_min, d2_min)
        d2 = generate_truncated_lognormal(mu, sigma, d2_min, d_max, 1)[0]
        d1_list.append(d1)
        d2_list.append(d2)
        d3_list.append(d3)
    elif distribution == 2:
        r = np.random.uniform(0.5,1.5,size=3)
        rr = np.floor(r)
        d_bin = d_max - d_min
        d1 = d_min + rr[0]*d_bin
        d2 = d_min + rr[1]*d_bin
        d3 = d_min + rr[2]*d_bin
        d1_list.append(d1)
        d2_list.append(d2)
        d3_list.append(d3)
    elif distribution == 3:
        d1, d3 = np.random.choice(
            generate_bimodal(mu1, sigma1, mu2, sigma2, weight1, d_min, d_max, 1000),
            size=2, replace=False
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
            d1_list.append(d1)
            d2_list.append(d2)
            d3_list.append(d3)

    Ec_norm, Ec = E_out_min(d1, d2, d3, rho, g)
    Ec_list1.append(Ec)
    Ec_norm_list1.append(Ec_norm)

    d3 = d2
    Ec_norm2, Ec = E_out_min(d1, d2, d3, rho, g)
    Ec_list2.append(Ec)
    Ec_norm_list2.append(Ec_norm2)

    d1 = d2
    Ec_norm3, Ec = E_out_min(d1, d2, d3, rho, g)
    Ec_list3.append(Ec)
    Ec_norm_list3.append(Ec_norm3)

Ec_norm_mean = [np.mean(Ec_norm_list1), np.mean(Ec_norm_list2), np.mean(Ec_norm_list3)]
d2_list_cubed = np.array(d2_list)**3
d_mean_cubed = (np.mean(d2_list_cubed))**(1/3)
d_mean = np.mean(d2_list)
d50 = np.percentile(d2_list, 50)
d90 = np.percentile(d2_list, 90)
m1 = np.pi*d_mean**3/6*rho
m1_cubed = np.pi*d_mean_cubed**3/6*rho
m1_50 = np.pi*d50**3/6*rho
m1_90 = np.pi*d90**3/6*rho
Ec0 = Ec_norm_mean[2]*m1*g*d_mean
Ec0_cubed = Ec_norm_mean[2]*m1_cubed*g*d_mean_cubed
Ec0_50 = Ec_norm_mean[2]*m1_50*g*d50
Ec0_90 = Ec_norm_mean[2]*m1_90*g*d90
Ec_mean = [np.mean(Ec_list2), np.mean(Ec_list3), Ec0_cubed, Ec0_50, Ec0_90, Ec0, np.mean(Ec_list1)]
Ec_mean_uni = Ec_mean/Ec_mean[6]
Eout_E = []
Eout_D = []
for i in range(7):
    Eout_mean, Eout_sigma = E_out(E_in, Ec_mean[i], e)
    Eout_E.append(Eout_mean)
    Eout_D.append(Eout_sigma)
Eout_E_uni = Eout_E/Eout_E[6]
Eout_D_uni = Eout_D/Eout_D[6]

if plot_type == 0:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(Ec_norm_list1, bins=50, density=True, alpha=1.0, label='d1, d2, d3', color='r')
    ax.hist(Ec_norm_list2, bins=50, density=True, alpha=0.5, label='d1, d2=d3', color='g')
    ax.axvline(Ec_norm_mean[0], color='r', linestyle='dashed', linewidth=2, label='d1, d2, d3')
    ax.axvline(Ec_norm_mean[1], color='g', linestyle='dashed', linewidth=2, label='d1, d2 = d3')
    ax.axvline(Ec_norm_mean[2], color='b', linestyle='dashed', linewidth=2, label='d1 = d2 = d3')
    ax.set_xlabel('Ec*')
    #ax.set_ylim(0,3000)
    ax.legend()

    # 插入插图
    inax = inset_axes(ax, width="50%", height="50%", loc='right')
    colors = ['k', 'g', 'b', 'c', 'm', 'y']
    labels = ['d1, d2=d3', 'd1=d2=d3', 'cubed', 'd50', 'd90', 'averaged']
    for i in range(len(Ec_mean_uni[:-1])):
        inax.plot(Ec_mean_uni[i], Eout_E_uni[i], marker='o', color=colors[i], linestyle='', label=labels[i])
        inax.plot(Ec_mean_uni[i], Eout_D_uni[i], marker='*', color=colors[i], linestyle='', label=labels[i])
    inax.plot(Ec_mean_uni[-1], Eout_E_uni[-1], marker='o', color='r', fillstyle='none', linestyle='', label='d1, d2, d3')
    inax.plot(Ec_mean_uni[-1], Eout_D_uni[-1], marker='*', color='r', fillstyle='none', linestyle='', label='d1, d2, d3')
    inax.set_xlabel('Unified Ec_mean')
    inax.set_ylabel('Unified Eout_mean/std')
    inax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
elif plot_type == 1:
    plt.hist(d1_list, bins=50, density=True, alpha=0.5, label='d1', color='r')
    plt.hist(d2_list, bins=50, density=True, alpha=0.5, label='d2', color='g')
    plt.hist(d3_list, bins=50, density=True, alpha=0.5, label='d3', color='b')
    plt.xlabel('d')
    plt.xlim(d_min, d_max)
    plt.legend()

plt.show()