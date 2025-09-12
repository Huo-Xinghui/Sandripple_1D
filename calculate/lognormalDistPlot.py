import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import truncnorm

def generate_truncated_lognormal(mu, sigma, dmin, dmax, size):
    """Generate truncated log-normal samples."""
    a = (np.log(dmin) - mu) / sigma
    b = (np.log(dmax) - mu) / sigma
    samples = truncnorm.rvs(a, b, loc=mu, scale=sigma, size=size)
    return np.exp(samples)


# 根据对数正态分布的均值和标准差求正态分布的参数
def get_normal_params(log_mean, log_std):
    # 解关于sigma的方程
    def equation(sigma):
        mu = math.log(log_mean) - 0.5 * sigma ** 2
        return log_std ** 2 - math.exp(2 * mu + sigma ** 2) * (math.exp(sigma ** 2) - 1)
    # 简单的二分法求根
    left, right = 0.001, 10
    while right - left > 1e-6:
        mid = (left + right) / 2
        if equation(mid) > 0:
            left = mid
        else:
            right = mid
    sigma = mid
    mu = math.log(log_mean) - 0.5 * sigma ** 2
    return mu, sigma


#log_mean = 2.4e-4
#log_std = 0.5e-4
xmin = 1e-4
xmax = 10e-4
xnum = 1000
dx = xmax / xnum
ddx = dx * 0.5
#mu, sigma = get_normal_params(log_mean, log_std)
mu = -7.5
sigma = 0.3246
log_std = np.sqrt((math.exp(sigma**2)-1)*math.exp(2*mu+sigma**2))
# 生成数据点
x_array = generate_truncated_lognormal(mu, sigma, xmin, xmax, 1000000)
x_mean = np.mean(x_array)
x = np.linspace(xmin+dx, xmax, xnum)
x = x - ddx
# 对数正态分布概率密度函数
y = (1 / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2))
# 归一化
y = y / y.sum()
# 计算x50
cum_y = np.cumsum(y)
idx = np.argmax(cum_y > 0.5)
x50 = x[idx] # 中位数
# 计算对数正态分布原始x50
x50_orig = np.exp(mu) * np.exp(0.5 * sigma**2)
# 计算x90
idx = np.argmax(cum_y > 0.9)
x90 = x[idx] # 90%分位数
idx = np.argmax(cum_y > 0.1)
x10 = x[idx] # 10%分位数
# 绘制曲线
plt.figure(1, figsize=(10, 6))
plt.plot(x, y, 'r-', label='lognormal')
plt.axvline(x50, color='r', linestyle='--', label=f'x50 = {x50:.4e}')
plt.axvline(x90, color='b', linestyle='--', label=f'x90 = {x90:.4e}')
plt.axvline(x10, color='b', linestyle=':', label=f'x10 = {x10:.4e}')
plt.axvline(x_mean, color='g', linestyle='--', label=f'xmean = {x_mean:.4e}')
plt.text(xmax, 0, f'$\\eta$ = {x90/x50:.4e}', ha='right', va='bottom')
plt.title(f'mu={mu:.4f}, sigma={sigma:.4f}, std={log_std:.4e}')
plt.xlabel('Value')
plt.ylabel('Probability Density')
plt.grid(True)
plt.legend()

plt.figure(2, figsize=(10, 6))
# 绘制累计分布曲线
#plt.plot(x, cum_y, 'b-', label='lognormal1')
plt.semilogx(x, cum_y, 'b-', label='cumulative')
plt.title(f'mu={mu:.4f}, sigma={sigma:.4f}, pdf({xmax:.4f})={y[-1]:.4e}')
plt.xlabel('Value')
plt.ylabel('Cumulative Probability')
plt.grid(True)
plt.legend()
plt.show()