import numpy as np
import matplotlib.pyplot as plt
import math


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


log_mean = 3.0e-4
log_std = 1.0e-4
xmin = 1e-4
xmax = 7e-4
xnum = 20
dx = xmax / xnum
ddx = dx * 0.5
mu, sigma = get_normal_params(log_mean, log_std)
mu1 = np.log(log_mean) - 0.5 * log_std ** 2
sigma1 = np.log(1 + log_std**2 / log_mean**2) ** 0.5
# 生成数据点
x = np.linspace(xmin+dx, xmax, xnum)
x = x - ddx
# 对数正态分布概率密度函数
y = (1 / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2))
y1 = (1 / (x * sigma1 * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu1) ** 2 / (2 * sigma1 ** 2))
# 归一化
y1 = y
y = y / y.sum()
y1 = y1 * dx
# 绘制曲线
plt.plot(x, y, 'ro-', label='lognormal')
plt.plot(x, y1, 'bo-', label='lognormal1')
plt.title(f'mu={mu:.4f}, sigma={sigma:.4f}, pdf({xmax:.4f})={y[-1]:.4e}')
plt.xlabel('Value')
plt.ylabel('Probability Density')
plt.grid(True)
plt.legend()
plt.show()