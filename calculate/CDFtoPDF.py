import numpy as np
from scipy.stats import norm

# 示例数据：已知的 x_i 和对应的累积概率 p_i
x_data = np.array([0.000150, 0.000180, 0.000212, 0.000250, 0.000300, 0.000355, 0.000425, 0.000500])  # 已知的 x 值
p_data = np.array([0.00785, 0.08593, 0.28992, 0.52017, 0.76283, 0.91208, 0.97593, 1.0])  # 对应的 P(X <= x_i)

# 步骤1：将累积概率转换为标准正态分位数 z_i
z_data = norm.ppf(p_data)

# 步骤2：对 ln(x_i) 和 z_i 进行线性回归
ln_x_data = np.log(x_data)

# 使用最小二乘法拟合：ln_x = mu + sigma * z
A = np.vstack([np.ones_like(z_data), z_data]).T
mu_opt, sigma_opt = np.linalg.lstsq(A, ln_x_data, rcond=None)[0]

print(f"拟合参数: μ = {mu_opt:.4f}, σ = {sigma_opt:.4f}")