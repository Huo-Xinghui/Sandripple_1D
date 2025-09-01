import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize

def objective(params, L, R, x_vals, G_vals):
    mu, sigma = params
    if sigma <= 0:  # 确保sigma为正
        return np.inf

    # 计算归一化常数
    F_L = norm.cdf((np.log(L) - mu) / sigma)
    F_R = norm.cdf((np.log(R) - mu) / sigma)
    norm_const = F_R - F_L

    if norm_const <= 0:  # 无效参数
        return np.inf

    # 计算理论CDF
    F_x = norm.cdf((np.log(x_vals) - mu) / sigma)
    G_theory = (F_x - F_L) / norm_const

    # 返回残差平方和
    return np.sum((G_vals - G_theory)**2)

# 示例数据
L = 0.5e-4  # 截断下限
R = 20e-4   # 截断上限
x_vals = np.array([1.36, 3.08, 5.29, 6.54, 7.94, 9.70, 13.0])*1.0e-4  # 分位数点
G_vals = np.array([0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9])           # 对应累积概率

# 初始猜测 (基于非截断分布)
initial_guess = [np.mean(np.log(x_vals)), np.std(np.log(x_vals))]

# 执行优化
result = minimize(objective, initial_guess, args=(L, R, x_vals, G_vals),
                  bounds=[(None, None), (1e-8, None)])  # 确保sigma>0

mu_hat, sigma_hat = result.x
print(f"拟合参数: mu = {mu_hat:.4f}, sigma = {sigma_hat:.4f}")