import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import truncnorm

def generate_bimodal(mu1, sigma1, mu2, sigma2, weight1, dmin, dmax, size):
    """生成截断的双峰分布样本"""
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

    combined = np.concatenate([samples1, samples2])
    np.random.shuffle(combined)
    return combined

def generate_constrained_triplets(
    mu1, sigma1, mu2, sigma2, weight1,
    dmin, dmax, num_samples, max_attempts=10000
):
    valid_samples = []
    attempts = 0

    while len(valid_samples) < num_samples and attempts < max_attempts:
        attempts += 1

        # 生成 d1 和 d3（同一双峰分布）
        d1, d3 = np.random.choice(
            generate_bimodal(mu1, sigma1, mu2, sigma2, weight1, dmin, dmax, 1000),
            size=2, replace=False
        )

        # 计算 d2 的最小允许值
        sqrt_term = np.sqrt((d1 + d3)**2 + 4 * d1 * d3)
        d2_min = (-(d1 + d3) + sqrt_term) / 2
        lower = max(dmin, d2_min)
        upper = dmax

        if lower > upper:
            continue  # 跳过无效组合

        # 生成 d2（同一双峰分布，截断到 [lower, upper]）
        d2_candidates = generate_bimodal(mu1, sigma1, mu2, sigma2, weight1, lower, upper, 100)
        if len(d2_candidates) == 0:
            continue

        d2 = np.random.choice(d2_candidates)

        # 验证所有变量在范围内
        if (dmin <= d1 <= dmax and
            dmin <= d2 <= dmax and
            dmin <= d3 <= dmax and
            d2 >= d2_min
        ):
            valid_samples.append((d1, d2, d3))

    return valid_samples[:num_samples]

# 取值范围
dmin = 1.0
dmax = 10.0

# 双峰分布参数
mu1 = (dmin + dmax) * 0.3  # 第一个峰靠左
mu2 = (dmin + dmax) * 0.7  # 第二个峰靠右
sigma1 = (dmax - dmin) * 0.1
sigma2 = (dmax - dmin) * 0.1
weight1 = 0.5 # 第一个峰的权重

# 生成大样本用于可视化
d1_samples = generate_bimodal(mu1, sigma1, mu2, sigma2, weight1, dmin, dmax, 10000)
d2_samples = [triplet[1] for triplet in generate_constrained_triplets(
    mu1, sigma1, mu2, sigma2, weight1, dmin, dmax, 10000
)]

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.hist(d1_samples, bins=50, density=True, alpha=0.6, label="d1/d3 Distribution")
plt.title("Original Bimodal Distribution")

plt.subplot(1, 2, 2)
plt.hist(d2_samples, bins=50, density=True, alpha=0.6, color="orange", label="Constrained d2")
plt.title("d2 Distribution after Constraints")
plt.show()