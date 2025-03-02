import numpy as np
import matplotlib.pyplot as plt

# 定义参数
d1 = 0.0003  # 入射颗粒直径
d2 = 0.001  # 床面颗粒直径
a = 0.1  # 入射角度
eN = 0.9  # 垂向恢复系数
eT = 0.0 # 切向恢复系数

# 计算中间变量
dd1 = d1 / (0.5 * (d1 + d2))
dd2 = d2 / (0.5 * (d1 + d2))
eta = eN * dd1**3 / (dd1**3 + eN * dd2**3)
alpha = (1.0 + eN) / (1.0 + eta) - 1.0
beta = 1.0 - (2.0 / 7.0) * (1.0 - eT) / (1.0 + eta)
gama = 4.0 / 9.0 * beta**2 / (alpha + beta)**2 / dd2

# 定义 PDF 函数
xmin = -a # x 的最小值
xmax = np.sqrt(a / gama) * 2 - a # x 的最大值
def pdf(x, a, gama):
    if xmin <= x <= xmax:
        return gama * (x + a) / a * np.log(2.0 * a / gama / (x + a)**2)
    else:
        return 0.0

# 生成 x 值
x = np.linspace(xmin+1.0e-8, xmax, 1000)

# 计算 PDF 值
y = np.array([pdf(xi, a, gama) for xi in x])

# 计算积分
integral = np.trapz(y, x)

# 绘制概率分布曲线
plt.figure(figsize=(10, 6))
plt.plot(x, y, label=f'$d1={d1}$, $d2={d2}$, $angIn={a}$')

# 设置标题和标签
plt.title(f'PDF of angOut, integral={integral:.4f}')
plt.xlabel('angOut')
plt.ylabel('p(x|a)')

# 显示图例
plt.legend()

# 显示图形
plt.show()