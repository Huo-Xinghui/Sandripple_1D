import numpy as np
import matplotlib.pyplot as plt

r1 = 1.5
r2 = 0.5
r3 = 1.5
r4 = 0.99
alpha_list = np.linspace(-np.pi/2, np.pi/2, 100)
alpha_array = np.array(alpha_list)
d13 = r1 + r3
d14 = r1 + r4
d23 = r2 + r3
d24 = r2 + r4
d34 = r3 + r4
cos_gamma = (d23**2 + d24**2 - d34**2)/(2*d23*d24)
gamma = np.arccos(cos_gamma)
K = d14**2 - d24**2 - d13**2 + d23*d24*cos_gamma
R = d23*d34
theta = np.arctan(-d24*np.sin(gamma)/(d23 - d24*cos_gamma))
target = 1/2*(theta + np.arccos(K/R))
sin_target = np.sin(target)
K = r1*(r2-r4) + r2*(r2+r4)
R = r2 + r4
if K < 0:
    target1 = 0
else:
    cos_beta = np.sqrt(K)/R
    beta = np.arccos(cos_beta)
    target1 = gamma - beta
tan_target2 = (np.sin(gamma)*np.sqrt(K) - np.cos(gamma)*np.sqrt(R**2 - K))/(np.cos(gamma)*np.sqrt(K) + np.sin(gamma)*np.sqrt(R**2 - K))
target2 = np.arctan(tan_target2)
# 坐标
O1O4_list = []
O2 = np.array([0, 0, 0])
for alpha in alpha_array:
    a = d23*np.cos(alpha)
    b = d23*np.sin(alpha)
    O3 = np.array([-a, -b, 0])
    beta = gamma - alpha
    c = d24*np.cos(beta)
    d = d24*np.sin(beta)
    O4 = np.array([-c, d, 0])
    e = d13*np.sqrt(1 - (b/d13)**2)
    O1 = np.array([-a, 0, e])
    O1O4 = np.linalg.norm(O1 - O4)
    O1O4_list.append(O1O4)
O1O4_array = np.array(O1O4_list)
alpha_array = alpha_array * 180 / np.pi  # 转换为角度制
target = target * 180 / np.pi  # 转换为角度制
target1 = target1 * 180 / np.pi  # 转换为角度制
target2 = target2 * 180 / np.pi  # 转换为角度制
plt.figure(figsize=(8, 8))
plt.plot(alpha_array, O1O4_array, label='O1O4 Distance')
plt.axvline(x=target, color='r', linestyle='-', label='Target Distance 1')
plt.axvline(x=target1, color='b', linestyle='-', label='Target Distance 2')
plt.axvline(x=target2, color='m', linestyle='--', label='Target Distance 3')
plt.axhline(y=d14, color='g', linestyle='--', label='d14 Distance')
plt.xlabel('Alpha (rad)')
plt.ylabel('Distance')
plt.legend()
plt.grid()
plt.show()
