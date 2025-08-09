import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import truncnorm
from tqdm import tqdm
from scipy.special import erfc
from scipy.integrate import quad
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.usetex'] = True
label_size = 20
ticks_size = 15
marker_size = 10
#-------------------------------------------------------------------
output_rebound_ratio = False # output rebound ratio
output_e = True # output e
output_ex = False # output e_x
output_ez = False # output e_z
output_theta2 = True # output rebound angle
output_theta2_dist = False # output rebound angle distribution
output_vn = False # output eject velocity
output_Nej = False # output eject number


variation_param = 1 # 0: v1, 1: theta1
shallow = False # shallow impact
simplify = False # first order approximation
impactor = 3 #impactor type: 0: bed, 1: coarse, 2: medium, 3: fine
#------------------------------------------------------------------
"""other's data"""
Beladjine07_v26= {
    'ang_in': [10, 20, 40, 60, 90],
    'ang_re': [21.21, 26.68, 33.87, 40.94, 90.01],
    'e': [0.77, 0.61, 0.43, 0.26, 0.22]
}
Zhou06_v_many= {
    'ang_in': [8, 11.5],
    'ang_re_mean': [46.97, 47.53],
    'ang_re_std': [0.59, 0.75],
    'e_mean': [0.63, 0.61],
    'e_std': [0.006, 0.009]
}
Rice95_v_many_coarse= {
    'ang_in': np.mean([13.94, 14.75, 14.73, 15.04]),
    'ang_in_std': np.std([13.94, 14.75, 14.73, 15.04]),
    'ang_re': np.mean([20.95, 23.03, 22.55, 25.63]),
    'ang_re_std': np.std([20.95, 23.03, 22.55, 25.63]),
    'e': np.mean([0.58, 0.56, 0.56, 0.58]),
    'e_std': np.std([0.58, 0.56, 0.56, 0.58]),
    'Nej': np.mean([5.7, 5.09, 5.12, 4.64]),
    'Nej_std': np.std([5.7, 5.09, 5.12, 4.64])
}
Rice95_v_many_medium= {
    'ang_in': np.mean([11.82, 11.47, 11.53, 11.36]),
    'ang_in_std': np.std([11.82, 11.47, 11.53, 11.36]),
    'ang_re': np.mean([29.95, 30.31, 31.06, 29.56]),
    'ang_re_std': np.std([29.95, 30.31, 31.06, 29.56]),
    'e': np.mean([0.57, 0.54, 0.59, 0.58]),
    'e_std': np.std([0.57, 0.54, 0.59, 0.58]),
    'Nej': np.mean([2.68, 3.43, 2.67, 2.76]),
    'Nej_std': np.std([2.68, 3.43, 2.67, 2.76])
}
Rice95_v_many_fine= {
    'ang_in': np.mean([10.85, 10.24, 10.46]),
    'ang_in_std': np.std([10.85, 10.24, 10.46]),
    'ang_re': np.mean([44.63, 38.03, 37.85]),
    'ang_re_std': np.std([44.63, 38.03, 37.85]),
    'e': np.mean([0.52, 0.56, 0.58]),
    'e_std': np.std([0.52, 0.56, 0.58]),
    'Nej': np.mean([1.86, 1.71, 1.58]),
    'Nej_std': np.std([1.86, 1.71, 1.58])
}
Chen18_v_many= {
    'ang_in': [23.2, 21.8, 21.9, 30.7, 30.6, 37.6, 47.1, 46.6, 46],
    'ang_re': [32.83, 35.59, 29.90, 43, 33.28, 44.16, 57.62, 50.97, 39.88],
    'e': [0.38, 0.4, 0.45, 0.32, 0.37, 0.27, 0.2, 0.19, 0.22]
}
Willetts89_v_many_coarse= {
    'ang_in': [12.7, 17.8, 23.2, 27.7],
    'ang_re': [19.1, 25.2, 21.4, 27.2],
    'e': [0.63, 0.57, 0.54, 0.46]
}
Willetts89_v_many_medium= {
    'ang_in': [11.7, 18.2, 21.4, 26.3],
    'ang_re': [24.9, 33.4, 33.3, 44.7],
    'e': [0.61, 0.53, 0.50, 0.40]
}
Willetts89_v_many_fine= {
    'ang_in': [9.5, 15.4, 19.7, 24.9],
    'ang_re': [38.8, 42, 42.2, 42.5],
    'e': [0.57, 0.50, 0.48, 0.46]
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

plt.figure(1, figsize=(8, 6))
"""read e data"""
#data_coarse = np.loadtxt('rebound_e_coarse_7390.txt')
#data_medium = np.loadtxt('rebound_e_medium_7390.txt')
#data_fine = np.loadtxt('rebound_e_fine_7390.txt')
#data_coarse = np.loadtxt('rebound_e_coarse_7692.txt')
#data_medium = np.loadtxt('rebound_e_medium_7692.txt')
#data_fine = np.loadtxt('rebound_e_fine_7692.txt')
data_coarse = np.loadtxt('rebound_e_coarse_7188.txt')
data_medium = np.loadtxt('rebound_e_medium_7188.txt')
data_fine = np.loadtxt('rebound_e_fine_7188.txt')
data_coarse_3D = np.loadtxt('rebound_e_coarse_6173_3D_randk.txt')
data_medium_3D = np.loadtxt('rebound_e_medium_6173_3D_randk.txt')
data_fine_3D = np.loadtxt('rebound_e_fine_6173_3D_randk.txt')
xc = data_coarse[:, 0]
xm = data_medium[:, 0]
xf = data_fine[:, 0]
e_coarse = data_coarse[:, 2]
e_medium = data_medium[:, 2]
e_fine = data_fine[:, 2]
e_coarse_1 = data_coarse[:, 3]
e_medium_1 = data_medium[:, 3]
e_fine_1 = data_fine[:, 3]
xc_3D = data_coarse_3D[:, 0]
xm_3D = data_medium_3D[:, 0]
xf_3D = data_fine_3D[:, 0]
e_coarse_3D = data_coarse_3D[:, 2]
e_medium_3D = data_medium_3D[:, 2]
e_fine_3D = data_fine_3D[:, 2]
e_coarse_3D_1 = data_coarse_3D[:, 3]
e_medium_3D_1 = data_medium_3D[:, 3]
e_fine_3D_1 = data_fine_3D[:, 3]
"""plot read e data"""
if variation_param == 0:
    plt.xlabel('$\\hat{v}_1$', fontsize=label_size)
    x_exp = 'v_in'
elif variation_param == 1:
    plt.xlabel('$\\theta$ (degree)', fontsize=label_size)
    x_exp = 'ang_in'
plt.ylabel('$\\overline{e}$', fontsize=label_size)
strlabel1 = 'Fine'
strlabel2 = 'Medium'
strlabel3 = 'Coarse'
pf = plt.plot(xf, e_fine, 'r--', label=strlabel1)
pm = plt.plot(xm, e_medium, 'g--', label=strlabel2)
pc = plt.plot(xc, e_coarse, 'b--', label=strlabel3)
pf1 = plt.plot(xf, e_fine_1, 'r:', label=strlabel1)
pm1 = plt.plot(xm, e_medium_1, 'g:', label=strlabel2)
pc1 = plt.plot(xc, e_coarse_1, 'b:', label=strlabel3)
pf_3D = plt.plot(xf_3D, e_fine_3D, 'r-', label=strlabel1)
pm_3D = plt.plot(xm_3D, e_medium_3D, 'g-', label=strlabel2)
pc_3D = plt.plot(xc_3D, e_coarse_3D, 'b-', label=strlabel3)
#pf1_3D = ax1.plot(xf_3D, e_fine_3D_1, 'r-.', label=strlabel1)
#pm1_3D = ax1.plot(xm_3D, e_medium_3D_1, 'g-.', label=strlabel2)
#pc1_3D = ax1.plot(xc_3D, e_coarse_3D_1, 'b-.', label=strlabel3)

plt.xlim(0, xc[-1])
plt.ylim(0.4, 0.65)
"""plot other's data"""
#plt.plot(Chen18_v_many[x_exp], Chen18_v_many['e'], 'kP', label='Chen18')
#plt.plot(Beladjine07_v26[x_exp], Beladjine07_v26['e'], 'k*', label='Beladjine07')
#plt.plot(Zhou06_v_many[x_exp], Zhou06_v_many['e_mean'], 'k^', label='Zhou06')
#plt.plot(Rioual20_v_many[x_exp], Rioual20_v_many['e'], 'k.', label='Rioual20')
#plt.plot(Gordon21_v_many[x_exp], Gordon21_v_many['e'], 'k1', label='Gordon21')
#plt.plot(Gordon09_v_many[x_exp], Gordon09_v_many['e'], 'k2', label='Gordon09')
ebf = plt.errorbar(
    Rice95_v_many_fine[x_exp],
    Rice95_v_many_fine['e'],
    yerr=Rice95_v_many_fine['e_std'],
    fmt='ro',
    markerfacecolor='none',
    markeredgecolor='r',
    markersize=marker_size,
    label='Ri95 fine',
    capsize=5
)
ebm = plt.errorbar(
    Rice95_v_many_medium[x_exp],
    Rice95_v_many_medium['e'],
    yerr=Rice95_v_many_medium['e_std'],
    fmt='g^',
    markerfacecolor='none',
    markeredgecolor='g',
    markersize=marker_size,
    label='Ri95 medium',
    capsize=5
)
ebc = plt.errorbar(
    Rice95_v_many_coarse[x_exp],
    Rice95_v_many_coarse['e'],
    yerr=Rice95_v_many_coarse['e_std'],
    fmt='bs',
    markerfacecolor='none',
    markeredgecolor='b',
    markersize=marker_size,
    label='Ri95 coarse',
    capsize=5
)
pwf = plt.plot(Willetts89_v_many_fine[x_exp], Willetts89_v_many_fine['e'], 'ro', markersize=marker_size, label='Willetts89 fine')
pwm = plt.plot(Willetts89_v_many_medium[x_exp], Willetts89_v_many_medium['e'], 'g^', markersize=marker_size, label='Willetts89 medium')
pwc = plt.plot(Willetts89_v_many_coarse[x_exp], Willetts89_v_many_coarse['e'], 'bs', markersize=marker_size, label='Willetts89 coarse')

plt.tick_params(labelsize=ticks_size)

plt.legend([pwf[0], pwm[0], pwc[0]], ['Fine', 'Medium', 'Coarse'], fontsize=ticks_size, loc='best')

plt.figure(2, figsize=(8, 6))
"""read th data"""
#data_coarse = np.loadtxt('rebound_th_coarse_7390.txt')
#data_medium = np.loadtxt('rebound_th_medium_7390.txt')
#data_fine = np.loadtxt('rebound_th_fine_7390.txt')
#data_coarse = np.loadtxt('rebound_th_coarse_7692.txt')
#data_medium = np.loadtxt('rebound_th_medium_7692.txt')
#data_fine = np.loadtxt('rebound_th_fine_7692.txt')
data_coarse = np.loadtxt('rebound_th_coarse_7188.txt')
data_medium = np.loadtxt('rebound_th_medium_7188.txt')
data_fine = np.loadtxt('rebound_th_fine_7188.txt')
data_coarse_3D = np.loadtxt('rebound_th_coarse_6173_3D_randk.txt')
data_medium_3D = np.loadtxt('rebound_th_medium_6173_3D_randk.txt')
data_fine_3D = np.loadtxt('rebound_th_fine_6173_3D_randk.txt')
xc = data_coarse[:, 0]
xm = data_medium[:, 0]
xf = data_fine[:, 0]
th_coarse = data_coarse[:, 2]
th_medium = data_medium[:, 2]
th_fine = data_fine[:, 2]
th_coarse_1 = data_coarse[:, 3]
th_medium_1 = data_medium[:, 3]
th_fine_1 = data_fine[:, 3]
xc_3D = data_coarse_3D[:, 0]
xm_3D = data_medium_3D[:, 0]
xf_3D = data_fine_3D[:, 0]
th_coarse_3D = data_coarse_3D[:, 2]
th_medium_3D = data_medium_3D[:, 2]
th_fine_3D = data_fine_3D[:, 2]
"""plot read th data"""
if variation_param == 0:
    plt.xlabel('$\\hat{v}_1$', fontsize=label_size)
elif variation_param == 1:
    plt.xlabel('$\\theta$ (degree)', fontsize=label_size)
plt.ylabel('$\\overline{\\theta \'}$(degree)', fontsize=label_size)
strlabel1 = 'Fine'
strlabel2 = 'Medium'
strlabel3 = 'Coarse'
pf = plt.plot(xf, th_fine, 'r--', label=strlabel1)
pm = plt.plot(xm, th_medium, 'g--', label=strlabel2)
pc = plt.plot(xc, th_coarse, 'b--', label=strlabel3)
pf1 = plt.plot(xf, th_fine_1, 'r:', label=strlabel1)
pm1 = plt.plot(xm, th_medium_1, 'g:', label=strlabel2)
pc1 = plt.plot(xc, th_coarse_1, 'b:', label=strlabel3)
pf_3D = plt.plot(xf_3D, th_fine_3D, 'r-', label=strlabel1)
pm_3D = plt.plot(xm_3D, th_medium_3D, 'g-', label=strlabel2)
pc_3D = plt.plot(xc_3D, th_coarse_3D, 'b-', label=strlabel3)

plt.xlim(0, xc[-1])
plt.ylim(0, 70)

"""plot other's data"""
#plt.plot(Chen18_v_many[x_exp], Chen18_v_many['ang_re'], 'kP', label='Chen et al. 2018')
#plt.plot(Beladjine07_v26[x_exp], Beladjine07_v26['ang_re'], 'k*', label='Beladjine et al. 2007')
#plt.plot(Zhou06_v_many[x_exp], Zhou06_v_many['ang_re_mean'], 'k^', label='Zhou et al. 2006')
plt.errorbar(
    Rice95_v_many_fine[x_exp],
    Rice95_v_many_fine['ang_re'],
    yerr=Rice95_v_many_fine['ang_re_std'],
    fmt='ro',
    markerfacecolor='none',
    markeredgecolor='r',
    markersize=marker_size,
    capsize=5
)
plt.errorbar(
    Rice95_v_many_medium[x_exp],
    Rice95_v_many_medium['ang_re'],
    yerr=Rice95_v_many_medium['ang_re_std'],
    fmt='g^',
    markerfacecolor='none',
    markeredgecolor='g',
    markersize=marker_size,
    capsize=5
)
plt.errorbar(
    Rice95_v_many_coarse[x_exp],
    Rice95_v_many_coarse['ang_re'],
    yerr=Rice95_v_many_coarse['ang_re_std'],
    fmt='bs',
    markerfacecolor='none',
    markeredgecolor='b',
    markersize=marker_size,
    capsize=5
)
pwf = plt.plot(Willetts89_v_many_fine[x_exp], Willetts89_v_many_fine['ang_re'], 'ro', markersize=marker_size)
pwm = plt.plot(Willetts89_v_many_medium[x_exp], Willetts89_v_many_medium['ang_re'], 'g^', markersize=marker_size)
pwc = plt.plot(Willetts89_v_many_coarse[x_exp], Willetts89_v_many_coarse['ang_re'], 'bs', markersize=marker_size)
plt.tick_params(labelsize=ticks_size)

#plt.legend([pwf[0], pwm[0], pwc[0]], ['Fine', 'Medium', 'Coarse'], fontsize=label_size, loc='best')

plt.show()