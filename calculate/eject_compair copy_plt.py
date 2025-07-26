import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.stats import truncnorm
import matplotlib as mpl

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.usetex'] = True
label_size = 20
ticks_size = 15

def generate_truncated_lognormal(mu, sigma, dmin, dmax, size):
    a = (np.log(dmin) - mu) / sigma
    b = (np.log(dmax) - mu) / sigma
    samples = truncnorm.rvs(a, b, loc=mu, scale=sigma, size=size)
    return np.exp(samples)
#-------------------------------------------------------------------
output_rebound_ratio = False # output rebound ratio
output_e = True # output e
output_theta2 = True # output rebound angle
output_vn = True # output eject velocity
output_Nej = True # output eject number

d_min = 1.5e-4
d_max = 6e-4
mu = -8.30271
sigma = 0.25778
num_samples = 10000
#------------------------------------------------------------------
# average bed diameter
d2_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, num_samples)
d2_mid = np.percentile(d2_array, 50)
# impactor
coarse_d1 = 1.4*d2_mid
medium_d1 = d2_mid
fine_d1 = 0.73*d2_mid
## 对d2_array中所有介于0.00015和0.00025之间的值求平均值
#fine_d1 = np.mean(d2_array[(d2_array >= 0.00015) & (d2_array <= 0.00025)])
#medium_d1 = np.mean(d2_array[(d2_array >= 0.0003) & (d2_array <= 0.000355)])
#coarse_d1 = np.mean(d2_array[(d2_array >= 0.000425) & (d2_array <= 0.0006)])


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
Zhou06_v_many= {
    'ang_in': [8, 11.5],
    'ang_re_mean': [46.97, 47.53],
    'ang_re_std': [0.59, 0.75],
    'e_mean': [0.63, 0.61],
    'e_std': [0.006, 0.009]
}
Rice95_v_many_coarse= {
    'd_in': [coarse_d1, coarse_d1, coarse_d1, coarse_d1],
    'v_in': [2.7323, 2.7013, 2.7070, 2.7979],
    'ang_in': [13.94, 14.75, 14.73, 15.04],
    'ang_re': [20.95, 23.03, 22.55, 25.63],
    'e': [0.58, 0.56, 0.56, 0.58],
    'Nej': [5.7, 5.09, 5.12, 4.64],
    'v_ej': [0.243, 0.2422, 0.2522, 0.2637]
}
Rice95_v_many_medium= {
    'd_in': [medium_d1, medium_d1, medium_d1, medium_d1],
    'v_in': [3.1649, 3.3638, 3.3604, 3.2963],
    'ang_in': [11.82, 11.47, 11.53, 11.36],
    'ang_re': [29.95, 30.31, 31.06, 29.56],
    'e': [0.57, 0.54, 0.59, 0.58],
    'Nej': [2.68, 3.43, 2.67, 2.76],
    'v_ej': [0.2544, 0.252, 0.2607, 0.295]
}
Rice95_v_many_fine= {
    'd_in': [fine_d1, fine_d1, fine_d1],
    'v_in': [3.8493, 3.7801, 3.7411],
    'ang_in': [10.85, 10.24, 10.46],
    'ang_re': [44.63, 38.03, 37.85],
    'e': [0.52, 0.56, 0.58],
    'Nej': [1.86, 1.71, 1.58],
    'v_ej': [0.2757, 0.2566, 0.2773]
}
Chen18_v_many= {
    'ang_in': [23.2, 21.8, 21.9, 30.7, 30.6, 37.6, 47.1, 46.6, 46],
    'ang_re': [32.83, 35.59, 29.90, 43, 33.28, 44.16, 57.62, 50.97, 39.88],
    'e': [0.38, 0.4, 0.45, 0.32, 0.37, 0.27, 0.2, 0.19, 0.22]
}
Willetts89_v_many_coarse= {
    'd_in': [coarse_d1, coarse_d1, coarse_d1, coarse_d1],
    'v_in': [3.38, 3.43, 3.18, 3.50],
    'ang_in': [12.7, 17.8, 23.2, 27.7],
    'ang_re': [19.1, 25.2, 21.4, 27.2],
    'e': [0.63, 0.57, 0.54, 0.46],
    'v_ej': [0.3718, 0.3773, 0.4134, 0.42]
}
Willetts89_v_many_medium= {
    'd_in': [medium_d1, medium_d1, medium_d1, medium_d1],
    'v_in': [3.56, 3.99, 4.02, 4.39],
    'ang_in': [11.7, 18.2, 21.4, 26.3],
    'ang_re': [24.9, 33.4, 33.3, 44.7],
    'e': [0.61, 0.53, 0.50, 0.40],
    'v_ej': [0.356, 0.399, 0.3618, 0.439]
}
Willetts89_v_many_fine= {
    'd_in': [fine_d1, fine_d1, fine_d1, fine_d1],
    'v_in': [3.61, 4.41, 4.26, 4.35],
    'ang_in': [9.5, 15.4, 19.7, 24.9],
    'ang_re': [38.8, 42, 42.2, 42.5],
    'e': [0.57, 0.50, 0.48, 0.46],
    'v_ej': [0.3249, 0.3969, 0.3834, 0.348]
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

d1_array = Rice95_v_many_coarse['d_in'] + Rice95_v_many_medium['d_in'] + Rice95_v_many_fine['d_in']
d1_array = d1_array + Willetts89_v_many_coarse['d_in'] + Willetts89_v_many_medium['d_in'] + Willetts89_v_many_fine['d_in']
v1_array = Rice95_v_many_coarse['v_in'] + Rice95_v_many_medium['v_in'] + Rice95_v_many_fine['v_in']
v1_array = v1_array + Willetts89_v_many_coarse['v_in'] + Willetts89_v_many_medium['v_in'] + Willetts89_v_many_fine['v_in']
theta1_array = Rice95_v_many_coarse['ang_in'] + Rice95_v_many_medium['ang_in'] + Rice95_v_many_fine['ang_in']
theta1_array = theta1_array + Willetts89_v_many_coarse['ang_in'] + Willetts89_v_many_medium['ang_in'] + Willetts89_v_many_fine['ang_in']
iteration_num = len(d1_array)

e0_bar_list_0 = np.loadtxt('e0_bar_list_0_old.txt').tolist()
e0_bar_list_1 = np.loadtxt('e0_bar_list_1_old.txt').tolist()
e0_bar_list_2 = np.loadtxt('e0_bar_list_2_old.txt').tolist()
theta20_bar_list_0 = np.loadtxt('theta20_bar_list_0_old.txt').tolist()
theta20_bar_list_1 = np.loadtxt('theta20_bar_list_1_old.txt').tolist()
theta20_bar_list_2 = np.loadtxt('theta20_bar_list_2_old.txt').tolist()
Nej_bar_list_0 = np.loadtxt('Nej_bar_list_0_old.txt').tolist()
Nej_bar_list_1 = np.loadtxt('Nej_bar_list_1_old.txt').tolist()
Nej_bar_list_2 = np.loadtxt('Nej_bar_list_2_old.txt').tolist()
vn_bar_list_0 = np.loadtxt('vn_bar_list_0_old.txt').tolist()
vn_bar_list_1 = np.loadtxt('vn_bar_list_1_old.txt').tolist()
vn_bar_list_2 = np.loadtxt('vn_bar_list_2_old.txt').tolist()
len_R95_c = len(Rice95_v_many_coarse['e'])
len_R95_m = len(Rice95_v_many_medium['e'])
len_R95_f = len(Rice95_v_many_fine['e'])
len_W89_c = len(Willetts89_v_many_coarse['e'])
len_W89_m = len(Willetts89_v_many_medium['e'])
len_W89_f = len(Willetts89_v_many_fine['e'])
color_map = 'jet'
size_c = 80
size_m = 50
size_f = 20
if output_e:
    plt.figure(1, figsize=(8, 6))
    colors = Rice95_v_many_coarse['ang_in']
    x_array = Rice95_v_many_coarse['e']
    y_array = e0_bar_list_0[0:len_R95_c]
    len0 = len_R95_c
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='^', label='Ri95 Coarse')
    min_xy = min(x_array + y_array)
    max_xy = max(x_array + y_array)

    colors = Rice95_v_many_medium['ang_in']
    x_array = Rice95_v_many_medium['e']
    y_array = e0_bar_list_0[len0:len0 + len_R95_m]
    len0 += len_R95_m
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='^', label='Ri95 Medium')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    colors = Rice95_v_many_fine['ang_in']
    x_array = Rice95_v_many_fine['e']
    y_array = e0_bar_list_0[len0:len0 + len_R95_f]
    len0 += len_R95_f
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='^', label='Ri95 Fine')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    colors = Willetts89_v_many_coarse['ang_in']
    x_array = Willetts89_v_many_coarse['e']
    y_array = e0_bar_list_0[len0:len0 + len_W89_c]
    len0 += len_W89_c
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='o', label='Wi89 Coarse')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    colors = Willetts89_v_many_medium['ang_in']
    x_array = Willetts89_v_many_medium['e']
    y_array = e0_bar_list_0[len0:len0 + len_W89_m]
    len0 += len_W89_m
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='o', label='Wi89 medium')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    colors = Willetts89_v_many_fine['ang_in']
    x_array = Willetts89_v_many_fine['e']
    y_array = e0_bar_list_0[len0:len0 + len_W89_f]
    len0 += len_W89_f
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='o', label='Wi89 fine')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    line_points = [0.4, 0.65]
    plt.plot(line_points, line_points, 'k--', label='y=x')
    plt.xlim(0.4, 0.65)
    plt.ylim(0.4, 0.65)

    plt.xlabel('$\overline{e}^{\mathrm{exp}}$', fontsize=label_size)
    plt.ylabel('$\overline{e}^{\mathrm{mod}}$', fontsize=label_size)
    plt.tick_params(labelsize=ticks_size)
    cbar = plt.colorbar(label='$\\theta$ (degrees)', extend='both')
    cbar.set_label('$\\theta$ (degrees)', rotation=90, labelpad=ticks_size, fontsize=label_size)
    cbar.ax.tick_params(labelsize=ticks_size)

    ms_c = np.sqrt(size_c)
    ms_m = np.sqrt(size_m)
    ms_f = np.sqrt(size_f)
    circle_c = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_c, label='Wi89 Coarse')
    circle_m = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_m, label='Wi89 Medium')
    circle_f = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_f, label='Wi89 Fine')
    triangle_c = mlines.Line2D([], [], color='k', marker='^', linestyle='None', markersize=ms_c, label='Ri95 Coarse')
    triangle_m = mlines.Line2D([], [], color='k', marker='^', linestyle='None', markersize=ms_m, label='Ri95 Medium')
    triangle_f = mlines.Line2D([], [], color='k', marker='^', linestyle='None', markersize=ms_f, label='Ri95 Fine')
    plt.legend(handles=[circle_c, circle_m, circle_f, triangle_c, triangle_m, triangle_f], loc='best', fontsize=ticks_size)

if output_theta2:
    plt.figure(3, figsize=(8, 6))
    colors = Rice95_v_many_coarse['ang_in']
    x_array = Rice95_v_many_coarse['ang_re']
    y_array = np.degrees(theta20_bar_list_0[0:len_R95_c])
    y_array = y_array.tolist()
    len0 = len_R95_c
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='^', label='Rice95 coarse')
    min_xy = min(x_array + y_array)
    max_xy = max(x_array + y_array)

    colors = Rice95_v_many_medium['ang_in']
    x_array = Rice95_v_many_medium['ang_re']
    y_array = np.degrees(theta20_bar_list_0[len0:len0 + len_R95_m])
    y_array = y_array.tolist()
    len0 += len_R95_m
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='^', label='Rice95 medium')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    colors = Rice95_v_many_fine['ang_in']
    x_array = Rice95_v_many_fine['ang_re']
    y_array = np.degrees(theta20_bar_list_0[len0:len0 + len_R95_f])
    y_array = y_array.tolist()
    len0 += len_R95_f
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='^', label='Rice95 fine')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    colors = Willetts89_v_many_coarse['ang_in']
    x_array = Willetts89_v_many_coarse['ang_re']
    y_array = np.degrees(theta20_bar_list_0[len0:len0 + len_W89_c])
    y_array = y_array.tolist()
    len0 += len_W89_c
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='o', label='Willetts89 coarse')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    colors = Willetts89_v_many_medium['ang_in']
    x_array = Willetts89_v_many_medium['ang_re']
    y_array = np.degrees(theta20_bar_list_0[len0:len0 + len_W89_m])
    y_array = y_array.tolist()
    len0 += len_W89_m
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='o', label='Willetts89 medium')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    colors = Willetts89_v_many_fine['ang_in']
    x_array = Willetts89_v_many_fine['ang_re']
    y_array = np.degrees(theta20_bar_list_0[len0:len0 + len_W89_f])
    y_array = y_array.tolist()
    len0 += len_W89_f
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='o', label='Willetts89 fine')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    line_points = [15, 45]
    plt.plot(line_points, line_points, 'k--', label='y=x')
    plt.xlim(15, 45)
    plt.ylim(15, 45)

    plt.xlabel('$\overline{\\theta \'}^{\mathrm{exp}}$ (degrees)', fontsize=label_size)
    plt.ylabel('$\overline{\\theta \'}^{\mathrm{mod}}$ (degrees)', fontsize=label_size)
    plt.tick_params(labelsize=ticks_size)
    cbar = plt.colorbar(label='$\\theta$ (degrees)', extend='both')
    cbar.set_label('$\\theta$ (degrees)', rotation=90, labelpad=ticks_size, fontsize=label_size)
    cbar.ax.tick_params(labelsize=ticks_size)

    ms_c = np.sqrt(size_c)
    ms_m = np.sqrt(size_m)
    ms_f = np.sqrt(size_f)
    circle_c = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_c, label='Wi89 Coarse')
    circle_m = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_m, label='Wi89 Medium')
    circle_f = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_f, label='Wi89 Fine')
    triangle_c = mlines.Line2D([], [], color='k', marker='^', linestyle='None', markersize=ms_c, label='Ri95 Coarse')
    triangle_m = mlines.Line2D([], [], color='k', marker='^', linestyle='None', markersize=ms_m, label='Ri95 Medium')
    triangle_f = mlines.Line2D([], [], color='k', marker='^', linestyle='None', markersize=ms_f, label='Ri95 Fine')
    plt.legend(handles=[circle_c, circle_m, circle_f, triangle_c, triangle_m, triangle_f], loc='best', fontsize=ticks_size)

    plt.figure(4, figsize=(8, 6))
    x_array = Rice95_v_many_coarse['ang_re'] + Rice95_v_many_medium['ang_re'] + Rice95_v_many_fine['ang_re']
    x_array = x_array + Willetts89_v_many_coarse['ang_re'] + Willetts89_v_many_medium['ang_re'] + Willetts89_v_many_fine['ang_re']
    y_array1 = np.degrees(theta20_bar_list_1)
    y_array2 = np.degrees(theta20_bar_list_2)
    y_array1 = y_array1.tolist()
    y_array2 = y_array2.tolist()
    min_xy = min(x_array + y_array1 + y_array2)
    max_xy = max(x_array + y_array1 + y_array2)
    line_points = [min_xy, max_xy]
    plt.plot(line_points, line_points, 'k--', label='y=x')
    plt.scatter(x_array, y_array1, c='k', marker='s', label='d2=d3')
    plt.scatter(x_array, y_array2, c='k', marker='*', label='d1=d2=d3')
    plt.xlabel('$\\theta_{re,exp}$ (degrees)')
    plt.ylabel('$\\theta_{re,sim}$ (degrees)')
    plt.legend()

if output_Nej:
    plt.figure(5, figsize=(8, 6))
    colors = Rice95_v_many_coarse['ang_in']
    x_array = Rice95_v_many_coarse['Nej']
    y_array = Nej_bar_list_0[0:len_R95_c]
    len0 = len_R95_c
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='o', label='Rice95 coarse')
    min_xy = min(x_array + y_array)
    max_xy = max(x_array + y_array)

    colors = Rice95_v_many_medium['ang_in']
    x_array = Rice95_v_many_medium['Nej']
    y_array = Nej_bar_list_0[len0:len0 + len_R95_m]
    len0 += len_R95_m
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='o', label='Rice95 medium')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    colors = Rice95_v_many_fine['ang_in']
    x_array = Rice95_v_many_fine['Nej']
    y_array = Nej_bar_list_0[len0:len0 + len_R95_f]
    len0 += len_R95_f
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='o', label='Rice95 fine')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    line_points = [min_xy, max_xy]
    plt.plot(line_points, line_points, 'k--', label='y=x')

    plt.xlabel('$N_{ej,exp}$')
    plt.ylabel('$N_{ej,sim}$')
    #plt.xlim(min_xy, max_xy)
    #plt.ylim(min_xy, max_xy)
    plt.colorbar(label='Impact angle (degrees)')
    ms_c = np.sqrt(size_c)
    ms_m = np.sqrt(size_m)
    ms_f = np.sqrt(size_f)
    circle_c = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_c, label='Rice95 coarse')
    circle_m = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_m, label='Rice95 medium')
    circle_f = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_f, label='Rice95 fine')
    plt.legend(handles=[circle_c, circle_m, circle_f], loc='best')

if output_vn:
    plt.figure(6, figsize=(8, 6))
    colors = Rice95_v_many_coarse['ang_in']
    x_array = Rice95_v_many_coarse['v_ej']
    y_array = vn_bar_list_0[0:len_R95_c]
    len0 = len_R95_c
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='o', label='Rice95 coarse')
    min_xy = min(x_array + y_array)
    max_xy = max(x_array + y_array)

    colors = Rice95_v_many_medium['ang_in']
    x_array = Rice95_v_many_medium['v_ej']
    y_array = vn_bar_list_0[len0:len0 + len_R95_m]
    len0 += len_R95_m
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='o', label='Rice95 medium')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    colors = Rice95_v_many_fine['ang_in']
    x_array = Rice95_v_many_fine['v_ej']
    y_array = vn_bar_list_0[len0:len0 + len_R95_f]
    len0 += len_R95_f
    plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='o', label='Rice95 fine')
    min_xy = min(min(x_array + y_array), min_xy)
    max_xy = max(max(x_array + y_array), max_xy)

    #colors = Willetts89_v_many_coarse['ang_in']
    #x_array = Willetts89_v_many_coarse['v_ej']
    #y_array = vn_bar_list_0[len0:len0 + len_W89_c]
    #len0 += len_W89_c
    #plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_c, marker='^', label='Willetts89 coarse')
    #min_xy = min(min(x_array + y_array), min_xy)
    #max_xy = max(max(x_array + y_array), max_xy)

    #colors = Willetts89_v_many_medium['ang_in']
    #x_array = Willetts89_v_many_medium['v_ej']
    #y_array = vn_bar_list_0[len0:len0 + len_W89_m]
    #len0 += len_W89_m
    #plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_m, marker='^', label='Willetts89 medium')
    #min_xy = min(min(x_array + y_array), min_xy)
    #max_xy = max(max(x_array + y_array), max_xy)

    #colors = Willetts89_v_many_fine['ang_in']
    #x_array = Willetts89_v_many_fine['v_ej']
    #y_array = vn_bar_list_0[len0:len0 + len_W89_f]
    #len0 += len_W89_f
    #plt.scatter(x_array, y_array, c=colors, cmap=color_map, s=size_f, marker='^', label='Willetts89 fine')
    #min_xy = min(min(x_array + y_array), min_xy)
    #max_xy = max(max(x_array + y_array), max_xy)

    line_points = [min_xy, max_xy]
    plt.plot(line_points, line_points, 'k--', label='y=x')

    plt.xlabel('$v_{n,exp}$ (m/s)')
    plt.ylabel('$v_{n,sim}$ (m/s)')
    #plt.xlim(min_xy, max_xy)
    #plt.ylim(min_xy, max_xy)
    plt.colorbar(label='Impact angle (degrees)')

    ms_c = np.sqrt(size_c)
    ms_m = np.sqrt(size_m)
    ms_f = np.sqrt(size_f)
    circle_c = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_c, label='Rice95 coarse')
    circle_m = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_m, label='Rice95 medium')
    circle_f = mlines.Line2D([], [], color='k', marker='o', linestyle='None', markersize=ms_f, label='Rice95 fine')
    plt.legend(handles=[circle_c, circle_m, circle_f], loc='best')

plt.show()