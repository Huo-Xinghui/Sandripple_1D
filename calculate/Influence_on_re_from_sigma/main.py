import numpy as np
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from data_process import print_time
from get_data import get_model_data_array

# Set up matplotlib parameters
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.usetex'] = True
A = 1/0.6 # amplification factor
label_size = 15*A
ticks_size = 10*A
marker_size = 10*A
marker_size_in = 6*A
marker_width = 1
marker_width_in = 1.0
linewidth = 1.5

# Mechanical properties
epsilon = 0.6442
nu = -0.7857
# Calculation parameters
perform_calculations = True
thd_list = [15, 30, 45] # incident angle
case_num = 100  # number of cases
# Bed PSD parameters
dist_params = {
    'd_min': 0.5e-4,
    'd_max': 10e-4,
    'mu': -8.30271,
    'sampling_num': 100
}
bed_type = {
    'simple': True,
    'monodisperse': False,
    'three_D': False
}
sigma_min = 0.01
sigma_max = 5.0

# *****************************************************************

total_steps = 9
current_step = 0
start_time = time.time()
last_time = start_time

if perform_calculations:
    # Get model data
    th = thd_list[0] * np.pi / 180  # convert to radians
    sigma_array = np.linspace(sigma_min, sigma_max, case_num)  # standard deviation array

    phi_array_15, e_array_15 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    th = thd_list[1] * np.pi / 180  # convert to radians

    phi_array_30, e_array_30 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    th = thd_list[2] * np.pi / 180  # convert to radians

    phi_array_45, e_array_45 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    rslt_dict_simple = {
        'sigma': sigma_array,
        'phid_15': np.degrees(phi_array_15),
        'e_15': e_array_15,
        'phid_30': np.degrees(phi_array_30),
        'e_30': e_array_30,
        'phid_45': np.degrees(phi_array_45),
        'e_45': e_array_45,
    }
    np.savez('rb_vs_sigma_simple.npz', **rslt_dict_simple)

    bed_type['simple'] = False  # switch to polydisperse bed
    bed_type['three_D'] = True  # switch to 3D bed

    # Get model data
    th = thd_list[0] * np.pi / 180  # convert to radians

    phi_array_15, e_array_15 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    th = thd_list[1] * np.pi / 180  # convert to radians

    phi_array_30, e_array_30 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    th = thd_list[2] * np.pi / 180  # convert to radians

    phi_array_45, e_array_45 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    rslt_dict_3D = {
        'sigma': sigma_array,
        'phid_15': np.degrees(phi_array_15),
        'e_15': e_array_15,
        'phid_30': np.degrees(phi_array_30),
        'e_30': e_array_30,
        'phid_45': np.degrees(phi_array_45),
        'e_45': e_array_45,
    }
    np.savez('rb_vs_sigma_3D.npz', **rslt_dict_3D)

    bed_type['simple'] = True  # switch to simple bed
    bed_type['three_D'] = False  # switch to 3D bed
    bed_type['monodisperse'] = True  # switch to monodisperse bed

    # Get model data
    th = thd_list[0] * np.pi / 180  # convert to radians

    phi_array_15, e_array_15 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    th = thd_list[1] * np.pi / 180  # convert to radians

    phi_array_30, e_array_30 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    th = thd_list[2] * np.pi / 180  # convert to radians

    phi_array_45, e_array_45 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    rslt_dict_mono = {
        'sigma': sigma_array,
        'phid_15': np.degrees(phi_array_15),
        'e_15': e_array_15,
        'phid_30': np.degrees(phi_array_30),
        'e_30': e_array_30,
        'phid_45': np.degrees(phi_array_45),
        'e_45': e_array_45,
    }
    np.savez('rb_vs_sigma_mono.npz', **rslt_dict_mono)

# Load results
rslt_dict_simple = np.load('rb_vs_sigma_simple.npz')
rslt_dict_3D = np.load('rb_vs_sigma_3D.npz')
rslt_dict_mono = np.load('rb_vs_sigma_mono.npz')
rslt_dict_simple_d10 = np.load('rb_vs_sigma_simple_d10.npz')
rslt_dict_simple_d90 = np.load('rb_vs_sigma_simple_d90.npz')
rslt_dict_3D_d10 = np.load('rb_vs_sigma_3D_d10.npz')
rslt_dict_3D_d90 = np.load('rb_vs_sigma_3D_d90.npz')
rslt_dict_mono_d10 = np.load('rb_vs_sigma_mono_d10.npz')
rslt_dict_mono_d90 = np.load('rb_vs_sigma_mono_d90.npz')

# Draw rebound angle data
fig = plt.figure(1, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

ax.plot(rslt_dict_simple['sigma'], rslt_dict_simple['phid_15'], 'r--', linewidth=linewidth)
ax.plot(rslt_dict_simple['sigma'], rslt_dict_simple['phid_30'], 'g--', linewidth=linewidth)
ax.plot(rslt_dict_simple['sigma'], rslt_dict_simple['phid_45'], 'b--', linewidth=linewidth)
ax.plot(rslt_dict_3D['sigma'], rslt_dict_3D['phid_15'], 'r-', linewidth=linewidth)
ax.plot(rslt_dict_3D['sigma'], rslt_dict_3D['phid_30'], 'g-', linewidth=linewidth)
ax.plot(rslt_dict_3D['sigma'], rslt_dict_3D['phid_45'], 'b-', linewidth=linewidth)
ax.plot(rslt_dict_mono['sigma'], rslt_dict_mono['phid_15'], 'r:', linewidth=linewidth)
ax.plot(rslt_dict_mono['sigma'], rslt_dict_mono['phid_30'], 'g:', linewidth=linewidth)
ax.plot(rslt_dict_mono['sigma'], rslt_dict_mono['phid_45'], 'b:', linewidth=linewidth)

#ax.set_xlim(0, 85)
#ax.set_ylim(0, 110)
ax.set_xlabel('$\\sigma$ (m)', fontsize=label_size)
ax.set_ylabel('$\\overline{\\theta \'}$(degree)', fontsize=label_size)
ax.tick_params(axis='both', labelsize=ticks_size)
# 创建自定义图例句柄
line_simple = mlines.Line2D([], [], color='k', linestyle='--', linewidth=linewidth, label='2D-SP')
line_3D = mlines.Line2D([], [], color='k', linestyle='-',linewidth=linewidth, label='3D-SP')
ax.legend([line_simple, line_3D], ['2D-SP', '3D-SP'], fontsize=ticks_size, loc='best', frameon=True)

# Draw restitution coefficient data
fig = plt.figure(2, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

ax.plot(rslt_dict_simple['sigma'], rslt_dict_simple['e_15'], 'r--', linewidth=linewidth)
ax.plot(rslt_dict_simple['sigma'], rslt_dict_simple['e_30'], 'g--', linewidth=linewidth)
ax.plot(rslt_dict_simple['sigma'], rslt_dict_simple['e_45'], 'b--', linewidth=linewidth)
ax.plot(rslt_dict_3D['sigma'], rslt_dict_3D['e_15'], 'r-', linewidth=linewidth)
ax.plot(rslt_dict_3D['sigma'], rslt_dict_3D['e_30'], 'g-', linewidth=linewidth)
ax.plot(rslt_dict_3D['sigma'], rslt_dict_3D['e_45'], 'b-', linewidth=linewidth)
ax.plot(rslt_dict_mono['sigma'], rslt_dict_mono['e_15'], 'r:', linewidth=linewidth)
ax.plot(rslt_dict_mono['sigma'], rslt_dict_mono['e_30'], 'g:', linewidth=linewidth)
ax.plot(rslt_dict_mono['sigma'], rslt_dict_mono['e_45'], 'b:', linewidth=linewidth)

#ax.set_xlim(0, 85)
#ax.set_ylim(0.15, 0.8)
ax.set_xlabel('$\\sigma$ (m)', fontsize=label_size)
ax.set_ylabel('$\\overline{e}$', fontsize=label_size)
ax.tick_params(axis='x', labelsize=ticks_size)
ax.tick_params(axis='y', labelsize=ticks_size)
# 创建自定义图例句柄
line_15 = mlines.Line2D([], [], color='r', linestyle='-', linewidth=linewidth, label='$\\theta = 15^{\\circ}$')
line_30 = mlines.Line2D([], [], color='k', linestyle='-',linewidth=linewidth, label='$\\theta = 30^{\\circ}$')
line_45 = mlines.Line2D([], [], color='b', linestyle='-',linewidth=linewidth, label='$\\theta = 45^{\\circ}$')
ax.legend([line_15, line_30, line_45], ['$\\theta = 15^{\\circ}$', '$\\theta = 30^{\\circ}$', '$\\theta = 45^{\\circ}$'], fontsize=ticks_size, loc='best', frameon=True)

plt.show()
