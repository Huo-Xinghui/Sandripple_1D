import numpy as np
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.ndimage import gaussian_filter1d
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from data_process import print_time
from get_data import get_model_data_array, get_d_data_array

# Set up matplotlib parameters
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amssymb}'
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
perform_calculations = False
thd_list = [15, 30, 45] # incident angle
case_num = 50  # number of cases
# Bed PSD parameters
dist_params = {
    'd_min': 0.5e-4,
    'd_max': 50e-4,
    'mu': -8.30271,
    'sampling_num': 1000
}
bed_type = {
    'simple': True,
    'monodisperse': True,
    'three_D': False
}
sigma_min = 0.01
sigma_max = 1.0

# *****************************************************************

total_steps = 9
current_step = 0
start_time = time.time()
last_time = start_time

if perform_calculations:
    # Get d data
    sigma_array = np.linspace(sigma_min, sigma_max, case_num)  # standard deviation array
    th = thd_list[0] * np.pi / 180  # convert to radians
    d1_array_15, d2_array_15, d3_array_15 = get_d_data_array(th, sigma_array, dist_params)
    th = thd_list[1] * np.pi / 180  # convert to radians
    d1_array_30, d2_array_30, d3_array_30 = get_d_data_array(th, sigma_array, dist_params)
    th = thd_list[2] * np.pi / 180  # convert to radians
    d1_array_45, d2_array_45, d3_array_45 = get_d_data_array(th, sigma_array, dist_params)
    rslt_dict_simple = {
        'sigma': sigma_array,
        'd1_15': d1_array_15,
        'd2_15': d2_array_15,
        'd3_15': d3_array_15,
        'd1_30': d1_array_30,
        'd2_30': d2_array_30,
        'd3_30': d3_array_30,
        'd1_45': d1_array_45,
        'd2_45': d2_array_45,
        'd3_45': d3_array_45
    }
    np.savez('d_vs_sigma.npz', **rslt_dict_simple)

    # Get model data
    th = thd_list[0] * np.pi / 180  # convert to radians
    sigma_array = np.linspace(sigma_min, sigma_max, case_num)  # standard deviation array

    phi_array_15, e_array_15, ecx_array_15, ecz_array_15, ez_array_15, ex_array_15 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    th = thd_list[1] * np.pi / 180  # convert to radians

    phi_array_30, e_array_30, ecx_array_30, ecz_array_30, ez_array_30, ex_array_30 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    th = thd_list[2] * np.pi / 180  # convert to radians

    phi_array_45, e_array_45, ecx_array_45, ecz_array_45, ez_array_45, ex_array_45 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    rslt_dict_simple = {
        'sigma': sigma_array,
        'phid_15': np.degrees(phi_array_15),
        'e_15': e_array_15,
        'ecx_15': ecx_array_15,
        'ecz_15': ecz_array_15,
        'ez_15': ez_array_15,
        'ex_15': ex_array_15,
        'phid_30': np.degrees(phi_array_30),
        'e_30': e_array_30,
        'ecx_30': ecx_array_30,
        'ecz_30': ecz_array_30,
        'ez_30': ez_array_30,
        'ex_30': ex_array_30,
        'phid_45': np.degrees(phi_array_45),
        'e_45': e_array_45,
        'ecx_45': ecx_array_45,
        'ecz_45': ecz_array_45,
        'ez_45': ez_array_45,
        'ex_45': ex_array_45,
    }
    #rslt_dict_simple = {
    #    'sigma': sigma_array,
    #    'phid_10': np.degrees(phi_array_15),
    #    'e_10': e_array_15,
    #    'ecx_10': ecx_array_15,
    #    'ecz_10': ecz_array_15,
    #    'ez_10': ez_array_15,
    #    'ex_10': ex_array_15,
    #    'phid_20': np.degrees(phi_array_30),
    #    'e_20': e_array_30,
    #    'ecx_20': ecx_array_30,
    #    'ecz_20': ecz_array_30,
    #    'ez_20': ez_array_30,
    #    'ex_20': ex_array_30,
    #    'phid_25': np.degrees(phi_array_45),
    #    'e_25': e_array_45,
    #    'ecx_25': ecx_array_45,
    #    'ecz_25': ecz_array_45,
    #    'ez_25': ez_array_45,
    #    'ex_25': ex_array_45,
    #}
    np.savez('rb_vs_sigma_simple_mono_d2.npz', **rslt_dict_simple)

    bed_type['simple'] = True  # switch to polydisperse bed
    bed_type['monodisperse'] = False  # switch to polydisperse bed
    bed_type['three_D'] = False  # switch to 3D bed

    # Get model data
    th = thd_list[0] * np.pi / 180  # convert to radians

    phi_array_15, e_array_15, ecx_array_15, ecz_array_15, ez_array_15, ex_array_15 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    th = thd_list[1] * np.pi / 180  # convert to radians

    phi_array_30, e_array_30, ecx_array_30, ecz_array_30, ez_array_30, ex_array_30 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    th = thd_list[2] * np.pi / 180  # convert to radians

    phi_array_45, e_array_45, ecx_array_45, ecz_array_45, ez_array_45, ex_array_45 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    rslt_dict_3D = {
        'sigma': sigma_array,
        'phid_15': np.degrees(phi_array_15),
        'e_15': e_array_15,
        'ecx_15': ecx_array_15,
        'ecz_15': ecz_array_15,
        'ez_15': ez_array_15,
        'ex_15': ex_array_15,
        'phid_30': np.degrees(phi_array_30),
        'e_30': e_array_30,
        'ecx_30': ecx_array_30,
        'ecz_30': ecz_array_30,
        'ez_30': ez_array_30,
        'ex_30': ex_array_30,
        'phid_45': np.degrees(phi_array_45),
        'e_45': e_array_45,
        'ecx_45': ecx_array_45,
        'ecz_45': ecz_array_45,
        'ez_45': ez_array_45,
        'ex_45': ex_array_45,
    }
    np.savez('rb_vs_sigma_simple_mono_d3.npz', **rslt_dict_3D)

    #bed_type['simple'] = True  # switch to simple bed
    #bed_type['monodisperse'] = True  # switch to monodisperse bed
    #bed_type['three_D'] = False  # switch to 3D bed

    ## Get model data
    #th = thd_list[0] * np.pi / 180  # convert to radians

    #phi_array_15, e_array_15, ecx_array_15, ecz_array_15, ez_array_15, ex_array_15 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    #current_time = time.time()
    #current_step += 1
    #print_time(current_time, last_time, start_time, current_step, total_steps)
    #last_time = current_time

    #th = thd_list[1] * np.pi / 180  # convert to radians

    #phi_array_30, e_array_30, ecx_array_30, ecz_array_30, ez_array_30, ex_array_30 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    #current_time = time.time()
    #current_step += 1
    #print_time(current_time, last_time, start_time, current_step, total_steps)
    #last_time = current_time

    #th = thd_list[2] * np.pi / 180  # convert to radians

    #phi_array_45, e_array_45, ecx_array_45, ecz_array_45, ez_array_45, ex_array_45 = get_model_data_array(epsilon, nu, th, sigma_array, dist_params, bed_type)
    #current_time = time.time()
    #current_step += 1
    #print_time(current_time, last_time, start_time, current_step, total_steps)
    #last_time = current_time

    #rslt_dict_mono = {
    #    'sigma': sigma_array,
    #    'phid_15': np.degrees(phi_array_15),
    #    'e_15': e_array_15,
    #    'ecx_15': ecx_array_15,
    #    'ecz_15': ecz_array_15,
    #    'ez_15': ez_array_15,
    #    'ex_15': ex_array_15,
    #    'phid_30': np.degrees(phi_array_30),
    #    'e_30': e_array_30,
    #    'ecx_30': ecx_array_30,
    #    'ecz_30': ecz_array_30,
    #    'ez_30': ez_array_30,
    #    'ex_30': ex_array_30,
    #    'phid_45': np.degrees(phi_array_45),
    #    'e_45': e_array_45,
    #    'ecx_45': ecx_array_45,
    #    'ecz_45': ecz_array_45,
    #    'ez_45': ez_array_45,
    #    'ex_45': ex_array_45,
    #}
    #np.savez('rb_vs_sigma_mono_dm.npz', **rslt_dict_mono)

# Load results
rslt_dict_2D = np.load('rb_vs_sigma_2D_poly.npz')
rslt_dict_3D = np.load('rb_vs_sigma_3D_poly.npz')
rslt_dict_d50 = np.load('rb_vs_sigma_mono_d50.npz')
rslt_dict_d90 = np.load('rb_vs_sigma_mono_d90.npz')
rslt_dict_dm = np.load('rb_vs_sigma_mono_dm.npz')
rslt_dict_monod2 = np.load('rb_vs_sigma_monod2.npz')
rslt_dict_monod3 = np.load('rb_vs_sigma_monod3.npz')
rslt_dict_smonod2 = np.load('rb_vs_sigma_simple_mono_d2.npz')
rslt_dict_smonod3 = np.load('rb_vs_sigma_simple_mono_d3.npz')
rslt_dict_3D_ex = np.load('rb_vs_sigma_3D_ex.npz')
rslt_dict_dm_ex = np.load('rb_vs_sigma_mono_dm_ex.npz')
rslt_dict_d90_ex = np.load('rb_vs_sigma_mono_d90_ex.npz')
d_dict = np.load('d_vs_sigma.npz')

skip = 4

ez_array_2D_15 = rslt_dict_2D['ez_15'][::skip]
ez_array_2D_30 = rslt_dict_2D['ez_30'][::skip]
ez_array_2D_45 = rslt_dict_2D['ez_45'][::skip]
ez_array_3D_15 = rslt_dict_3D['ez_15'][::skip]
ez_array_3D_30 = rslt_dict_3D['ez_30'][::skip]
ez_array_3D_45 = rslt_dict_3D['ez_45'][::skip]
ez_array_d2_15 = rslt_dict_monod2['ez_15'][::skip]
ez_array_d2_30 = rslt_dict_monod2['ez_30'][::skip]
ez_array_d2_45 = rslt_dict_monod2['ez_45'][::skip]
ez_array_d3_15 = rslt_dict_monod3['ez_15'][::skip]
ez_array_d3_30 = rslt_dict_monod3['ez_30'][::skip]
ez_array_d3_45 = rslt_dict_monod3['ez_45'][::skip]
ez_array_3D_ex_10 = rslt_dict_3D_ex['ez_10'][::skip]
ez_array_3D_ex_20 = rslt_dict_3D_ex['ez_20'][::skip]
ez_array_3D_ex_25 = rslt_dict_3D_ex['ez_25'][::skip]

ex_array_2D_15 = rslt_dict_2D['ex_15'][::skip]
ex_array_2D_30 = rslt_dict_2D['ex_30'][::skip]
ex_array_2D_45 = rslt_dict_2D['ex_45'][::skip]
ex_array_3D_15 = rslt_dict_3D['ex_15'][::skip]
ex_array_3D_30 = rslt_dict_3D['ex_30'][::skip]
ex_array_3D_45 = rslt_dict_3D['ex_45'][::skip]
ex_array_d2_15 = rslt_dict_monod2['ex_15'][::skip]
ex_array_d2_30 = rslt_dict_monod2['ex_30'][::skip]
ex_array_d2_45 = rslt_dict_monod2['ex_45'][::skip]
ex_array_d3_15 = rslt_dict_monod3['ex_15'][::skip]
ex_array_d3_30 = rslt_dict_monod3['ex_30'][::skip]
ex_array_d3_45 = rslt_dict_monod3['ex_45'][::skip]
ex_array_3D_ex_10 = rslt_dict_3D_ex['ex_10'][::skip]
ex_array_3D_ex_20 = rslt_dict_3D_ex['ex_20'][::skip]
ex_array_3D_ex_25 = rslt_dict_3D_ex['ex_25'][::skip]

sigma_array_2D = rslt_dict_2D['sigma'][::skip]
sigma_array_3D = rslt_dict_3D['sigma'][::skip]
sigma_array_d2 = rslt_dict_monod2['sigma'][::skip]

# smooth ez
#ez_array_2D_15 = gaussian_filter1d(ez_array_2D_15, sigma=1)
#ez_array_2D_30 = gaussian_filter1d(ez_array_2D_30, sigma=1)
#ez_array_2D_45 = gaussian_filter1d(ez_array_2D_45, sigma=1)
#ez_array_3D_15 = gaussian_filter1d(ez_array_3D_15, sigma=1)
#ez_array_3D_30 = gaussian_filter1d(ez_array_3D_30, sigma=1)
#ez_array_3D_45 = gaussian_filter1d(ez_array_3D_45, sigma=1)

# Draw rebound angle data
fig = plt.figure(1, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

th = 15/180*np.pi
d1_15 = d_dict['d1_15']
d2_15 = d_dict['d2_15']
d3_15 = d_dict['d3_15']
d13_15 = (d1_15 + d3_15)/(d1_15 + d2_15)
d23_15 = (d2_15 + d3_15)/(d1_15 + d2_15)
xmin_15 = d13_15/np.sin(th) - d23_15

th = 30/180*np.pi
d1_30 = d_dict['d1_30']
d2_30 = d_dict['d2_30']
d3_30 = d_dict['d3_30']
d13_30 = (d1_30 + d3_30)/(d1_30 + d2_30)
d23_30 = (d2_30 + d3_30)/(d1_30 + d2_30)
xmin_30 = d13_30/np.sin(th) - d23_30

th = 45/180*np.pi
d1_45 = d_dict['d1_45']
d2_45 = d_dict['d2_45']
d3_45 = d_dict['d3_45']
d13_45 = (d1_45 + d3_45)/(d1_45 + d2_45)
d23_45 = (d2_45 + d3_45)/(d1_45 + d2_45)
xmin_45 = d13_45/np.sin(th) - d23_45

#ax.plot(d_dict['sigma'], d3_15, 'r--')
#ax.plot(d_dict['sigma'], d3_30, 'g--')
#ax.plot(d_dict['sigma'], d3_45, 'b--')

ax.plot(d_dict['sigma'], xmin_15, 'r-')
ax.plot(d_dict['sigma'], xmin_30, 'g-')
ax.plot(d_dict['sigma'], xmin_45, 'b-')

#ax.set_xlim(0, 85)
#ax.set_ylim(0, 110)
ax.set_xlabel('$d_2$ (m)', fontsize=label_size)
ax.set_ylabel('$d_3$ (m)', fontsize=label_size)
ax.tick_params(axis='both', labelsize=ticks_size)
## 创建自定义图例句柄
#line_simple = mlines.Line2D([], [], color='k', linestyle='--', linewidth=linewidth, label='2D-SP')
#line_3D = mlines.Line2D([], [], color='k', linestyle='-',linewidth=linewidth, label='3D-SP')
#ax.legend([line_simple, line_3D], ['2D-SP', '3D-SP'], fontsize=ticks_size, loc='best', frameon=True)

# Draw restitution coefficient data
fig = plt.figure(2, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

ax.plot(rslt_dict_2D['sigma'], rslt_dict_smonod3['e_15'], 'r--', linewidth=linewidth)
ax.plot(rslt_dict_2D['sigma'], rslt_dict_smonod3['e_30'], 'g--', linewidth=linewidth)
ax.plot(rslt_dict_2D['sigma'], rslt_dict_smonod3['e_45'], 'b--', linewidth=linewidth)
ax.plot(rslt_dict_3D['sigma'], rslt_dict_3D['e_15'], 'r-', linewidth=linewidth)
ax.plot(rslt_dict_3D['sigma'], rslt_dict_3D['e_30'], 'g-', linewidth=linewidth)
ax.plot(rslt_dict_3D['sigma'], rslt_dict_3D['e_45'], 'b-', linewidth=linewidth)
ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm['e_15'], 'r:', linewidth=linewidth)
ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm['e_30'], 'g:', linewidth=linewidth)
ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm['e_45'], 'b:', linewidth=linewidth)

#ax.set_xlim(0, 85)
#ax.set_ylim(0.15, 0.8)
ax.set_xlabel('$\\sigma$', fontsize=label_size)
ax.set_ylabel('$\\overline{e}$', fontsize=label_size)
ax.tick_params(axis='x', labelsize=ticks_size)
ax.tick_params(axis='y', labelsize=ticks_size)
# 创建自定义图例句柄
line_15 = mlines.Line2D([], [], color='r', linestyle='-', linewidth=linewidth, label='$\\theta = 15^{\\circ}$')
line_30 = mlines.Line2D([], [], color='k', linestyle='-',linewidth=linewidth, label='$\\theta = 30^{\\circ}$')
line_45 = mlines.Line2D([], [], color='b', linestyle='-',linewidth=linewidth, label='$\\theta = 45^{\\circ}$')
ax.legend([line_15, line_30, line_45], ['$\\theta = 15^{\\circ}$', '$\\theta = 30^{\\circ}$', '$\\theta = 45^{\\circ}$'], fontsize=ticks_size, loc='best', frameon=True)

# Draw ez data
fig = plt.figure(3, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

ax.plot(sigma_array_2D, ez_array_2D_15, 'ro', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
#ax.plot(sigma_array_2D, ez_array_2D_30, 'g^', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(sigma_array_2D, ez_array_2D_45, 'bs', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(sigma_array_3D, ez_array_3D_15, 'ro', markersize=marker_size, markeredgewidth=marker_width)
ax.plot(sigma_array_3D, ez_array_3D_ex_20, 'c^', markersize=marker_size, markeredgewidth=marker_width)
ax.plot(sigma_array_3D, ez_array_3D_30, 'yh', markersize=marker_size, markeredgewidth=marker_width)
ax.plot(sigma_array_3D, ez_array_3D_45, 'bs', markersize=marker_size, markeredgewidth=marker_width)
ax.plot(rslt_dict_d50['sigma'], rslt_dict_d50['ez_15'], 'r-', linewidth=linewidth)
#ax.plot(rslt_dict_d50['sigma'], rslt_dict_d50['ez_30'], 'g-', linewidth=linewidth)
ax.plot(rslt_dict_d50['sigma'], rslt_dict_d50['ez_45'], 'b-', linewidth=linewidth)
ax.plot(rslt_dict_d90['sigma'], rslt_dict_d90['ez_15'], 'r--', linewidth=linewidth)
#ax.plot(rslt_dict_d90['sigma'], rslt_dict_d90['ez_30'], 'g--', linewidth=linewidth)
ax.plot(rslt_dict_d90['sigma'], rslt_dict_d90['ez_45'], 'b--', linewidth=linewidth)
ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm['ez_15'], 'k-', linewidth=linewidth*0.5)
ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm_ex['ez_20'], 'k-', linewidth=linewidth*0.5)
ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm['ez_30'], 'k-', linewidth=linewidth*0.5)
ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm['ez_45'], 'k-', linewidth=linewidth*0.5)
#ax.plot(sigma_array_3D, ez_array_d2_15, 'r*', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
#ax.plot(sigma_array_3D, ez_array_d2_45, 'b*', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
#ax.plot(sigma_array_3D, ez_array_d3_15, 'r^', markersize=marker_size, markeredgewidth=marker_width)
#ax.plot(sigma_array_3D, ez_array_d3_45, 'b^', markersize=marker_size, markeredgewidth=marker_width)


ax.set_xlim(0, 1)
ax.set_ylim(0.1, 2)
ax.set_xlabel('$\\sigma_d$', fontsize=label_size)
ax.set_ylabel('$\\overline{e_z}$', fontsize=label_size)
ax.tick_params(axis='x', labelsize=ticks_size)
ax.tick_params(axis='y', labelsize=ticks_size)

inset_ax = inset_axes(ax, width='45%', height='45%', loc='upper left')
inset_ax.plot(sigma_array_d2, ez_array_d2_15, 'mx-', markersize=marker_size_in, markerfacecolor='none', markeredgewidth=marker_width_in)
inset_ax.plot(sigma_array_d2, ez_array_d2_45, 'gx', markersize=marker_size_in, markerfacecolor='none', markeredgewidth=marker_width_in)
inset_ax.plot(sigma_array_3D, ez_array_3D_15, 'ro', markersize=marker_size_in, markeredgewidth=marker_width_in)
inset_ax.plot(sigma_array_3D, ez_array_3D_45, 'bs', markersize=marker_size_in, markeredgewidth=marker_width_in)
#inset_ax.plot(rslt_dict_smonod2['sigma'], rslt_dict_smonod2['ez_15'], 'm--', linewidth=linewidth)
#inset_ax.plot(rslt_dict_smonod2['sigma'], rslt_dict_smonod2['ez_45'], 'g--', linewidth=linewidth)
inset_ax.plot(rslt_dict_smonod3['sigma'], rslt_dict_smonod3['ez_15'], 'm-', linewidth=linewidth)
inset_ax.plot(rslt_dict_smonod3['sigma'], rslt_dict_smonod3['ez_45'], 'g-', linewidth=linewidth)
inset_ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm['ez_15'], 'k-', linewidth=linewidth*0.5)
inset_ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm['ez_45'], 'k-', linewidth=linewidth*0.5)
inset_ax.set_xlim(0, 1)
inset_ax.set_ylim(0.0, 1.2)
#inset_ax.tick_params(axis='x', labelsize=ticks_size)
#inset_ax.tick_params(axis='y', labelsize=ticks_size)
## 把y的tick移到右边
inset_ax.yaxis.tick_right()
inset_ax.yaxis.set_label_position("right")
inset_ax.tick_params(axis='both', length=0)
inset_ax.set_xticks([])
inset_ax.set_yticks([])
# 创建自定义图例句柄
solid_line_m = inset_ax.plot([], [], 'm-', linewidth=linewidth)[0]
x_line_m = inset_ax.plot([], [], 'mx-', linewidth=linewidth)[0]
solid_line_g = inset_ax.plot([], [], 'g-', linewidth=linewidth)[0]
x_line_g = inset_ax.plot([], [], 'gx-', linewidth=linewidth)[0]
inset_ax.legend([solid_line_m, x_line_m, solid_line_g, x_line_g],
                ['$\\theta = 15^{\\circ}$, $d_3 \\equiv d_{50}$', '$\\theta = 15^{\\circ}$, Sampled $d_3$', '$\\theta = 45^{\\circ}$, $d_3 \\equiv d_{50}$', '$\\theta = 45^{\\circ}$, Sampled $d_3$'],
                fontsize=ticks_size,
                loc='upper right',
                bbox_to_anchor=(2.0, 1.06),
                frameon=True)


# Draw ex data
fig = plt.figure(4, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

ax.plot(sigma_array_2D, ex_array_2D_15, 'ro', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(sigma_array_2D, ex_array_2D_30, 'yh', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(sigma_array_2D, ex_array_2D_45, 'bs', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(sigma_array_3D, ex_array_3D_15, 'ro', markersize=marker_size, markeredgewidth=marker_width)
ax.plot(sigma_array_3D, ex_array_3D_30, 'yh', markersize=marker_size, markeredgewidth=marker_width)
ax.plot(sigma_array_3D, ex_array_3D_45, 'bs', markersize=marker_size, markeredgewidth=marker_width)
ax.plot(rslt_dict_d50['sigma'], rslt_dict_d50['ex_15'], 'r-', linewidth=linewidth)
ax.plot(rslt_dict_d50['sigma'], rslt_dict_d50['ex_30'], 'y-', linewidth=linewidth)
ax.plot(rslt_dict_d50['sigma'], rslt_dict_d50['ex_45'], 'b-', linewidth=linewidth)
#ax.plot(rslt_dict_d50['sigma'], rslt_dict_d90_ex['ex_10'], 'c--', linewidth=linewidth)
ax.plot(rslt_dict_d90['sigma'], rslt_dict_d90['ex_15'], 'r--', linewidth=linewidth)
ax.plot(rslt_dict_d90['sigma'], rslt_dict_d90['ex_30'], 'y--', linewidth=linewidth)
ax.plot(rslt_dict_d90['sigma'], rslt_dict_d90['ex_45'], 'b--', linewidth=linewidth)
ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm['ex_15'], 'k-', linewidth=linewidth*0.5)
ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm['ex_30'], 'k-', linewidth=linewidth*0.5)
ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm['ex_45'], 'k-', linewidth=linewidth*0.5)
#ax.plot(sigma_array_3D, ex_array_d2_15, 'r*', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
#ax.plot(sigma_array_3D, ex_array_d2_45, 'b*', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
#ax.plot(sigma_array_3D, ex_array_d3_15, 'r^', markersize=marker_size, markeredgewidth=marker_width)
#ax.plot(sigma_array_3D, ex_array_d3_45, 'b^', markersize=marker_size, markeredgewidth=marker_width)

ax.set_xlim(0, 1)
ax.set_ylim(0.0, 0.6)
ax.set_xlabel('$\\sigma_d$', fontsize=label_size)
ax.set_ylabel('$\\overline{e_x}$', fontsize=label_size)
ax.tick_params(axis='x', labelsize=ticks_size)
ax.tick_params(axis='y', labelsize=ticks_size)
# 创建自定义图例句柄
thin_line = mlines.Line2D([], [], color='k', linestyle='-', linewidth=linewidth*0.5)
dashed_line = mlines.Line2D([], [], color='k', linestyle='-', linewidth=linewidth)
dotted_line = mlines.Line2D([], [], color='k', linestyle='--', linewidth=linewidth)
ax.legend([thin_line, dashed_line, dotted_line],
          ['$d = \\mathbb{E}[d]$', '$d = d_{50}$', '$d = d_{90}$'],
          fontsize=ticks_size,
          loc='lower right',
          bbox_to_anchor=(1.01, 0.0),
          frameon=True)

inset_ax = inset_axes(ax, width='45%', height='45%', loc='lower left')
inset_ax.plot(sigma_array_d2, ex_array_d2_15, 'mx-', markersize=marker_size_in, markerfacecolor='none', markeredgewidth=marker_width_in)
inset_ax.plot(sigma_array_d2, ex_array_d2_45, 'gx-', markersize=marker_size_in, markerfacecolor='none', markeredgewidth=marker_width_in)
inset_ax.plot(sigma_array_3D, ex_array_3D_15, 'ro', markersize=marker_size_in, markeredgewidth=marker_width_in)
inset_ax.plot(sigma_array_3D, ex_array_3D_45, 'bs', markersize=marker_size_in, markeredgewidth=marker_width_in)
#inset_ax.plot(rslt_dict_smonod2['sigma'], rslt_dict_smonod2['ex_15'], 'm--', linewidth=linewidth)
#inset_ax.plot(rslt_dict_smonod2['sigma'], rslt_dict_smonod2['ex_45'], 'g--', linewidth=linewidth)
inset_ax.plot(rslt_dict_smonod3['sigma'], rslt_dict_smonod3['ex_15'], 'm-', linewidth=linewidth)
inset_ax.plot(rslt_dict_smonod3['sigma'], rslt_dict_smonod3['ex_45'], 'g-', linewidth=linewidth)
inset_ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm['ex_15'], 'k-', linewidth=linewidth*0.5)
inset_ax.plot(rslt_dict_dm['sigma'], rslt_dict_dm['ex_45'], 'k-', linewidth=linewidth*0.5)
inset_ax.set_xlim(0, 1)
inset_ax.set_ylim(0.2, 0.6)
#inset_ax.tick_params(axis='x', labelsize=ticks_size)
#inset_ax.tick_params(axis='y', labelsize=ticks_size)
## 把y的tick移到右边
inset_ax.yaxis.tick_right()
inset_ax.yaxis.set_label_position("right")
inset_ax.tick_params(axis='both', length=0)
inset_ax.set_xticks([])
inset_ax.set_yticks([])
# 创建自定义图例句柄
symbol_15 = mlines.Line2D([], [], color='r', marker='o', markersize=marker_size, markeredgewidth=marker_width)
symbol_20 = mlines.Line2D([], [], color='c', marker='^', markersize=marker_size, markeredgewidth=marker_width)
symbol_30 = mlines.Line2D([], [], color='y', marker='h', markersize=marker_size, markeredgewidth=marker_width)
symbol_45 = mlines.Line2D([], [], color='b', marker='s', markersize=marker_size, markeredgewidth=marker_width)
inset_ax.legend([symbol_15, symbol_20, symbol_30, symbol_45],
          ['$\\theta = 15^{\\circ}$', '$\\theta = 20^{\\circ}$', '$\\theta = 30^{\\circ}$', '$\\theta = 45^{\\circ}$'],
          fontsize=ticks_size,
          loc='upper right',
          bbox_to_anchor=(1.66, 0.65),
          frameon=True)

plt.show()
