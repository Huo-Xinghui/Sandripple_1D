import numpy as np
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.ndimage import gaussian_filter1d
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from data_process import print_time
from get_data import get_model_data_array

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
rho = 2650
rhof = 1.263
s = rho/rhof
g = 9.8
physical_dict = {
    'epsilon': 0.6821,
    'nu': -0.8606,
    'rho': rho,
    'g': g,
    'gamma': 0.042
}
# Calculation parameters
perform_calculations = False
v1_list = [30, 50, 100] # incident angle
thd = 15
th = thd * np.pi / 180  # convert to radians
case_num = 50  # number of cases
# Bed PSD parameters
dist_params = {
    'd_min': 0.5e-4,
    'd_max': 50e-4,
    'mu': -8.30271,
    'sigma': 0.5,
    'sampling_num': 1000
}
bed_type = {
    'monodisperse': True,
    'zc': 0,
    'd2': 0
}
sigma_min = 0.01
sigma_max = 1.0

# *****************************************************************

total_steps = 9
current_step = 0
start_time = time.time()
last_time = start_time

if perform_calculations:
    # Get model data
    v1 = v1_list[0]
    sigma_array = np.linspace(sigma_min, sigma_max, case_num)  # standard deviation array

    N_array_0, v_array_0, E_array_0, E1_array_0 = get_model_data_array(v1, th, sigma_array, dist_params, bed_type, physical_dict)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    v1 = v1_list[1]

    N_array_1, v_array_1, E_array_1, E1_array_1 = get_model_data_array(v1, th, sigma_array, dist_params, bed_type, physical_dict)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    v1 = v1_list[2]

    N_array_2, v_array_2, E_array_2, E1_array_2 = get_model_data_array(v1, th, sigma_array, dist_params, bed_type, physical_dict)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    rslt_dict_simple = {
        'sigma': sigma_array,
        'Nej0': N_array_0,
        'vej0': v_array_0,
        'Eej0': E_array_0,
        'E10': E1_array_0,
        'Nej1': N_array_1,
        'vej1': v_array_1,
        'Eej1': E_array_1,
        'E11': E1_array_1,
        'Nej2': N_array_2,
        'vej2': v_array_2,
        'Eej2': E_array_2,
        'E12': E1_array_2
    }
    np.savez('ej_vs_sigma_mono.npz', **rslt_dict_simple)

    bed_type['monodisperse'] = False  # switch to polydisperse bed

    # Get model data
    v1 = v1_list[0]
    sigma_array = np.linspace(sigma_min, sigma_max, case_num)  # standard deviation array

    N_array_0, v_array_0, E_array_0, E1_array_0 = get_model_data_array(v1, th, sigma_array, dist_params, bed_type, physical_dict)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    v1 = v1_list[1]

    N_array_1, v_array_1, E_array_1, E1_array_1 = get_model_data_array(v1, th, sigma_array, dist_params, bed_type, physical_dict)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    v1 = v1_list[2]

    N_array_2, v_array_2, E_array_2, E1_array_2 = get_model_data_array(v1, th, sigma_array, dist_params, bed_type, physical_dict)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    rslt_dict_simple = {
        'sigma': sigma_array,
        'Nej0': N_array_0,
        'vej0': v_array_0,
        'Eej0': E_array_0,
        'E10': E1_array_0,
        'Nej1': N_array_1,
        'vej1': v_array_1,
        'Eej1': E_array_1,
        'E11': E1_array_1,
        'Nej2': N_array_2,
        'vej2': v_array_2,
        'Eej2': E_array_2,
        'E12': E1_array_2
    }
    np.savez('ej_vs_sigma.npz', **rslt_dict_simple)

# Load results
rslt_dict_mono = np.load('ej_vs_sigma_mono.npz')
rslt_dict = np.load('ej_vs_sigma.npz')

skip = 4

sigma_skp = rslt_dict['sigma'][::skip]
N_array_0 = rslt_dict['Nej0'][::skip]
N_array_1 = rslt_dict['Nej1'][::skip]
N_array_2 = rslt_dict['Nej2'][::skip]
N_array_m_0 = rslt_dict_mono['Nej0'][::skip]
N_array_m_1 = rslt_dict_mono['Nej1'][::skip]
N_array_m_2 = rslt_dict_mono['Nej2'][::skip]
v_array_0 = rslt_dict['vej0'][::skip]
v_array_1 = rslt_dict['vej1'][::skip]
v_array_2 = rslt_dict['vej2'][::skip]
E_array_0 = rslt_dict['Eej0'][::skip]
E_array_1 = rslt_dict['Eej1'][::skip]
E_array_2 = rslt_dict['Eej2'][::skip]
E_array_m_0 = rslt_dict_mono['Eej0'][::skip]
E_array_m_1 = rslt_dict_mono['Eej1'][::skip]
E_array_m_2 = rslt_dict_mono['Eej2'][::skip]
E1_array_0 = rslt_dict['E10'][::skip]
E1_array_1 = rslt_dict['E11'][::skip]
E1_array_2 = rslt_dict['E12'][::skip]
Et_array_0 = E_array_0 * N_array_0
Et_array_1 = E_array_1 * N_array_1
Et_array_2 = E_array_2 * N_array_2
Et_array_m_0 = E_array_m_0 * N_array_m_0
Et_array_m_1 = E_array_m_1 * N_array_m_1
Et_array_m_2 = E_array_m_2 * N_array_m_2
Et_array_mono_0 = rslt_dict_mono['Eej0']*rslt_dict_mono['Nej0']
Et_array_mono_1 = rslt_dict_mono['Eej1']*rslt_dict_mono['Nej1']
Et_array_mono_2 = rslt_dict_mono['Eej2']*rslt_dict_mono['Nej2']

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

ax.plot(sigma_skp, E_array_0/E1_array_0, 'ro', markersize=marker_size)
ax.plot(sigma_skp, E_array_1/E1_array_1, 'g^', markersize=marker_size)
ax.plot(sigma_skp, E_array_2/E1_array_2, 'bs', markersize=marker_size)
ax.plot(rslt_dict_mono['sigma'], rslt_dict_mono['Eej0']/rslt_dict_mono['E10'], 'r--', linewidth=linewidth)
ax.plot(rslt_dict_mono['sigma'], rslt_dict_mono['Eej1']/rslt_dict_mono['E11'], 'g--', linewidth=linewidth)
ax.plot(rslt_dict_mono['sigma'], rslt_dict_mono['Eej2']/rslt_dict_mono['E12'], 'b--', linewidth=linewidth)

ax.set_xlim(0, 1)
#ax.set_ylim(0, 0.6)
ax.set_xlabel('$\\sigma_d$', fontsize=label_size)
ax.set_ylabel('$\\overline{E_{ej}}/E_{i}$', fontsize=label_size)
ax.tick_params(axis='x', labelsize=label_size)
ax.tick_params(axis='y', labelsize=label_size)

inset_ax = inset_axes(ax, width='50%', height='50%', loc='upper left')
inset_ax.plot(sigma_skp, Et_array_0/E1_array_0, 'ro', markersize=marker_size_in, markeredgewidth=marker_width_in)
inset_ax.plot(sigma_skp, Et_array_1/E1_array_1, 'g^', markersize=marker_size_in, markeredgewidth=marker_width_in)
inset_ax.plot(sigma_skp, Et_array_2/E1_array_2, 'bs', markersize=marker_size_in, markeredgewidth=marker_width_in)
inset_ax.plot(rslt_dict_mono['sigma'], Et_array_mono_0/rslt_dict_mono['E10'], 'r--', linewidth=linewidth)
inset_ax.plot(rslt_dict_mono['sigma'], Et_array_mono_1/rslt_dict_mono['E11'], 'g--', linewidth=linewidth)
inset_ax.plot(rslt_dict_mono['sigma'], Et_array_mono_2/rslt_dict_mono['E12'], 'b--', linewidth=linewidth)
inset_ax.set_xticks([0, 0.5, 1])
inset_ax.set_ylabel('$\\overline{E_{ej}} \\cdot  \\overline{N_{ej}}/E_{i}$', fontsize=label_size)

inset_ax.tick_params(axis='x', labelsize=ticks_size)
inset_ax.tick_params(axis='y', labelsize=ticks_size)
# 把y的tick移到右边
inset_ax.yaxis.tick_right()
inset_ax.yaxis.set_label_position("right")

fig = plt.figure(2, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

ax.plot(sigma_skp, N_array_0, 'ro', markersize=marker_size)
ax.plot(sigma_skp, N_array_1, 'g^', markersize=marker_size)
ax.plot(sigma_skp, N_array_2, 'bs', markersize=marker_size)
ax.plot(rslt_dict_mono['sigma'], rslt_dict_mono['Nej0'], 'r--', linewidth=linewidth)
ax.plot(rslt_dict_mono['sigma'], rslt_dict_mono['Nej1'], 'g--', linewidth=linewidth)
ax.plot(rslt_dict_mono['sigma'], rslt_dict_mono['Nej2'], 'b--', linewidth=linewidth)
## y = 1 等值线
#ax.axhline(1, color='k', linestyle='-', linewidth=linewidth)
ax.set_xlim(0, 1)
#ax.set_ylim(0, 10)
ax.set_xlabel('$\\sigma_d$', fontsize=label_size)
ax.set_ylabel('$\\overline{N_{ej}}$', fontsize=label_size)
ax.tick_params(axis='x', labelsize=label_size)
ax.tick_params(axis='y', labelsize=label_size)
# 创建自定义图例句柄
symbol_1 = mlines.Line2D([], [], color='r', marker='o', markersize=marker_size, markeredgewidth=marker_width)
symbol_2 = mlines.Line2D([], [], color='g', marker='^', markersize=marker_size, markeredgewidth=marker_width)
symbol_3 = mlines.Line2D([], [], color='b', marker='s', markersize=marker_size, markeredgewidth=marker_width)
ax.legend([symbol_1, symbol_2, symbol_3],
          ['$v_1/\\sqrt{g \\mathbb{E}[d]} = 30$', '$v_1/\\sqrt{g \\mathbb{E}[d]} = 50$', '$v_1/\\sqrt{g \\mathbb{E}[d]} = 100$'],
          fontsize=label_size,
          loc='upper right',
          bbox_to_anchor=(1.0, 0.95),
          frameon=True)

plt.show()
