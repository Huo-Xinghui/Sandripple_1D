import numpy as np
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerTuple
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from data_process import generate_truncated_lognormal
from get_data import get_exp_data_array, get_model_data_array

# Set up matplotlib parameters
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.usetex'] = True
A = 1/0.5 # amplification factor
label_size = 12.5*A
ticks_size = 10*A
marker_size = 8*A
marker_size_in = 6*A
marker_width = 2
marker_width_in = 2
linewidth = 1.5
color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

# Bed PSD parameters
dist_params = {
    'd_min': 1.5e-4,
    'd_max': 6e-4,
    'mu': -8.30271,
    'sigma': 0.25778,
    'sampling_num': 100
}
# Average bed diameter
sampling_num = 100000
mu = dist_params['mu']
sigma = dist_params['sigma']
d_min = dist_params['d_min']
d_max = dist_params['d_max']
d2_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, sampling_num)
d2_mid = np.percentile(d2_array, 50)
d2_mean = np.mean(d2_array)
# Bed type parameters
bed_type = {
    'three_D': False,  # 3D bed
    'monodisperse': True,  # monodisperse bed
    'd50': d2_mid,  # average bed diameter
    'd_mean': d2_mean
}
# Impactor diameters
d1_dict = {
    'coarse': 1.4*d2_mean,
    'medium': d2_mean,
    'fine': 0.73*d2_mean
}
#d1_coarse = np.percentile(d2_array[(d2_array > 3.55e-4) & (d2_array <= d_max)], 50)
#d1_medium = np.percentile(d2_array[(d2_array > 2.5e-4) & (d2_array <= 3.55e-4)], 50)
#d1_fine = np.percentile(d2_array[(d2_array > d_min) & (d2_array <= 2.5e-4)], 50)
#r_c = d1_coarse/d2_mid
#r_m = d1_medium/d2_mid
#r_f = d1_fine/d2_mid
#d1_dict = {
#    'coarse': d1_coarse,
#    'medium': d1_medium,
#    'fine': d1_fine
#}
# Mechanical properties
epsilon = 0.6442
nu = -0.7857
epsilon_mono = epsilon  # restitution coefficient for monodisperse bed
nu_mono = nu  # friction coefficient for monodisperse bed
epsilon_2D = epsilon  # restitution coefficient for 2D bed
nu_2D = nu  # friction coefficient for 2D bed
epsilon_3D = epsilon  # restitution coefficient for 3D bed
nu_3D = nu  # friction coefficient for 3D bed
# Calculation parameters
perform_calculations = False
thd_min = 0.1  # minimum impact angle
thd_max = 85.0  # maximum impact angle
case_num = 100  # number of cases

# *****************************************************************

start_time = time.time()

# Get experimental data
exp_dicts = get_exp_data_array(d1_dict['coarse'], d1_dict['medium'], d1_dict['fine'])

if perform_calculations:
    current_time = time.time()
    cost_time = current_time - start_time
    print('/n Experimental data loaded (1/10)')
    print(f'Cost time: {cost_time:.2f} s')

    # Get model data
    thd_array = np.linspace(thd_min, thd_max, case_num)  # impact angle array
    th_array = np.radians(thd_array)  # convert to radians
    len_thd = len(th_array)
    d1_array_f = np.array([d1_dict['fine']] * len_thd)
    d1_array_m = np.array([d1_dict['medium']] * len_thd)
    d1_array_c = np.array([d1_dict['coarse']] * len_thd)

    phi_array_f, e_array_f = get_model_data_array(epsilon_mono, nu_mono, d1_array_f, th_array, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n Mono-fine model data calculated (2/10)')
    print(f'Cost time: {cost_time:.2f} s')
    phi_array_m, e_array_m = get_model_data_array(epsilon_mono, nu_mono, d1_array_m, th_array, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n Mono-medium model data calculated (3/10)')
    print(f'Cost time: {cost_time:.2f} s')
    phi_array_c, e_array_c = get_model_data_array(epsilon_mono, nu_mono, d1_array_c, th_array, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n Mono-coarse model data calculated (4/10)')
    print(f'Cost time: {cost_time:.2f} s')

    rslt_dict_mono = {
        'th': th_array,
        'thd': thd_array,
        'phi_f': phi_array_f,
        'phid_f': np.degrees(phi_array_f),
        'e_f': e_array_f,
        'phi_m': phi_array_m,
        'phid_m': np.degrees(phi_array_m),
        'e_m': e_array_m,
        'phi_c': phi_array_c,
        'phid_c': np.degrees(phi_array_c),
        'e_c': e_array_c,
    }
    np.savez('model_results_mono.npz', **rslt_dict_mono)

    bed_type['monodisperse'] = False  # switch to polydisperse bed

    phi_array_f, e_array_f = get_model_data_array(epsilon_2D, nu_2D, d1_array_f, th_array, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n 2D-fine model data calculated (5/10)')
    print(f'Cost time: {cost_time:.2f} s')
    phi_array_m, e_array_m = get_model_data_array(epsilon_2D, nu_2D, d1_array_m, th_array, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n 2D-medium model data calculated (6/10)')
    print(f'Cost time: {cost_time:.2f} s')
    phi_array_c, e_array_c = get_model_data_array(epsilon_2D, nu_2D, d1_array_c, th_array, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n 2D-coarse model data calculated (7/10)')
    print(f'Cost time: {cost_time:.2f} s')

    rslt_dict_2D = {
        'th': th_array,
        'thd': thd_array,
        'phi_f': phi_array_f,
        'phid_f': np.degrees(phi_array_f),
        'e_f': e_array_f,
        'phi_m': phi_array_m,
        'phid_m': np.degrees(phi_array_m),
        'e_m': e_array_m,
        'phi_c': phi_array_c,
        'phid_c': np.degrees(phi_array_c),
        'e_c': e_array_c,
    }
    np.savez('model_results_2D.npz', **rslt_dict_2D)

    bed_type['three_D'] = True  # switch to 3D bed

    phi_array_f, e_array_f = get_model_data_array(epsilon_3D, nu_3D, d1_array_f, th_array, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n 3D-fine model data calculated (8/10)')
    print(f'Cost time: {cost_time:.2f} s')
    phi_array_m, e_array_m = get_model_data_array(epsilon_3D, nu_3D, d1_array_m, th_array, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n 3D-medium model data calculated (9/10)')
    print(f'Cost time: {cost_time:.2f} s')
    phi_array_c, e_array_c = get_model_data_array(epsilon_3D, nu_3D, d1_array_c, th_array, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n 3D-coarse model data calculated (10/10)')
    print(f'Cost time: {cost_time:.2f} s')

    rslt_dict_3D = {
        'th': th_array,
        'thd': thd_array,
        'phi_f': phi_array_f,
        'phid_f': np.degrees(phi_array_f),
        'e_f': e_array_f,
        'phi_m': phi_array_m,
        'phid_m': np.degrees(phi_array_m),
        'e_m': e_array_m,
        'phi_c': phi_array_c,
        'phid_c': np.degrees(phi_array_c),
        'e_c': e_array_c,
    }
    np.savez('model_results_3D.npz', **rslt_dict_3D)

# Load results
rslt_dict_mono = np.load('model_results_mono_same.npz')
rslt_dict_2D = np.load('model_results_2D_same.npz')
rslt_dict_3D = np.load('model_results_3D_old.npz')
rslt_dict_mono_in = np.load('model_results_mono_old.npz')
rslt_dict_2D_in = np.load('model_results_2D_old.npz')
rslt_dict_3D_in = np.load('model_results_3D_old.npz')
rslt_dict_mono_good = np.load('model_results_mono_good.npz')

# Draw rebound angle data
fig = plt.figure(1, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

ax.scatter(
    exp_dicts['Ri95_f']['thd'],
    exp_dicts['Ri95_f']['phid'],
    marker='o',
    edgecolor=color[0],
    facecolor='none',
    s=marker_size**2,
    linewidth=marker_width
)
ax.scatter(
    exp_dicts['Ri95_m']['thd'],
    exp_dicts['Ri95_m']['phid'],
    marker='^',
    edgecolor=color[1],
    facecolor='none',
    s=marker_size**2,
    linewidth=marker_width
)
ax.scatter(
    exp_dicts['Ri95_c']['thd'],
    exp_dicts['Ri95_c']['phid'],
    marker='s',
    edgecolor=color[2],
    facecolor='none',
    s=marker_size**2,
    linewidth=marker_width
)
ax.scatter(
    exp_dicts['Wi89_f']['thd'],
    exp_dicts['Wi89_f']['phid'],
    marker='o',
    color=color[0],
    s=marker_size**2,
    linewidth=marker_width
)
ax.scatter(
    exp_dicts['Wi89_m']['thd'],
    exp_dicts['Wi89_m']['phid'],
    marker='^',
    color=color[1],
    s=marker_size**2,
    linewidth=marker_width
)
ax.scatter(
    exp_dicts['Wi89_c']['thd'],
    exp_dicts['Wi89_c']['phid'],
    marker='s',
    color=color[2],
    s=marker_size**2,
    linewidth=marker_width
)

ax.plot(rslt_dict_mono['thd'], rslt_dict_mono['phid_f'], color[0], linestyle=':', linewidth=linewidth)
ax.plot(rslt_dict_mono['thd'], rslt_dict_mono['phid_m'], color[1], linestyle=':', linewidth=linewidth)
ax.plot(rslt_dict_mono['thd'], rslt_dict_mono['phid_c'], color[2], linestyle=':', linewidth=linewidth)
#ax.plot(rslt_dict_2D['thd'], rslt_dict_2D['phid_f'], 'r--', linewidth=linewidth)
#ax.plot(rslt_dict_2D['thd'], rslt_dict_2D['phid_m'], 'g--', linewidth=linewidth)
#ax.plot(rslt_dict_2D['thd'], rslt_dict_2D['phid_c'], 'b--', linewidth=linewidth)
ax.plot(rslt_dict_3D['thd'], rslt_dict_3D['phid_f'], color[0], linestyle='-', linewidth=linewidth)
ax.plot(rslt_dict_3D['thd'], rslt_dict_3D['phid_m'], color[1], linestyle='-', linewidth=linewidth)
ax.plot(rslt_dict_3D['thd'], rslt_dict_3D['phid_c'], color[2], linestyle='-', linewidth=linewidth)
#ax.plot(rslt_dict_mono_good['thd'], rslt_dict_mono_good['phid_f'], 'r-', linewidth=0.8)
#ax.plot(rslt_dict_mono_good['thd'], rslt_dict_mono_good['phid_m'], 'g-', linewidth=0.8)
#ax.plot(rslt_dict_mono_good['thd'], rslt_dict_mono_good['phid_c'], 'b-', linewidth=0.8)

ax.set_xlim(0, 85)
ax.set_ylim(0, 110)
ax.set_xlabel('$\\theta$ ($^{\\circ}$)', fontsize=label_size)
ax.set_ylabel('$\\overline{\\theta \'}$($^{\\circ}$)', fontsize=label_size)
ax.tick_params(axis='both', labelsize=ticks_size)
# 创建自定义图例句柄
line_mono = mlines.Line2D([], [], color='k', linestyle=':', linewidth=linewidth, label='2D-P')
line_2D = mlines.Line2D([], [], color='k', linestyle='--',linewidth=linewidth, label='2D-M')
line_3D = mlines.Line2D([], [], color='k', linestyle='-',linewidth=linewidth, label='3D-2D')
#line_2DPM = mlines.Line2D([], [], color='k', linestyle='-',linewidth=0.8, label='2D-PM')
#ax.legend([line_2DP, line_2DM, line_3D2D, line_2DPM], ['2D-P', '2D-M', '3D-2D', '2D-PM'], fontsize=ticks_size, loc='best', frameon=True)
ax.legend([line_3D, line_2D, line_mono],
          ['3D', '2D', 'Monodisperse'],
          fontsize=ticks_size,
          loc='upper right',
          bbox_to_anchor=(0.9, 1),
          frameon=True)

inset_ax = inset_axes(ax, width='45%', height='45%', loc='upper left')
inset_ax.scatter(
    exp_dicts['Ri95_f']['thd'],
    exp_dicts['Ri95_f']['phid'],
    marker='o',
    edgecolor=color[0],
    facecolor='none',
    s=marker_size_in**2,
    linewidth=marker_width_in
)
inset_ax.scatter(
    exp_dicts['Ri95_m']['thd'],
    exp_dicts['Ri95_m']['phid'],
    marker='^',
    edgecolor=color[1],
    facecolor='none',
    s=marker_size_in**2,
    linewidth=marker_width_in
)
inset_ax.scatter(
    exp_dicts['Ri95_c']['thd'],
    exp_dicts['Ri95_c']['phid'],
    marker='s',
    edgecolor=color[2],
    facecolor='none',
    s=marker_size_in**2,
    linewidth=marker_width_in
)
inset_ax.scatter(
    exp_dicts['Wi89_f']['thd'],
    exp_dicts['Wi89_f']['phid'],
    marker='o',
    color=color[0],
    s=marker_size_in**2
)
inset_ax.scatter(
    exp_dicts['Wi89_m']['thd'],
    exp_dicts['Wi89_m']['phid'],
    marker='^',
    color=color[1],
    s=marker_size_in**2
)
inset_ax.scatter(
    exp_dicts['Wi89_c']['thd'],
    exp_dicts['Wi89_c']['phid'],
    marker='s',
    color=color[2],
    s=marker_size_in**2
)
inset_ax.plot(rslt_dict_mono_in['thd'], rslt_dict_mono_in['phid_f'], color[0], linestyle=':', linewidth=linewidth)
inset_ax.plot(rslt_dict_mono_in['thd'], rslt_dict_mono_in['phid_m'], color[1], linestyle=':', linewidth=linewidth)
inset_ax.plot(rslt_dict_mono_in['thd'], rslt_dict_mono_in['phid_c'], color[2], linestyle=':', linewidth=linewidth)
inset_ax.plot(rslt_dict_2D_in['thd'], rslt_dict_2D_in['phid_f'], color[0], linestyle='--', linewidth=linewidth)
inset_ax.plot(rslt_dict_2D_in['thd'], rslt_dict_2D_in['phid_m'], color[1], linestyle='--', linewidth=linewidth)
inset_ax.plot(rslt_dict_2D_in['thd'], rslt_dict_2D_in['phid_c'], color[2], linestyle='--', linewidth=linewidth)
inset_ax.plot(rslt_dict_3D_in['thd'], rslt_dict_3D_in['phid_f'], color[0], linestyle='-', linewidth=linewidth)
inset_ax.plot(rslt_dict_3D_in['thd'], rslt_dict_3D_in['phid_m'], color[1], linestyle='-', linewidth=linewidth)
inset_ax.plot(rslt_dict_3D_in['thd'], rslt_dict_3D_in['phid_c'], color[2], linestyle='-', linewidth=linewidth)
inset_ax.set_xlim(0, 42)
inset_ax.set_ylim(0, 62)
inset_ax.tick_params(axis='both', labelsize=ticks_size)
# 把y的tick移到右边
inset_ax.yaxis.tick_right()
inset_ax.yaxis.set_label_position("right")

#plt.tight_layout()


# Draw restitution coefficient data
fig = plt.figure(2, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

ax.plot(
    exp_dicts['Ri95_f']['thd'],
    exp_dicts['Ri95_f']['e'],
    'C0o',
    markerfacecolor='none',
    label='Rice95 fine',
    markersize=marker_size,
    markeredgewidth=marker_width
)
ax.plot(
    exp_dicts['Ri95_m']['thd'],
    exp_dicts['Ri95_m']['e'],
    'C1^',
    markerfacecolor='none',
    label='Rice95 medium',
    markersize=marker_size,
    markeredgewidth=marker_width
)
ax.plot(
    exp_dicts['Ri95_c']['thd'],
    exp_dicts['Ri95_c']['e'],
    'C2s',
    markerfacecolor='none',
    label='Rice95 coarse',
    markersize=marker_size,
    markeredgewidth=marker_width
)
wi89_f = ax.plot(
    exp_dicts['Wi89_f']['thd'],
    exp_dicts['Wi89_f']['e'],
    'C0o',
    label='Willetts89 fine',
    #markerfacecolor='none',
    markersize=marker_size
)
wi89_m = ax.plot(
    exp_dicts['Wi89_m']['thd'],
    exp_dicts['Wi89_m']['e'],
    'C1^',
    label='Willetts89 medium',
    #markerfacecolor='none',
    markersize=marker_size
)
wi89_c = ax.plot(
    exp_dicts['Wi89_c']['thd'],
    exp_dicts['Wi89_c']['e'],
    'C2s',
    label='Willetts89 coarse',
    #markerfacecolor='none',
    markersize=marker_size
)

ax.plot(rslt_dict_mono['thd'], rslt_dict_mono['e_f'], 'C0:', linewidth=linewidth)
ax.plot(rslt_dict_mono['thd'], rslt_dict_mono['e_m'], 'C1:', linewidth=linewidth)
ax.plot(rslt_dict_mono['thd'], rslt_dict_mono['e_c'], 'C2:', linewidth=linewidth)
#ax.plot(rslt_dict_2D['thd'], rslt_dict_2D['e_f'], 'r--', linewidth=linewidth)
#ax.plot(rslt_dict_2D['thd'], rslt_dict_2D['e_m'], 'g--', linewidth=linewidth)
#ax.plot(rslt_dict_2D['thd'], rslt_dict_2D['e_c'], 'b--', linewidth=linewidth)
ax.plot(rslt_dict_3D['thd'], rslt_dict_3D['e_f'], 'C0-', linewidth=linewidth)
ax.plot(rslt_dict_3D['thd'], rslt_dict_3D['e_m'], 'C1-', linewidth=linewidth)
ax.plot(rslt_dict_3D['thd'], rslt_dict_3D['e_c'], 'C2-', linewidth=linewidth)
#ax.plot(rslt_dict_mono_good['thd'], rslt_dict_mono_good['e_f'], 'r-', linewidth=0.8)
#ax.plot(rslt_dict_mono_good['thd'], rslt_dict_mono_good['e_m'], 'g-', linewidth=0.8)
#ax.plot(rslt_dict_mono_good['thd'], rslt_dict_mono_good['e_c'], 'b-', linewidth=0.8)

ax.set_xlim(0, 85)
ax.set_ylim(0.15, 0.8)
ax.set_xlabel('$\\theta$ ($^{\\circ}$)', fontsize=label_size)
ax.set_ylabel('$\\overline{e}$', fontsize=label_size)
ax.tick_params(axis='x', labelsize=ticks_size)
ax.tick_params(axis='y', labelsize=ticks_size)
proxy_ro = plt.Line2D([0], [0],
					 color='C0',
					 marker='o',
					 linestyle='',
					 markersize=marker_size
					)
proxy_ro_h = plt.Line2D([0], [0],
					 color='C0',
					 marker='o',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
proxy_gt = plt.Line2D([0], [0],
					 color='C1',
					 marker='^',
					 linestyle='',
					 markersize=marker_size
					)
proxy_gt_h = plt.Line2D([0], [0],
					 color='C1',
					 marker='^',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
proxy_bs = plt.Line2D([0], [0],
					 color='C2',
					 marker='s',
					 linestyle='',
					 markersize=marker_size
					)
proxy_bs_h = plt.Line2D([0], [0],
					 color='C2',
					 marker='s',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
ax.legend([(proxy_ro, proxy_ro_h),
		   (proxy_gt, proxy_gt_h),
		   (proxy_bs, proxy_bs_h)],
          ['Fine', 'Medium', 'Coarse'],
		  handler_map={tuple: HandlerTuple(ndivide=None)},
          fontsize=ticks_size,
          loc='upper left',
          bbox_to_anchor=(0.16, 1),
          frameon=True,
		  framealpha=1,
		  ncol=1
		  )

inset_ax = inset_axes(ax, width='45%', height='45%', loc='upper right')
inset_ax.plot(
    exp_dicts['Ri95_f']['thd'],
    exp_dicts['Ri95_f']['e'],
    'C0o',
    markerfacecolor='none',
    label='Rice95 fine',
    markersize=marker_size_in,
    markeredgewidth=marker_width_in
)
inset_ax.plot(
    exp_dicts['Ri95_m']['thd'],
    exp_dicts['Ri95_m']['e'],
    'C1^',
    markerfacecolor='none',
    label='Rice95 medium',
    markersize=marker_size_in,
    markeredgewidth=marker_width_in
)
inset_ax.plot(
    exp_dicts['Ri95_c']['thd'],
    exp_dicts['Ri95_c']['e'],
    'C2s',
    markerfacecolor='none',
    label='Rice95 coarse',
    markersize=marker_size_in,
    markeredgewidth=marker_width_in
)
wi89_f = inset_ax.plot(
    exp_dicts['Wi89_f']['thd'],
    exp_dicts['Wi89_f']['e'],
    'C0o',
    label='Willetts89 fine',
    #markerfacecolor='none',
    markersize=marker_size_in
)
wi89_m = inset_ax.plot(
    exp_dicts['Wi89_m']['thd'],
    exp_dicts['Wi89_m']['e'],
    'C1^',
    label='Willetts89 medium',
    #markerfacecolor='none',
    markersize=marker_size_in
)
wi89_c = inset_ax.plot(
    exp_dicts['Wi89_c']['thd'],
    exp_dicts['Wi89_c']['e'],
    'C2s',
    label='Willetts89 coarse',
    #markerfacecolor='none',
    markersize=marker_size_in
)

inset_ax.plot(rslt_dict_mono_in['thd'], rslt_dict_mono_in['e_f'], 'C0:', linewidth=linewidth)
inset_ax.plot(rslt_dict_mono_in['thd'], rslt_dict_mono_in['e_m'], 'C1:', linewidth=linewidth)
inset_ax.plot(rslt_dict_mono_in['thd'], rslt_dict_mono_in['e_c'], 'C2:', linewidth=linewidth)
inset_ax.plot(rslt_dict_2D_in['thd'], rslt_dict_2D_in['e_f'], 'C0--', linewidth=linewidth)
inset_ax.plot(rslt_dict_2D_in['thd'], rslt_dict_2D_in['e_m'], 'C1--', linewidth=linewidth)
inset_ax.plot(rslt_dict_2D_in['thd'], rslt_dict_2D_in['e_c'], 'C2--', linewidth=linewidth)
inset_ax.plot(rslt_dict_3D_in['thd'], rslt_dict_3D_in['e_f'], 'C0-', linewidth=linewidth)
inset_ax.plot(rslt_dict_3D_in['thd'], rslt_dict_3D_in['e_m'], 'C1-', linewidth=linewidth)
inset_ax.plot(rslt_dict_3D_in['thd'], rslt_dict_3D_in['e_c'], 'C2-', linewidth=linewidth)

inset_ax.set_xlim(0, 42)
inset_ax.set_ylim(0.44, 0.64)
inset_ax.tick_params(axis='x', labelsize=ticks_size)
inset_ax.tick_params(axis='y', labelsize=ticks_size)
## 把y的tick移到右边
#inset_ax.yaxis.tick_right()
#inset_ax.yaxis.set_label_position("right")
#inset_ax.xaxis.tick_top()
#inset_ax.xaxis.set_label_position("top")

# Draw vertical restitution coefficient data
plt.figure(3, figsize=(8, 6))
ez_array_f_mono = rslt_dict_mono['e_f'] * np.sin(rslt_dict_mono['phi_f']) / np.sin(rslt_dict_mono['th'])
ez_array_m_mono = rslt_dict_mono['e_m'] * np.sin(rslt_dict_mono['phi_m']) / np.sin(rslt_dict_mono['th'])
ez_array_c_mono = rslt_dict_mono['e_c'] * np.sin(rslt_dict_mono['phi_c']) / np.sin(rslt_dict_mono['th'])
ez_array_f_2D = rslt_dict_2D['e_f'] * np.sin(rslt_dict_2D['phi_f']) / np.sin(rslt_dict_2D['th'])
ez_array_m_2D = rslt_dict_2D['e_m'] * np.sin(rslt_dict_2D['phi_m']) / np.sin(rslt_dict_2D['th'])
ez_array_c_2D = rslt_dict_2D['e_c'] * np.sin(rslt_dict_2D['phi_c']) / np.sin(rslt_dict_2D['th'])
ez_array_f_3D = rslt_dict_3D['e_f'] * np.sin(rslt_dict_3D['phi_f']) / np.sin(rslt_dict_3D['th'])
ez_array_m_3D = rslt_dict_3D['e_m'] * np.sin(rslt_dict_3D['phi_m']) / np.sin(rslt_dict_3D['th'])
ez_array_c_3D = rslt_dict_3D['e_c'] * np.sin(rslt_dict_3D['phi_c']) / np.sin(rslt_dict_3D['th'])

plt.plot(rslt_dict_mono['thd'], ez_array_f_mono, 'r:', label='Model mono fine')
plt.plot(rslt_dict_mono['thd'], ez_array_m_mono, 'g:', label='Model mono medium')
plt.plot(rslt_dict_mono['thd'], ez_array_c_mono, 'b:', label='Model mono coarse')
plt.plot(rslt_dict_2D['thd'], ez_array_f_2D, 'r--', label='Model 2D fine')
plt.plot(rslt_dict_2D['thd'], ez_array_m_2D, 'g--', label='Model 2D medium')
plt.plot(rslt_dict_2D['thd'], ez_array_c_2D, 'b--', label='Model 2D coarse')
plt.plot(rslt_dict_3D['thd'], ez_array_f_3D, 'r-', label='Model 3D fine')
plt.plot(rslt_dict_3D['thd'], ez_array_m_3D, 'g-', label='Model 3D medium')
plt.plot(rslt_dict_3D['thd'], ez_array_c_3D, 'b-', label='Model 3D coarse')

plt.xlim(0, 90)
plt.ylim(0, 2.5)
plt.xlabel('$\\theta$ (degree)', fontsize=label_size)
plt.ylabel('$\\overline{e_z}$', fontsize=label_size)
plt.xticks(fontsize=ticks_size)
plt.yticks(fontsize=ticks_size)
plt.tight_layout()
plt.show()
