import numpy as np
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from data_process import generate_truncated_lognormal, print_time
from get_data import get_exp_data_array, get_model_data_array_v1, get_model_data_array_th

# Set up matplotlib parameters
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.usetex'] = True
A = 1/0.6 # amplification factor
label_size = 20*A
ticks_size = 15*A
marker_size = 10*A
marker_width = 2
# Bed PSD parameters
dist_params = {
    'd_min': 1.5e-4,
    'd_max': 6e-4,
    'mu': -8.30271,
    'sigma': 0.25778,
    'sampling_num': 1000
}
# Average bed diameter
sampling_num = 100000
mu = dist_params['mu']
sigma = dist_params['sigma']
d_min = dist_params['d_min']
d_max = dist_params['d_max']
d2_array = generate_truncated_lognormal(mu, sigma, d_min, d_max, sampling_num)
d2_mid = np.percentile(d2_array, 50)
d90 = np.percentile(d2_array, 90)
# Bed type parameters
bed_type = {
    'three_D': False,  # 3D bed
    'monodisperse': True,  # monodisperse bed
    'd50': d2_mid,  # average bed diameter
    'd90': d90  # 90th percentile diameter
}
# Impactor diameters
d1_dict = {
    'coarse': 1.4*d2_mid,
    'medium': d2_mid,
    'fine': 0.73*d2_mid
}
#d1_coarse = np.mean(d2_array[(d2_array > 3.55e-4) & (d2_array <= d_max)])
#d1_medium = np.mean(d2_array[(d2_array > 2.5e-4) & (d2_array <= 3.55e-4)])
#d1_fine = np.mean(d2_array[(d2_array > d_min) & (d2_array <= 2.5e-4)])
#r_c = d1_coarse/d2_mid
#r_m = d1_medium/d2_mid
#r_f = d1_fine/d2_mid
#d1_dict = {
#    'coarse': d1_coarse,
#    'medium': d1_medium,
#    'fine': d1_fine
#}
# Mechanical properties
rho = 2650
rhof = 1.263
s = rho/rhof
g = 9.8*(1 - 1/s)
epsilon_2DM = 1.0435  # restitution coefficient for monodisperse bed
nu_2DM = -1.2743  # friction coefficient for monodisperse bed
epsilon_2DP = 0.8570  # restitution coefficient for 2D bed
nu_2DP = -1.0726  # friction coefficient for 2D bed
epsilon_3D = 0.8154  # restitution coefficient for 3D bed
nu_3D = -1.0574  # friction coefficient for 3D bed
gamma = 0.047
# Calculation parameters
perform_calculations = False
thd_min = 0.1  # minimum impact angle
thd_max = 85.0  # maximum impact angle
v1_min = 1.0
v1_max = 10.0
case_num = 50  # number of cases

# *****************************************************************

total_steps = 19
current_step = 0
start_time = time.time()
last_time = start_time

# Get experimental data
exp_dicts = get_exp_data_array(d1_dict['coarse'], d1_dict['medium'], d1_dict['fine'])

if perform_calculations:
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    # Get model data vs thd
    thd_array = np.linspace(thd_min, thd_max, case_num)  # impact angle array
    th_array = np.radians(thd_array)  # convert to radians
    len_thd = len(th_array)
    d1_array_f = np.array([d1_dict['fine']] * len_thd)
    d1_array_m = np.array([d1_dict['medium']] * len_thd)
    d1_array_c = np.array([d1_dict['coarse']] * len_thd)
    v1_f = np.mean(exp_dicts['Ri95_f']['v1'])
    v1_m = np.mean(exp_dicts['Ri95_m']['v1'])
    v1_c = np.mean(exp_dicts['Ri95_c']['v1'])

    physical_dict = {
        'epsilon': epsilon_2DM,
        'nu': nu_2DM,
        'rho': rho,
        'g': g
    }

    Nej_array_f, vn_array_f = get_model_data_array_th(th_array, v1_f, gamma, d1_array_f, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    Nej_array_m, vn_array_m = get_model_data_array_th(th_array, v1_m, gamma, d1_array_m, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    Nej_array_c, vn_array_c = get_model_data_array_th(th_array, v1_c, gamma, d1_array_c, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    rslt_dict_2DM_th = {
        'th': th_array,
        'thd': thd_array,
        'Nej_f': Nej_array_f,
        'vn_f': vn_array_f,
        'Nej_m': Nej_array_m,
        'vn_m': vn_array_m,
        'Nej_c': Nej_array_c,
        'vn_c': vn_array_c,
        'v1_f': v1_f,
        'v1_m': v1_m,
        'v1_c': v1_c,
        'd2_mid': d2_mid,
        'd1_f': d1_dict['fine'],
        'd1_m': d1_dict['medium'],
        'd1_c': d1_dict['coarse']
    }
    np.savez('eject_2DM_vs_th.npz', **rslt_dict_2DM_th)
    #np.savez('eject_2DM_vs_th_d50.npz', **rslt_dict_2DM_th)

    bed_type['monodisperse'] = False  # switch to polydisperse bed

    physical_dict = {
        'epsilon': epsilon_2DP,
        'nu': nu_2DP,
        'rho': rho,
        'g': g
    }

    Nej_array_f, vn_array_f = get_model_data_array_th(th_array, v1_f, gamma, d1_array_f, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    Nej_array_m, vn_array_m = get_model_data_array_th(th_array, v1_m, gamma, d1_array_m, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    Nej_array_c, vn_array_c = get_model_data_array_th(th_array, v1_c, gamma, d1_array_c, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    rslt_dict_2DP_th = {
        'th': th_array,
        'thd': thd_array,
        'Nej_f': Nej_array_f,
        'vn_f': vn_array_f,
        'Nej_m': Nej_array_m,
        'vn_m': vn_array_m,
        'Nej_c': Nej_array_c,
        'vn_c': vn_array_c,
        'v1_f': v1_f,
        'v1_m': v1_m,
        'v1_c': v1_c,
        'd2_mid': d2_mid,
        'd1_f': d1_dict['fine'],
        'd1_m': d1_dict['medium'],
        'd1_c': d1_dict['coarse']
    }
    np.savez('eject_2DP_vs_th.npz', **rslt_dict_2DP_th)
    #np.savez('eject_2DP_vs_th_d50.npz', **rslt_dict_2DP_th)

    bed_type['three_D'] = True  # switch to 3D bed

    physical_dict = {
        'epsilon': epsilon_3D,
        'nu': nu_3D,
        'rho': rho,
        'g': g
    }

    Nej_array_f, vn_array_f = get_model_data_array_th(th_array, v1_f, gamma, d1_array_f, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    Nej_array_m, vn_array_m = get_model_data_array_th(th_array, v1_m, gamma, d1_array_m, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    Nej_array_c, vn_array_c = get_model_data_array_th(th_array, v1_c, gamma, d1_array_c, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    rslt_dict_3D_th = {
        'th': th_array,
        'thd': thd_array,
        'Nej_f': Nej_array_f,
        'vn_f': vn_array_f,
        'Nej_m': Nej_array_m,
        'vn_m': vn_array_m,
        'Nej_c': Nej_array_c,
        'vn_c': vn_array_c,
        'v1_f': v1_f,
        'v1_m': v1_m,
        'v1_c': v1_c,
        'd2_mid': d2_mid,
        'd1_f': d1_dict['fine'],
        'd1_m': d1_dict['medium'],
        'd1_c': d1_dict['coarse']
    }
    np.savez('eject_3D_vs_th.npz', **rslt_dict_3D_th)
    #np.savez('eject_3D_vs_th_d50.npz', **rslt_dict_3D_th)

    # Get model data vs v1
    v1_array = np.linspace(v1_min, v1_max, case_num)  # impact angle array
    len_v1 = len(v1_array)
    d1_array_f = np.array([d1_dict['fine']] * len_v1)
    d1_array_m = np.array([d1_dict['medium']] * len_v1)
    d1_array_c = np.array([d1_dict['coarse']] * len_v1)
    thd_f = np.mean(exp_dicts['Ri95_f']['thd'])
    thd_m = np.mean(exp_dicts['Ri95_m']['thd'])
    thd_c = np.mean(exp_dicts['Ri95_c']['thd'])
    th_f = np.radians(thd_f)
    th_m = np.radians(thd_m)
    th_c = np.radians(thd_c)

    bed_type['monodisperse'] = True  # switch to monodisperse bed
    bed_type['three_D'] = False  # switch to 2D bed

    physical_dict = {
        'epsilon': epsilon_2DM,
        'nu': nu_2DM,
        'rho': rho,
        'g': g
    }

    Nej_array_f, vn_array_f = get_model_data_array_v1(th_f, v1_array, gamma, d1_array_f, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    Nej_array_m, vn_array_m = get_model_data_array_v1(th_m, v1_array, gamma, d1_array_m, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    Nej_array_c, vn_array_c = get_model_data_array_v1(th_c, v1_array, gamma, d1_array_c, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    rslt_dict_2DM_v1 = {
        'v1': v1_array,
        'Nej_f': Nej_array_f,
        'vn_f': vn_array_f,
        'Nej_m': Nej_array_m,
        'vn_m': vn_array_m,
        'Nej_c': Nej_array_c,
        'vn_c': vn_array_c,
        'thd_f': thd_f,
        'thd_m': thd_m,
        'thd_c': thd_c,
        'd2_mid': d2_mid,
        'd1_f': d1_dict['fine'],
        'd1_m': d1_dict['medium'],
        'd1_c': d1_dict['coarse']
    }
    np.savez('eject_2DM_vs_v1.npz', **rslt_dict_2DM_v1)
    #np.savez('eject_2DM_vs_v1_d50.npz', **rslt_dict_2DM_v1)

    bed_type['monodisperse'] = False  # switch to polydisperse bed

    physical_dict = {
        'epsilon': epsilon_2DP,
        'nu': nu_2DP,
        'rho': rho,
        'g': g
    }

    Nej_array_f, vn_array_f = get_model_data_array_v1(th_f, v1_array, gamma, d1_array_f, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    Nej_array_m, vn_array_m = get_model_data_array_v1(th_m, v1_array, gamma, d1_array_m, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    Nej_array_c, vn_array_c = get_model_data_array_v1(th_c, v1_array, gamma, d1_array_c, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    rslt_dict_2DP_v1 = {
        'v1': v1_array,
        'Nej_f': Nej_array_f,
        'vn_f': vn_array_f,
        'Nej_m': Nej_array_m,
        'vn_m': vn_array_m,
        'Nej_c': Nej_array_c,
        'vn_c': vn_array_c,
        'thd_f': thd_f,
        'thd_m': thd_m,
        'thd_c': thd_c,
        'd2_mid': d2_mid,
        'd1_f': d1_dict['fine'],
        'd1_m': d1_dict['medium'],
        'd1_c': d1_dict['coarse']
    }
    np.savez('eject_2DP_vs_v1.npz', **rslt_dict_2DP_v1)
    #np.savez('eject_2DP_vs_v1_d50.npz', **rslt_dict_2DP_v1)

    bed_type['three_D'] = True  # switch to 3D bed

    physical_dict = {
        'epsilon': epsilon_3D,
        'nu': nu_3D,
        'rho': rho,
        'g': g
    }

    Nej_array_f, vn_array_f = get_model_data_array_v1(th_f, v1_array, gamma, d1_array_f, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    Nej_array_m, vn_array_m = get_model_data_array_v1(th_m, v1_array, gamma, d1_array_m, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    Nej_array_c, vn_array_c = get_model_data_array_v1(th_c, v1_array, gamma, d1_array_c, physical_dict, bed_type, dist_params)
    current_time = time.time()
    current_step += 1
    print_time(current_time, last_time, start_time, current_step, total_steps)
    last_time = current_time

    rslt_dict_3D_v1 = {
        'v1': v1_array,
        'Nej_f': Nej_array_f,
        'vn_f': vn_array_f,
        'Nej_m': Nej_array_m,
        'vn_m': vn_array_m,
        'Nej_c': Nej_array_c,
        'vn_c': vn_array_c,
        'thd_f': thd_f,
        'thd_m': thd_m,
        'thd_c': thd_c,
        'd2_mid': d2_mid,
        'd1_f': d1_dict['fine'],
        'd1_m': d1_dict['medium'],
        'd1_c': d1_dict['coarse']
    }
    np.savez('eject_3D_vs_v1.npz', **rslt_dict_3D_v1)
    #np.savez('eject_3D_vs_v1_d50.npz', **rslt_dict_3D_v1)

# Load results
rslt_2DM_vs_th = np.load('eject_2DM_vs_th.npz')
rslt_2DP_vs_th = np.load('eject_2DP_vs_th.npz')
rslt_3D_vs_th = np.load('eject_3D_vs_th.npz')
rslt_2DM_vs_v1 = np.load('eject_2DM_vs_v1.npz')
rslt_2DP_vs_v1 = np.load('eject_2DP_vs_v1.npz')
rslt_3D_vs_v1 = np.load('eject_3D_vs_v1.npz')

rslt_2DM_vs_th_d50 = np.load('eject_2DM_vs_th_d50.npz')
rslt_2DP_vs_th_d50 = np.load('eject_2DP_vs_th_d50.npz')
rslt_3D_vs_th_d50 = np.load('eject_3D_vs_th_d50.npz')
rslt_2DM_vs_v1_d50 = np.load('eject_2DM_vs_v1_d50.npz')
rslt_2DP_vs_v1_d50 = np.load('eject_2DP_vs_v1_d50.npz')
rslt_3D_vs_v1_d50 = np.load('eject_3D_vs_v1_d50.npz')

# Draw Nej vs thd
plt.figure(1, figsize=(8, 6))
p1=plt.plot(
    exp_dicts['Ri95_f']['thd'],
    exp_dicts['Ri95_f']['Nej'],
    'ro',
    label='Rice95 fine',
    markersize=marker_size
)
p2=plt.plot(
    exp_dicts['Ri95_m']['thd'],
    exp_dicts['Ri95_m']['Nej'],
    'g^',
    label='Rice95 medium',
    markersize=marker_size
)

p3=plt.plot(
    exp_dicts['Ri95_c']['thd'],
    exp_dicts['Ri95_c']['Nej'],
    'bs',
    label='Rice95 coarse',
    markersize=marker_size
)

plt.plot(
    exp_dicts['Yi21_f_N_t']['thd'],
    exp_dicts['Yi21_f_N_t']['Nej'],
    'ro',
    label='Yin21 fine',
    markerfacecolor='none',
    markeredgewidth=marker_width,
    markersize=marker_size
)

plt.plot(
    exp_dicts['Yi21_m_N_t']['thd'],
    exp_dicts['Yi21_m_N_t']['Nej'],
    'g^',
    label='Yin21 medium',
    markerfacecolor='none',
    markeredgewidth=marker_width,
    markersize=marker_size
)

plt.plot(
    exp_dicts['Yi21_c_N_t']['thd'],
    exp_dicts['Yi21_c_N_t']['Nej'],
    'bs',
    label='Yin21 coarse',
    markerfacecolor='none',
    markeredgewidth=marker_width,
    markersize=marker_size
)

plt.plot(rslt_2DM_vs_th['thd'], rslt_2DM_vs_th['Nej_f'], 'r:', label='Model mono fine')
plt.plot(rslt_2DM_vs_th['thd'], rslt_2DM_vs_th['Nej_m'], 'g:', label='Model mono medium')
plt.plot(rslt_2DM_vs_th['thd'], rslt_2DM_vs_th['Nej_c'], 'b:', label='Model mono coarse')
plt.plot(rslt_3D_vs_th_d50['thd'], rslt_3D_vs_th_d50['Nej_f'], 'r--', label='Model 3D fine d50')
plt.plot(rslt_3D_vs_th_d50['thd'], rslt_3D_vs_th_d50['Nej_m'], 'g--', label='Model 3D medium d50')
plt.plot(rslt_3D_vs_th_d50['thd'], rslt_3D_vs_th_d50['Nej_c'], 'b--', label='Model 3D coarse d50')
plt.plot(rslt_3D_vs_th['thd'], rslt_3D_vs_th['Nej_f'], 'r-', label='Model 3D fine d90')
plt.plot(rslt_3D_vs_th['thd'], rslt_3D_vs_th['Nej_m'], 'g-', label='Model 3D medium d90')
plt.plot(rslt_3D_vs_th['thd'], rslt_3D_vs_th['Nej_c'], 'b-', label='Model 3D coarse d90')

plt.xlim(0, 80)
plt.ylim(0, 12)
plt.xlabel('$\\theta$ (degree)', fontsize=label_size)
plt.ylabel('$N_{e}$', fontsize=label_size)
plt.xticks(fontsize=ticks_size)
plt.yticks(fontsize=ticks_size)
plt.legend([p1[0], p2[0], p3[0]],
           ['Fine', 'Medium', 'Coarse'],
           fontsize=ticks_size,
           loc='upper right',
           frameon=True)
plt.tight_layout()

# Draw Nej vs v1
plt.figure(2, figsize=(8, 6))
sgd = np.sqrt(g*d2_mid)
sgd1 = np.sqrt(g*3.2e-4)
p1=plt.plot(
    np.array(exp_dicts['Ri95_f']['v1'])/sgd,
    np.array(exp_dicts['Ri95_f']['Nej']),
    'ro',
    label='Rice95 fine',
    markersize=marker_size
)
p2=plt.plot(
    np.array(exp_dicts['Ri95_m']['v1'])/sgd,
    np.array(exp_dicts['Ri95_m']['Nej']),
    'g^',
    label='Rice95 medium',
    markersize=marker_size
)

p3=plt.plot(
    np.array(exp_dicts['Ri95_c']['v1'])/sgd,
    np.array(exp_dicts['Ri95_c']['Nej']),
    'bs',
    label='Rice95 coarse',
    markersize=marker_size
)

plt.plot(
    np.array(exp_dicts['Yi21_f_N_v']['v1'])/sgd1,
    np.array(exp_dicts['Yi21_f_N_v']['Nej']),
    'ro',
    label='Yin21 fine',
    markerfacecolor='none',
    markeredgewidth=marker_width,
    markersize=marker_size
)

plt.plot(
    np.array(exp_dicts['Yi21_m_N_v']['v1'])/sgd1,
    np.array(exp_dicts['Yi21_m_N_v']['Nej']),
    'g^',
    label='Yin21 medium',
    markerfacecolor='none',
    markeredgewidth=marker_width,
    markersize=marker_size
)

plt.plot(
    np.array(exp_dicts['Yi21_c_N_v']['v1'])/sgd1,
    np.array(exp_dicts['Yi21_c_N_v']['Nej']),
    'bs',
    label='Yin21 coarse',
    markerfacecolor='none',
    markeredgewidth=marker_width,
    markersize=marker_size
)


plt.plot(rslt_2DM_vs_v1['v1']/sgd, rslt_2DM_vs_v1['Nej_f'], 'r:', label='Model mono fine')
plt.plot(rslt_2DM_vs_v1['v1']/sgd, rslt_2DM_vs_v1['Nej_m'], 'g:', label='Model mono medium')
plt.plot(rslt_2DM_vs_v1['v1']/sgd, rslt_2DM_vs_v1['Nej_c'], 'b:', label='Model mono coarse')
plt.plot(rslt_3D_vs_v1_d50['v1']/sgd, rslt_3D_vs_v1_d50['Nej_f'], 'r--', label='Model 3D fine d50')
plt.plot(rslt_3D_vs_v1_d50['v1']/sgd, rslt_3D_vs_v1_d50['Nej_m'], 'g--', label='Model 3D medium d50')
plt.plot(rslt_3D_vs_v1_d50['v1']/sgd, rslt_3D_vs_v1_d50['Nej_c'], 'b--', label='Model 3D coarse d50')
plt.plot(rslt_3D_vs_v1['v1']/sgd, rslt_3D_vs_v1['Nej_f'], 'r-', label='Model 3D fine d90')
plt.plot(rslt_3D_vs_v1['v1']/sgd, rslt_3D_vs_v1['Nej_m'], 'g-', label='Model 3D medium d90')
plt.plot(rslt_3D_vs_v1['v1']/sgd, rslt_3D_vs_v1['Nej_c'], 'b-', label='Model 3D coarse d90')

plt.xlim(20, 150)
plt.ylim(0, 30)
plt.xlabel('$v_1/\\sqrt{\\hat{g} d_{50}}$', fontsize=label_size)
plt.ylabel('$N_{e}$', fontsize=label_size)
plt.xticks(fontsize=ticks_size)
plt.yticks(fontsize=ticks_size)
# 创建自定义图例句柄
line_2DP = mlines.Line2D([], [], color='k', linestyle=':', label='2D-P')
line_2DM = mlines.Line2D([], [], color='k', linestyle='--', label='2D-M')
line_3D2D = mlines.Line2D([], [], color='k', linestyle='-', label='3D-2D')
plt.legend([line_2DP, line_2DM, line_3D2D], ['2D-P', '3D-2D $z_c=d_{50}$', '3D-2D $z_c=d_{90}$'], fontsize=ticks_size, loc='best', frameon=True)
plt.tight_layout()

# Draw vn vs thd
plt.figure(3, figsize=(8, 6))
sgd = np.sqrt(g*d2_mid)
p1=plt.plot(
    exp_dicts['Ri95_f']['thd'],
    exp_dicts['Ri95_f']['vn']/sgd,
    'ro',
    label='Rice95 fine',
    markersize=marker_size
)
p2=plt.plot(
    exp_dicts['Ri95_m']['thd'],
    exp_dicts['Ri95_m']['vn']/sgd,
    'g^',
    label='Rice95 medium',
    markersize=marker_size
)

p3=plt.plot(
    exp_dicts['Ri95_c']['thd'],
    exp_dicts['Ri95_c']['vn']/sgd,
    'bs',
    label='Rice95 coarse',
    markersize=marker_size
)

#plt.plot(
#    exp_dicts['Wi89_f']['thd'],
#    exp_dicts['Wi89_f']['vn']/sgd,
#    'ro',
#    label='Yin21 fine',
#    markerfacecolor='none',
#    markersize=marker_size
#)
#
#plt.plot(
#    exp_dicts['Wi89_m']['thd'],
#    exp_dicts['Wi89_m']['vn']/sgd,
#    'g^',
#    label='Yin21 medium',
#    markerfacecolor='none',
#    markersize=marker_size
#)
#
#plt.plot(
#    exp_dicts['Wi89_c']['thd'],
#    exp_dicts['Wi89_c']['vn']/sgd,
#    'bs',
#    label='Yin21 coarse',
#    markerfacecolor='none',
#    markersize=marker_size
#)

plt.plot(rslt_2DM_vs_th['thd'], rslt_2DM_vs_th['vn_f']/sgd, 'r:', label='Model mono fine')
plt.plot(rslt_2DM_vs_th['thd'], rslt_2DM_vs_th['vn_m']/sgd, 'g:', label='Model mono medium')
plt.plot(rslt_2DM_vs_th['thd'], rslt_2DM_vs_th['vn_c']/sgd, 'b:', label='Model mono coarse')
plt.plot(rslt_3D_vs_th_d50['thd'], rslt_3D_vs_th_d50['vn_f']/sgd, 'r--', label='Model 3D fine d50')
plt.plot(rslt_3D_vs_th_d50['thd'], rslt_3D_vs_th_d50['vn_m']/sgd, 'g--', label='Model 3D medium d50')
plt.plot(rslt_3D_vs_th_d50['thd'], rslt_3D_vs_th_d50['vn_c']/sgd, 'b--', label='Model 3D coarse d50')
plt.plot(rslt_3D_vs_th['thd'], rslt_3D_vs_th['vn_f']/sgd, 'r-', label='Model 3D fine D90')
plt.plot(rslt_3D_vs_th['thd'], rslt_3D_vs_th['vn_m']/sgd, 'g-', label='Model 3D medium D90')
plt.plot(rslt_3D_vs_th['thd'], rslt_3D_vs_th['vn_c']/sgd, 'b-', label='Model 3D coarse D90')

plt.xlim(0, 80)
plt.ylim(3, 6.5)
plt.xlabel('$\\theta$ (degree)', fontsize=label_size)
plt.ylabel('$\\overline{v_{e}}/\\sqrt{\\hat{g} d_{50}}$', fontsize=label_size)
plt.xticks(fontsize=ticks_size)
plt.yticks(fontsize=ticks_size)
#plt.legend([p1[0], p2[0], p3[0]],
#           ['Fine', 'Medium', 'Coarse'],
#           fontsize=ticks_size,
#           loc='upper right',
#           frameon=True)
plt.tight_layout()

# Draw vn vs v1
plt.figure(4, figsize=(8, 6))
sgd = np.sqrt(g*d2_mid)
sgd1 = np.sqrt(g*3.2e-4)
p1=plt.plot(
    np.array(exp_dicts['Ri95_f']['v1'])/sgd,
    np.array(exp_dicts['Ri95_f']['vn'])/sgd,
    'ro',
    label='Rice95 fine',
    markersize=marker_size
)
p2=plt.plot(
    np.array(exp_dicts['Ri95_m']['v1'])/sgd,
    np.array(exp_dicts['Ri95_m']['vn'])/sgd,
    'g^',
    label='Rice95 medium',
    markersize=marker_size
)

p3=plt.plot(
    np.array(exp_dicts['Ri95_c']['v1'])/sgd,
    np.array(exp_dicts['Ri95_c']['vn'])/sgd,
    'bs',
    label='Rice95 coarse',
    markersize=marker_size
)

plt.plot(
    np.array(exp_dicts['Yi21_vn_v1']['v1'])/sgd1,
    np.array(exp_dicts['Yi21_vn_v1']['vn'])/sgd1,
    'go',
    label='Yin21',
    markerfacecolor='none',
    markeredgewidth=marker_width,
    markersize=marker_size
)

#plt.plot(
#    exp_dicts['Wi89_f']['v1']/sgd,
#    exp_dicts['Wi89_f']['vn']/sgd,
#    'ro',
#    label='Yin21 fine',
#    markerfacecolor='none',
#    markersize=marker_size
#)
#
#plt.plot(
#    exp_dicts['Wi89_m']['v1']/sgd,
#    exp_dicts['Wi89_m']['vn']/sgd,
#    'g^',
#    label='Yin21 medium',
#    markerfacecolor='none',
#    markersize=marker_size
#)
#
#plt.plot(
#    exp_dicts['Wi89_c']['v1']/sgd,
#    exp_dicts['Wi89_c']['vn']/sgd,
#    'bs',
#    label='Yin21 coarse',
#    markerfacecolor='none',
#    markersize=marker_size
#)

plt.plot(rslt_2DM_vs_v1['v1']/sgd, rslt_2DM_vs_v1['vn_f']/sgd, 'r:', label='Model mono fine')
plt.plot(rslt_2DM_vs_v1['v1']/sgd, rslt_2DM_vs_v1['vn_m']/sgd, 'g:', label='Model mono medium')
plt.plot(rslt_2DM_vs_v1['v1']/sgd, rslt_2DM_vs_v1['vn_c']/sgd, 'b:', label='Model mono coarse')
plt.plot(rslt_3D_vs_v1_d50['v1']/sgd, rslt_3D_vs_v1_d50['vn_f']/sgd, 'r--', label='Model 2D fine d50')
plt.plot(rslt_3D_vs_v1_d50['v1']/sgd, rslt_3D_vs_v1_d50['vn_m']/sgd, 'g--', label='Model 2D medium d50')
plt.plot(rslt_3D_vs_v1_d50['v1']/sgd, rslt_3D_vs_v1_d50['vn_c']/sgd, 'b--', label='Model 2D coarse d50')
plt.plot(rslt_3D_vs_v1['v1']/sgd, rslt_3D_vs_v1['vn_f']/sgd, 'r-', label='Model 3D fine d90')
plt.plot(rslt_3D_vs_v1['v1']/sgd, rslt_3D_vs_v1['vn_m']/sgd, 'g-', label='Model 3D medium d90')
plt.plot(rslt_3D_vs_v1['v1']/sgd, rslt_3D_vs_v1['vn_c']/sgd, 'b-', label='Model 3D coarse d90')

plt.xlim(20, 150)
plt.ylim(3, 6.5)
plt.xlabel('$v_1/\\sqrt{\\hat{g} d_{50}}$', fontsize=label_size)
plt.ylabel('$\\overline{v_{e}}/\\sqrt{\\hat{g} d_{50}}$', fontsize=label_size)
plt.xticks(fontsize=ticks_size)
plt.yticks(fontsize=ticks_size)
## 创建自定义图例句柄
#line_2DP = mlines.Line2D([], [], color='k', linestyle=':', label='2D-P')
#line_2DM = mlines.Line2D([], [], color='k', linestyle='--', label='2D-M')
#line_3D2D = mlines.Line2D([], [], color='k', linestyle='-', label='3D-2D')
#plt.legend([line_2DP, line_2DM, line_3D2D], ['2D-P', '3D-2D $z_c=d_{50}$', '3D-2D $z_c=d_{90}$'], fontsize=ticks_size, loc='best', frameon=True)
plt.tight_layout()

plt.show()