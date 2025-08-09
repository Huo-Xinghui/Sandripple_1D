import numpy as np
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
from data_process import generate_truncated_lognormal
from get_data import get_exp_data_array, get_model_data_array

# Set up matplotlib parameters
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.usetex'] = True
label_size = 20
ticks_size = 15
marker_size = 10
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
# Bed type parameters
bed_type = {
    'three_D': False,  # 3D bed
    'monodisperse': True,  # monodisperse bed
    'd50': d2_mid,  # average bed diameter
}
# Impactor diameters
d1_coarse = np.mean(d2_array[(d2_array > 3.55e-4) & (d2_array <= d_max)])
d1_medium = np.mean(d2_array[(d2_array > 2.5e-4) & (d2_array <= 3.55e-4)])
d1_fine = np.mean(d2_array[(d2_array > d_min) & (d2_array <= 2.5e-4)])
#d1_dict = {
#    'coarse': 1.4*d2_mid,
#    'medium': d2_mid,
#    'fine': 0.73*d2_mid
#}
r_c = d1_coarse/d2_mid
r_m = d1_medium/d2_mid
r_f = d1_fine/d2_mid
d1_dict = {
    'coarse': d1_coarse,
    'medium': d1_medium,
    'fine': d1_fine
}
# Mechanical properties
epsilon_mono = 1.0435  # restitution coefficient for monodisperse bed
nu_mono = -1.2743  # friction coefficient for monodisperse bed
epsilon_2D = 0.8570  # restitution coefficient for 2D bed
nu_2D = -1.0726  # friction coefficient for 2D bed
epsilon_3D = 0.8154  # restitution coefficient for 3D bed
nu_3D = -1.0574  # friction coefficient for 3D bed
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
else:
    # Load results
    rslt_dict_mono = np.load('model_results_mono.npz')
    rslt_dict_2D = np.load('model_results_2D.npz')
    rslt_dict_3D = np.load('model_results_3D.npz')

# Draw rebound angle data
plt.figure(1, figsize=(8, 6))
plt.errorbar(
    exp_dicts['Ri95_f']['thd'],
    exp_dicts['Ri95_f']['phid'], 
    xerr=exp_dicts['Ri95_f']['thd_std'], 
    yerr=exp_dicts['Ri95_f']['phid_std'],
    fmt='ro',
    markerfacecolor='none',
    label='Rice95 fine',
    markersize=marker_size,
    capsize=5
)
plt.errorbar(
    exp_dicts['Ri95_m']['thd'],
    exp_dicts['Ri95_m']['phid'], 
    xerr=exp_dicts['Ri95_m']['thd_std'], 
    yerr=exp_dicts['Ri95_m']['phid_std'],
    fmt='g^',
    markerfacecolor='none',
    label='Rice95 medium',
    markersize=marker_size,
    capsize=5
)
plt.errorbar(
    exp_dicts['Ri95_c']['thd'],
    exp_dicts['Ri95_c']['phid'], 
    xerr=exp_dicts['Ri95_c']['thd_std'], 
    yerr=exp_dicts['Ri95_c']['phid_std'],
    fmt='bs',
    markerfacecolor='none',
    label='Rice95 coarse',
    markersize=marker_size,
    capsize=5
)
wi89_f = plt.plot(
    exp_dicts['Wi89_f']['thd'],
    exp_dicts['Wi89_f']['phid'],
    'ro',
    label='Willetts89 fine',
    markersize=marker_size
)
wi89_m = plt.plot(
    exp_dicts['Wi89_m']['thd'],
    exp_dicts['Wi89_m']['phid'],
    'g^',
    label='Willetts89 medium',
    markersize=marker_size
)
wi89_c = plt.plot(
    exp_dicts['Wi89_c']['thd'],
    exp_dicts['Wi89_c']['phid'],
    'bs',
    label='Willetts89 coarse',
    markersize=marker_size
)

plt.plot(rslt_dict_mono['thd'], rslt_dict_mono['phid_f'], 'r:', label='Model mono fine')
plt.plot(rslt_dict_mono['thd'], rslt_dict_mono['phid_m'], 'g:', label='Model mono medium')
plt.plot(rslt_dict_mono['thd'], rslt_dict_mono['phid_c'], 'b:', label='Model mono coarse')
plt.plot(rslt_dict_2D['thd'], rslt_dict_2D['phid_f'], 'r--', label='Model 2D fine')
plt.plot(rslt_dict_2D['thd'], rslt_dict_2D['phid_m'], 'g--', label='Model 2D medium')
plt.plot(rslt_dict_2D['thd'], rslt_dict_2D['phid_c'], 'b--', label='Model 2D coarse')
plt.plot(rslt_dict_3D['thd'], rslt_dict_3D['phid_f'], 'r-', label='Model 3D fine')
plt.plot(rslt_dict_3D['thd'], rslt_dict_3D['phid_m'], 'g-', label='Model 3D medium')
plt.plot(rslt_dict_3D['thd'], rslt_dict_3D['phid_c'], 'b-', label='Model 3D coarse')

plt.xlim(0, thd_max)
#plt.ylim(0, 60)
plt.xlabel('$\\theta$ (degree)', fontsize=label_size)
plt.ylabel('$\\overline{\\theta \'}$(degree)', fontsize=label_size)
plt.xticks(fontsize=ticks_size)
plt.yticks(fontsize=ticks_size)
plt.legend([wi89_f[0], wi89_m[0], wi89_c[0]], ['Fine', 'Medium', 'Coarse'], fontsize=ticks_size, loc='best')

# Draw restitution coefficient data
plt.figure(2, figsize=(8, 6))
plt.errorbar(
    exp_dicts['Ri95_f']['thd'],
    exp_dicts['Ri95_f']['e'], 
    xerr=exp_dicts['Ri95_f']['thd_std'], 
    yerr=exp_dicts['Ri95_f']['e_std'],
    fmt='ro',
    markerfacecolor='none',
    label='Rice95 fine',
    markersize=marker_size,
    capsize=5
)
plt.errorbar(
    exp_dicts['Ri95_m']['thd'],
    exp_dicts['Ri95_m']['e'], 
    xerr=exp_dicts['Ri95_m']['thd_std'], 
    yerr=exp_dicts['Ri95_m']['e_std'],
    fmt='g^',
    markerfacecolor='none',
    label='Rice95 medium',
    markersize=marker_size,
    capsize=5
)
plt.errorbar(
    exp_dicts['Ri95_c']['thd'],
    exp_dicts['Ri95_c']['e'], 
    xerr=exp_dicts['Ri95_c']['thd_std'], 
    yerr=exp_dicts['Ri95_c']['e_std'],
    fmt='bs',
    markerfacecolor='none',
    label='Rice95 coarse',
    markersize=marker_size,
    capsize=5
)
plt.plot(
    exp_dicts['Wi89_f']['thd'],
    exp_dicts['Wi89_f']['e'],
    'ro',
    label='Willetts89 fine',
    markersize=marker_size
)
plt.plot(
    exp_dicts['Wi89_m']['thd'],
    exp_dicts['Wi89_m']['e'],
    'g^',
    label='Willetts89 medium',
    markersize=marker_size
)
plt.plot(
    exp_dicts['Wi89_c']['thd'],
    exp_dicts['Wi89_c']['e'],
    'bs',
    label='Willetts89 coarse',
    markersize=marker_size
)

plt.plot(rslt_dict_mono['thd'], rslt_dict_mono['e_f'], 'r:', label='Model mono fine')
plt.plot(rslt_dict_mono['thd'], rslt_dict_mono['e_m'], 'g:', label='Model mono medium')
plt.plot(rslt_dict_mono['thd'], rslt_dict_mono['e_c'], 'b:', label='Model mono coarse')
plt.plot(rslt_dict_2D['thd'], rslt_dict_2D['e_f'], 'r--', label='Model 2D fine')
plt.plot(rslt_dict_2D['thd'], rslt_dict_2D['e_m'], 'g--', label='Model 2D medium')
plt.plot(rslt_dict_2D['thd'], rslt_dict_2D['e_c'], 'b--', label='Model 2D coarse')
plt.plot(rslt_dict_3D['thd'], rslt_dict_3D['e_f'], 'r-', label='Model 3D fine')
plt.plot(rslt_dict_3D['thd'], rslt_dict_3D['e_m'], 'g-', label='Model 3D medium')
plt.plot(rslt_dict_3D['thd'], rslt_dict_3D['e_c'], 'b-', label='Model 3D coarse')

plt.xlim(0, thd_max)
#plt.ylim(0.4, 0.65)
plt.xlabel('$\\theta$ (degree)', fontsize=label_size)
plt.ylabel('$\\overline{e}$', fontsize=label_size)
plt.xticks(fontsize=ticks_size)
plt.yticks(fontsize=ticks_size)
plt.legend([wi89_f[0], wi89_m[0], wi89_c[0]], ['Fine', 'Medium', 'Coarse'], fontsize=ticks_size, loc='best')

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

plt.show()
