import numpy as np
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
from data_process import generate_truncated_lognormal
from get_data import get_model_data_array

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
sampling_num = 10000
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
# Mechanical properties
epsilon_mono = 0.8229  # restitution coefficient for monodisperse bed
nu_mono = -1.0049  # friction coefficient for monodisperse bed
epsilon_2D = 0.6814  # restitution coefficient for 2D bed
nu_2D = -0.88  # friction coefficient for 2D bed
epsilon_3D = 0.6405  # restitution coefficient for 3D bed
nu_3D = -0.8078  # friction coefficient for 3D bed
# Calculation parameters
perform_calculations = True
thd_small = 5.0  # small impact angle
thd_medium = 45.0  # medium impact angle
thd_large = 85.0  # large impact angle
d1_min = 0.1*d2_mid # minimum impactor diameter
d1_max = 10*d2_mid # maximum impactor diameter
case_num = 100  # number of cases

# *****************************************************************

start_time = time.time()

if perform_calculations:
    # Get model data
    #d1_array = np.linspace(d1_min, d1_max, case_num)
    # d1_array 以指数方式分布
    d1_array = np.logspace(np.log10(d1_min), np.log10(d1_max), case_num)
    th_small = np.radians(thd_small)
    th_medium = np.radians(thd_medium)
    th_large = np.radians(thd_large)
    len_d1 = len(d1_array)
    th_array_s = np.array([th_small] * len_d1)
    th_array_m = np.array([th_medium] * len_d1)
    th_array_l = np.array([th_large] * len_d1)

    phi_array_s, e_array_s = get_model_data_array(epsilon_mono, nu_mono, d1_array, th_array_s, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n Mono-small model data calculated (2/10)')
    print(f'Cost time: {cost_time:.2f} s')
    phi_array_m, e_array_m = get_model_data_array(epsilon_mono, nu_mono, d1_array, th_array_m, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n Mono-medium model data calculated (3/10)')
    print(f'Cost time: {cost_time:.2f} s')
    phi_array_l, e_array_l = get_model_data_array(epsilon_mono, nu_mono, d1_array, th_array_l, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n Mono-large model data calculated (4/10)')
    print(f'Cost time: {cost_time:.2f} s')

    rslt_dict_mono = {
        'd1': d1_array,
        'phi_s': phi_array_s,
        'phid_s': np.degrees(phi_array_s),
        'e_s': e_array_s,
        'phi_m': phi_array_m,
        'phid_m': np.degrees(phi_array_m),
        'e_m': e_array_m,
        'phi_l': phi_array_l,
        'phid_l': np.degrees(phi_array_l),
        'e_l': e_array_l,
    }
    np.savez('rslt_vs_d1_mono.npz', **rslt_dict_mono)

    bed_type['monodisperse'] = False  # switch to polydisperse bed

    phi_array_s, e_array_s = get_model_data_array(epsilon_2D, nu_2D, d1_array, th_array_s, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n 2D-small model data calculated (5/10)')
    print(f'Cost time: {cost_time:.2f} s')
    phi_array_m, e_array_m = get_model_data_array(epsilon_2D, nu_2D, d1_array, th_array_m, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n 2D-medium model data calculated (6/10)')
    print(f'Cost time: {cost_time:.2f} s')
    phi_array_l, e_array_l = get_model_data_array(epsilon_2D, nu_2D, d1_array, th_array_l, dist_params, bed_type)
    current_time = time.time()
    cost_time = current_time - start_time
    print(f'/n 2D-large model data calculated (7/10)')
    print(f'Cost time: {cost_time:.2f} s')

    rslt_dict_2D = {
        'd1': d1_array,
        'phi_s': phi_array_s,
        'phid_s': np.degrees(phi_array_s),
        'e_s': e_array_s,
        'phi_m': phi_array_m,
        'phid_m': np.degrees(phi_array_m),
        'e_m': e_array_m,
        'phi_l': phi_array_l,
        'phid_l': np.degrees(phi_array_l),
        'e_l': e_array_l,
    }
    np.savez('rslt_vs_d1_2D.npz', **rslt_dict_2D)

    #bed_type['three_D'] = True  # switch to 3D bed

    #phi_array_s, e_array_s = get_model_data_array(epsilon_3D, nu_3D, d1_array, th_array_s, dist_params, bed_type)
    #current_time = time.time()
    #cost_time = current_time - start_time
    #print(f'/n 3D-fine model data calculated (8/10)')
    #print(f'Cost time: {cost_time:.2f} s')
    #phi_array_m, e_array_m = get_model_data_array(epsilon_3D, nu_3D, d1_array, th_array_m, dist_params, bed_type)
    #current_time = time.time()
    #cost_time = current_time - start_time
    #print(f'/n 3D-medium model data calculated (9/10)')
    #print(f'Cost time: {cost_time:.2f} s')
    #phi_array_l, e_array_l = get_model_data_array(epsilon_3D, nu_3D, d1_array, th_array_l, dist_params, bed_type)
    #current_time = time.time()
    #cost_time = current_time - start_time
    #print(f'/n 3D-coarse model data calculated (10/10)')
    #print(f'Cost time: {cost_time:.2f} s')

    #rslt_dict_3D = {
    #    'd1': d1_array,
    #    'phi_s': phi_array_s,
    #    'phid_s': np.degrees(phi_array_s),
    #    'e_s': e_array_s,
    #    'phi_m': phi_array_m,
    #    'phid_m': np.degrees(phi_array_m),
    #    'e_m': e_array_m,
    #    'phi_l': phi_array_l,
    #    'phid_l': np.degrees(phi_array_l),
    #    'e_l': e_array_l,
    #}
    #np.savez('rslt_vs_d1_3D.npz', **rslt_dict_3D)
else:
    # Load results
    rslt_dict_mono = np.load('rslt_vs_d1_mono.npz')
    rslt_dict_2D = np.load('rslt_vs_d1_2D.npz')
    #rslt_dict_3D = np.load('rslt_vs_d1_3D.npz')

# Draw rebound angle data
plt.figure(1, figsize=(8, 6))
l1=plt.semilogx(rslt_dict_mono['d1']/d2_mid, rslt_dict_mono['phid_s'], 'r:', label='Model mono small')
l2=plt.semilogx(rslt_dict_mono['d1']/d2_mid, rslt_dict_mono['phid_m'], 'g:', label='Model mono medium')
l3=plt.semilogx(rslt_dict_mono['d1']/d2_mid, rslt_dict_mono['phid_l'], 'b:', label='Model mono large')
plt.semilogx(rslt_dict_2D['d1']/d2_mid, rslt_dict_2D['phid_s'], 'r--', label='Model 2D small')
plt.semilogx(rslt_dict_2D['d1']/d2_mid, rslt_dict_2D['phid_m'], 'g--', label='Model 2D medium')
plt.semilogx(rslt_dict_2D['d1']/d2_mid, rslt_dict_2D['phid_l'], 'b--', label='Model 2D large')
#plt.plot(rslt_dict_3D['d1']/d2_mid, rslt_dict_3D['phid_s'], 'r-', label='Model 3D small')
#plt.plot(rslt_dict_3D['d1']/d2_mid, rslt_dict_3D['phid_m'], 'g-', label='Model 3D medium')
#plt.plot(rslt_dict_3D['d1']/d2_mid, rslt_dict_3D['phid_l'], 'b-', label='Model 3D large')

plt.xlim(d1_min/d2_mid, d1_max/d2_mid)
#plt.ylim(0, 60)
plt.xlabel('$d_1/d_{50}$', fontsize=label_size)
plt.ylabel('$\\overline{\\theta \'}$(degree)', fontsize=label_size)
plt.xticks(fontsize=ticks_size)
plt.yticks(fontsize=ticks_size)
plt.legend([l1, l2, l3], ['$\\theta = 5 \\degree$', '$\\theta = 45 \\degree$', '$\\theta = 85 \\degree$'], fontsize=ticks_size, loc='best')

# Draw restitution coefficient data
plt.figure(2, figsize=(8, 6))
l1=plt.semilogx(rslt_dict_mono['d1']/d2_mid, rslt_dict_mono['e_s'], 'r:', label='Model mono fine')
plt.semilogx(rslt_dict_mono['d1']/d2_mid, rslt_dict_mono['e_m'], 'g:', label='Model mono medium')
plt.semilogx(rslt_dict_mono['d1']/d2_mid, rslt_dict_mono['e_l'], 'b:', label='Model mono coarse')
l2=plt.semilogx(rslt_dict_2D['d1']/d2_mid, rslt_dict_2D['e_s'], 'r--', label='Model 2D fine')
plt.semilogx(rslt_dict_2D['d1']/d2_mid, rslt_dict_2D['e_m'], 'g--', label='Model 2D medium')
plt.semilogx(rslt_dict_2D['d1']/d2_mid, rslt_dict_2D['e_l'], 'b--', label='Model 2D coarse')
#l3=plt.plot(rslt_dict_3D['d1']/d2_mid, rslt_dict_3D['e_s'], 'r-', label='Model 3D fine')
#plt.plot(rslt_dict_3D['d1']/d2_mid, rslt_dict_3D['e_m'], 'g-', label='Model 3D medium')
#plt.plot(rslt_dict_3D['d1']/d2_mid, rslt_dict_3D['e_l'], 'b-', label='Model 3D coarse')

plt.xlim(d1_min/d2_mid, d1_max/d2_mid)
#plt.ylim(0.4, 0.65)
plt.xlabel('$d_1/d_{50}$', fontsize=label_size)
plt.ylabel('$\\overline{e}$', fontsize=label_size)
plt.xticks(fontsize=ticks_size)
plt.yticks(fontsize=ticks_size)
plt.legend([l1, l2], ['2D_M', '2D_P'], fontsize=ticks_size, loc='best')

## Draw vertical restitution coefficient data
#plt.figure(3, figsize=(8, 6))
#ez_array_s_mono = rslt_dict_mono['e_s'] * np.sin(rslt_dict_mono['phi_s']) / np.sin(th_small)
#ez_array_m_mono = rslt_dict_mono['e_m'] * np.sin(rslt_dict_mono['phi_m']) / np.sin(th_medium)
#ez_array_l_mono = rslt_dict_mono['e_l'] * np.sin(rslt_dict_mono['phi_l']) / np.sin(th_large)
#ez_array_s_2D = rslt_dict_2D['e_s'] * np.sin(rslt_dict_2D['phi_s']) / np.sin(th_small)
#ez_array_m_2D = rslt_dict_2D['e_m'] * np.sin(rslt_dict_2D['phi_m']) / np.sin(th_medium)
#ez_array_l_2D = rslt_dict_2D['e_l'] * np.sin(rslt_dict_2D['phi_l']) / np.sin(th_large)
##ez_array_s_3D = rslt_dict_3D['e_s'] * np.sin(rslt_dict_3D['phi_s']) / np.sin(th_small)
##ez_array_m_3D = rslt_dict_3D['e_m'] * np.sin(rslt_dict_3D['phi_m']) / np.sin(th_medium)
##ez_array_l_3D = rslt_dict_3D['e_l'] * np.sin(rslt_dict_3D['phi_l']) / np.sin(th_large)
#
#plt.plot(rslt_dict_mono['d1']/d2_mid, ez_array_s_mono, 'r:', label='Model mono fine')
#plt.plot(rslt_dict_mono['d1']/d2_mid, ez_array_m_mono, 'g:', label='Model mono medium')
#plt.plot(rslt_dict_mono['d1']/d2_mid, ez_array_l_mono, 'b:', label='Model mono coarse')
#plt.plot(rslt_dict_2D['d1']/d2_mid, ez_array_s_2D, 'r--', label='Model 2D fine')
#plt.plot(rslt_dict_2D['d1']/d2_mid, ez_array_m_2D, 'g--', label='Model 2D medium')
#plt.plot(rslt_dict_2D['d1']/d2_mid, ez_array_l_2D, 'b--', label='Model 2D coarse')
##plt.plot(rslt_dict_3D['d1']/d2_mid, ez_array_s_3D, 'r-', label='Model 3D fine')
##plt.plot(rslt_dict_3D['d1']/d2_mid, ez_array_m_3D, 'g-', label='Model 3D medium')
##plt.plot(rslt_dict_3D['d1']/d2_mid, ez_array_l_3D, 'b-', label='Model 3D coarse')
#
#plt.xlim(d1_min/d2_mid, d1_max/d2_mid)
#plt.ylim(0, 2.5)
#plt.xlabel('$d_1/d_{50}$', fontsize=label_size)
#plt.ylabel('$\\overline{e_z}$', fontsize=label_size)
#plt.xticks(fontsize=ticks_size)
#plt.yticks(fontsize=ticks_size)

plt.show()
