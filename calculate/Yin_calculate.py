import numpy as np
import matplotlib.pyplot as plt

theta_in = 11.5
theta_bed_list = [0, 5, 10, 12.5, 15]
A_list = [0.29, 0.29, 0.29, 0.29, 0.29]
C_list = [1.48, 1.74, 1.88, 1.69, 1.53]
L_list = [3.3518e-4, 2.9602e-4, 2.7397e-4, 2.8809e-4, 2.11e-4]
din = 3e-4
theta_bed_array = np.array(theta_bed_list)
A_array = np.array(A_list)
C_array = np.array(C_list)
L_array = np.array(L_list)
tan_phi_array = A_array + np.sin(C_array)/((din/L_array)**3 - np.cos(C_array))
phi_array = np.arctan(tan_phi_array)
phi_array = np.rad2deg(phi_array)
plt.figure(figsize=(8, 6))
plt.plot(theta_bed_array, phi_array, marker='o')
plt.xlabel(r'$\theta_{bed}$ (degrees)')
plt.ylabel(r'$\tan(\phi)$')
plt.show()