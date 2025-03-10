import numpy as np
import matplotlib.pyplot as plt

def evx_calculate(x):
    res = -(alpha + beta)*(1 - (x*np.sin(theta))**2)**(3/2) + x*np.cos(theta)*(-3*alpha + (alpha + beta)*(x*np.sin(theta))**2)
    return res

def evz_calculate(x):
    csc = 1/np.sin(theta)
    cot = 1/np.tan(theta)
    res = alpha*x - 1/3*(alpha + beta)*(x**3*np.sin(theta)**2 + cot*csc*(1 - (x*np.sin(theta))**2)**(3/2))
    return res

def equation(x, theta, rhs):
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    return (x*sin_theta)**2 - cos_theta*x*np.sqrt(1 - (x*sin_theta)**2) - rhs

def find_xmax(theta, xmin, rhs, num_points=100):
    res_max = 1/np.sin(theta)
    res_min = xmin
    res_values = np.linspace(res_min, res_max, num_points)
    valid_res = []
    for res in res_values:
        if equation(res, theta, rhs) > 0:
            valid_res.append(res)
    return valid_res

angle = np.linspace(0.01, np.pi/2, num=100).tolist()
d1 = 3e-4
d2 = 3e-4
epsilon = 0.78
nu = -0.13
dd1 = 2/(1+d2/d1)
dd2 = 2/(1+d1/d2)
mu = ((dd2/dd1)**3+1/epsilon)**(-1)
alpha = (1 + epsilon)/(1 + mu) - 1
beta = 1 - (2 / 7) * (1 - nu)/(1 + mu)
e_list = []
ez_list = []
theta_out_list = []
for theta in angle:
    csc = 1/np.sin(theta)
    cot = 1/np.tan(theta)
    if 2*np.sin(theta) < dd2:
        xmin = csc - dd2
    else:
        xmin = cot*np.sqrt(1 - (dd2/2)**2) - dd2/2
    rhs = 1/(1 + beta/alpha)
    xmaxs = find_xmax(theta, xmin, rhs)
    xmax = xmaxs[0]
    evxmax = evx_calculate(xmax)
    evxmin = evx_calculate(xmin)
    evx = 1/(xmax - xmin)/3*(evxmax - evxmin)
    evzmax = evz_calculate(xmax)
    evzmin = evz_calculate(xmin)
    evz = np.sin(theta)/(xmax - xmin)*(evzmax - evzmin)
    e = np.sqrt(evx**2 + evz**2)
    ez = evz*(np.sin(theta))**(-1)
    theta_out = np.arcsin(ez*np.sin(theta)/e)
    e_list.append(e)
    ez_list.append(ez)
    theta_out_list.append(theta_out)
angle_deg = [i*180/np.pi for i in angle]
angle_out_deg = [i*180/np.pi for i in theta_out_list]
#plt.plot(angle_deg, e_list, 'r-', label='e')
plt.plot(angle_deg, ez_list, 'b-', label='ez')
#plt.plot(angle_deg, angle_out_deg, 'g-', label='theta_out')
plt.xlabel('Angle')
plt.ylabel('Value')
#plt.xlim(0, 90)
#plt.ylim(0, 90)
plt.grid(True)
plt.legend()
plt.show()