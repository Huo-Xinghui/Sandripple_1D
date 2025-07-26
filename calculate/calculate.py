import numpy as np
import matplotlib.pyplot as plt

d1 = 3e-4
d2 = 3e-4
epsilon = 0.78
nu = -0.13
dd1 = d1/(0.5*(d1+d2))
dd2 = d2/(0.5*(d1+d2))
mu = epsilon*dd1**3/(dd1**3 + epsilon*dd2**3)
alpha = (1 + epsilon)/(1 + mu) - 1
beta = 1 - (2 / 7) * (1 - nu)/(1 + mu)
#print('alpha =', alpha)
#print('beta =', beta)

Q1 = 0.6939146e-02
Q2 = 0.7682588E-02
rhop = 2650
rho = 1.263
s = rhop/rho
ghat = 9.8*(1 - 1/s)
Qstar1 = Q1/(rhop*d1*(s*ghat*d1)**0.5)
Qstar2 = Q2/(rhop*d2*(s*ghat*d2)**0.5)
print('Qstar =', Qstar1, Qstar2)

d13 = np.linspace(0.1, 1, 10)
d23 = 1
thetac = np.zeros_like(d13)
for i, d in enumerate(d13):
    cos1 = (1 + d**2 - d23**2)/(2*d)
    cos2 = (1 + d23**2 - d**2)/(2*d23)
    thetac[i] = np.arccos(cos1) + np.arccos(cos2) - np.pi/2
    thetac[i] = thetac[i] * 180/np.pi
plt.plot(d13, thetac, 'o')
plt.show()