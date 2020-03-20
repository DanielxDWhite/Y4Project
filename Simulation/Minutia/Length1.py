print('Length calculation')
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

Vmax = 500
Vmin = 20

m = 1.44316060e-25
h = sc.hbar
wvL = 780.241209686e-9
#d = 15*10e6
d = np.linspace(10e7,10e8,1000)

k = 2*np.pi/wvL
xi = 2*d/np.pi

L = m*(Vmax**2-Vmin**2)/(h*k*xi)

plt.plot(d/10**6,L)
plt.title("Tube Length as f(Detuning) from velocity {} --> {} ms-'".format(Vmax,Vmin))
plt.xlabel('Detuning / MHz')
plt.ylabel('Length')