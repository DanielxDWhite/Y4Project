import matplotlib.pyplot as plt
import numpy as np
c = 3*10**8

#d = 30*10**6
d = np.linspace(0.5 * 10**7, 11**7, 15)

f0 = 382.230484468*10**12
f = f0 + d

a = f/f0
xx = 1-a
xy = np.negative(2)
xz = 1+a

disc = xy**2-4*xz*xx
Bp = (xy+np.sqrt(disc))/2*xx
Bm = (xy-np.sqrt(disc))/2*xx

Vp = c * Bp
Vm = c * Bm


print(Vp,Vm)

dplot = d*10**(-6)
Vpplot = Vp *10**(7)

plt.plot(dplot, Vpplot)
plt.xlabel('Detuning / MHz')
plt.ylabel('Velocity of Mirror $x10^{-7}$ / $ms^{-1}$')
plt.title('Doppler shift via a moving mirror')
plt.grid(color='grey', linestyle='-', linewidth=0.2)
plt.show()
plt.close










































