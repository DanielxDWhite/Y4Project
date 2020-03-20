# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 13:59:01 2019

@author: Daniel White
"""

import matplotlib.pyplot as plt
import numpy as np

x = 10

t = np.array(range(-x, x))

d = 300             # Detuning Value, delta    
O = 50 + 40j     # Rabi Frequnecy (COMPLEX), Omega

dt = d * t * 0.5 
#H = [[-d,O],
#    [np.conj(O),d]]
def square(list):
    return [i**2 for i in list]
              

E1 = np.sqrt(np.abs(O)**2 + square(dt))
E2 = - np.abs(np.sqrt(np.abs(O)**2 + square(dt)))

Lo1 = [1] * len(t)
E1p = np.multiply(E1,Lo1)
E2p = np.multiply(E2,Lo1)
dtp = np.multiply(d, Lo1)

plt.plot(t, E1p, label = '$E_+$')
plt.plot(t, E2p, label = '$E_-$')
plt.plot(t, dtp, label = 'Detuning', color = 'grey', linewidth='0.4')
plt.plot(t, -dtp, color = 'grey', linewidth='0.4')

plt.xlabel('Time',fontsize = 23)
plt.ylabel('Energy', fontsize = 23)
plt.title('Energies of states allowed in 2x2 simple light dressed Hamiltonian',fontsize = 23)
plt.legend()
plt.show()
print (O)




