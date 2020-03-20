# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 17:11:37 2019

@author: Daniel White
"""

import matplotlib.pyplot as plt
import numpy as np


d = 1000         #Detuning, detla
O = np.sqrt(1.5)*d

Ob, Or = O+d, O-d  # Rabi Freq, Omega ******* PLACEHOLDER ********** 

x = 10
t = np.array(range(-x,x))
dt = d*t*0.5

H3 = np.array([[d, Or, 0],
               [Or, 0, Ob],
               [0, Ob, -d]])

(E1,E2,E3) = np.linalg.eigvals(H3)

Lo1 = [1] * len(t)
E1p = np.multiply(Lo1, E1)
E2p = np.multiply(Lo1, E2)
E3p = np.multiply(Lo1, E3)
dtp = np.multiply(Lo1, d)

plt.plot(t, E1p, label = '$E_1$', color='red')
plt.plot(t, E2p, label = '$E_2$', color='orangered')
plt.plot(t, E3p, label = '$E_3$', color='gold')
plt.plot(t, dtp, label = 'Detuning', color = 'grey', linewidth='0.4')
plt.plot(t, -dtp, color = 'grey', linewidth='0.4')
plt.show
plt.xlabel('time')
plt.ylabel('Energy/hbar')
plt.title('Energies of states allowed in 3x3 simple light dressed Hamiltonian')
plt.legend()
print (O)


