# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 13:59:01 2019

@author: Daniel White
"""

import matplotlib.pyplot as plt
import numpy as np

x = 10

t = np.linspace(0,100,100)

d = 7            # Detuning Value, delta    
def O(i):
    Oo = 299*i
    return float(Oo)              # Rabi Frequnecy (COMPLEX), Omega

dt = d * t * 0.5 
#H = [[-d,O],
#    [np.conj(O),d]]
def H2(i):
    Haux = np.array([[d,O],
                     [O,-d]])#33           #Ox/2 is from EdbifExp
    return np.linalg.eigvals(Haux)
    
E1, E2 = np.zeros((1,0)),np.zeros((1,0))
for i in range(int(-0.5*(len(t)),len(t)))):     # C U R R E N T L Y   F U C K E D     !!!!!!!!!!!!!!!!!!!!!!!!!!
    i1, i2 = np.split(H2(i), 2)
    
E1, E2 = np.append(E1, i1), np.append(E2, i2)

Lo1 = [1] * len(t)
E1p = np.multiply(E1,Lo1)
E2p = np.multiply(E2,Lo1)
dtp = np.multiply(d, Lo1)
plt.close('all')
plt.plot(t, E1, label = '$E_+$', color='blue')
plt.plot(t, E2, label = '$E_-$', color='blueviolet')
plt.plot(t, dtp, label = 'Detuning from 0', color = 'grey', linewidth='0.4')
plt.plot(t, -dtp, color = 'grey', linewidth='0.4')
plt.show
plt.xlabel('time')
plt.ylabel('Energy/hbar')
plt.title('Energies of states allowed in 2x2 simple light dressed Hamiltonian')
plt.legend(title='Detuning = {}  &  Rabi Freq = {}'.format(d,O))
print (O)




