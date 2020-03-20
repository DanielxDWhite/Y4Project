# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 17:59:19 2019

@author: Daniel White
"""

import matplotlib.pyplot as plt
import numpy as np


d = 1000         #Detuning, detla
O = np.sqrt(1.5)*d

#Ob, Or = O+d, O-d  # Rabi Freq, Omega ******* PLACEHOLDER ********** 

t = np.linspace(0,20,200)
dt = d*t*0.5
#A=np.zeros(len(t))
#for i in range(0,len(t)):
#    a=np.sin(t[i])
#    A[i]+=a
    
#Or = (O+d) # Rabi Freq, Omega ******* PLACEHOLDER ********** 
#Ob = O-d

E1,E2,E3,E4,E5,E6,E7 = np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0))

def H7(Or, Ob, d, t):
    Haux = np.array([[3*d,Or,0,0,0,0,0],
                     [Or,2*d,Ob,0,0,0,0],
                     [0,Ob,d,Or,0,0,0],
                     [0,0,Or,0,Ob,0,0],
                     [0,0,0,Ob,-d,Or,0],
                     [0,0,0,0,Or,-2*d,Ob],
                     [0,0,0,0,0,Ob,-3*d]])
    return np.linalg.eigvals(Haux)
    
    


for i in range(0,len(t)):   #                Or                    Ob            d
    i1,i2,i3,i4,i5,i6,i7 = np.split(H7(100*np.sin(i*9), 100*np.sin(90+i*8), 100, t), 7)
    E1 = np.append(E1, i1)
    E2 = np.append(E2, i2)
    E3 = np.append(E3, i3)
    E4 = np.append(E4, i4)
    E5 = np.append(E5, i5)
    E6 = np.append(E6, i6)
    E7 = np.append(E7, i7)
    
    

print (H7(np.cos(i*20), np.sin(i**4), 1000, t))
#(E1,E2,E3,E4,E5,E6,E7) = H7(Or, 2, 1000, t)
             
    
#H7(A, Ob, d)
#(E1,E2,E3,E4,E5,E6,E7) = np.linalg.eigvals(H7)

Lo1 = [1] * len(t)

#dtp = np.multiply(Lo1, d)

plt.plot(t, E1, label = '$E_1$', color='red')
plt.plot(t, E2, label = '$E_2$', color='crimson')
plt.plot(t, E3, label = '$E_3$', color='deeppink')
plt.plot(t, E4, label = '$E_4$', color='mediumorchid')
plt.plot(t, E5, label = '$E_5$', color='blueviolet')
plt.plot(t, E6, label = '$E_6$', color='blue')
plt.plot(t, E7, label = '$E_7$', color='midnightblue')

#plt.plot(t, dtp, label = 'Detuning', color = 'grey', linewidth='0.4')
#plt.plot(t, -dtp, color = 'grey', linewidth='0.4')


plt.show
plt.xlabel('time')
plt.ylabel('Energy/hbar')
plt.title('Energies of states allowed in 7x7 bichromatic dressed Hamiltonian')
plt.legend()



