'''
~Created on Thu Oct 10 17:59:19 2019
@author: Daniel White
'''

import matplotlib.pyplot as plt
import numpy as np


d = 1000         #Detuning, detla
O = np.sqrt(1.5)*d



t = np.linspace(0,100,70)
dt = d*t*0.5


E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13 = np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0))

def H13(Or, Ob, d, t):
    Haux = np.array([[6*d,Ob/2,0,0,0,0,0,0,0,0,0,0,0],
                     [Ob/2,5*d,Or/2,0,0,0,0,0,0,0,0,0,0],
                     [0,Or/2,4*d,Ob/2,0,0,0,0,0,0,0,0,0],
                     [0,0,Ob/2,3*d,Or/2,0,0,0,0,0,0,0,0],
                     [0,0,0,Or/2,2*d,Ob/2,0,0,0,0,0,0,0],
                     [0,0,0,0,Ob/2,d,Or/2,0,0,0,0,0,0],
                     [0,0,0,0,0,Or/2,0,Ob/2,0,0,0,0,0],
                     [0,0,0,0,0,0,Ob/2,-d,Or/2,0,0,0,0],
                     [0,0,0,0,0,0,0,Or/2,-2*d,Ob/2,0,0,0],
                     [0,0,0,0,0,0,0,0,Ob/2,-3*d,Or/2,0,0],
                     [0,0,0,0,0,0,0,0,0,Or/2,-4*d,Ob/2,0],
                     [0,0,0,0,0,0,0,0,0,0,Ob/2,-5*d,Or/2],
                     [0,0,0,0,0,0,0,0,0,0,0,Or/2,-6*d]])   #Ox/2 is from EdbifExp
    return np.linalg.eigvals(Haux)
    
#print( H13(1,2,3,4) )


for i in range(0,len(t)):   #                                  Or                    Ob            d
    i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13 = np.split(H13(28*np.sin(i*3), 28*np.sin(19+i*5), 15, t), 13)
    E1 = np.append(E1, i1)
    E2 = np.append(E2, i2)
    E3 = np.append(E3, i3)
    E4 = np.append(E4, i4)
    E5 = np.append(E5, i5)
    E6 = np.append(E6, i6)
    E7 = np.append(E7, i7)
    E8 = np.append(E8, i8)
    E9 = np.append(E9, i9)
    E10 = np.append(E10, i10)
    E11 = np.append(E11, i11)
    E12 = np.append(E12, i12)
    E13 = np.append(E13, i13)
    

Lo1 = [1] * len(t)
th = 3
#plt.plot(t, E1, label = '$E_1$', color='brown')
#plt.plot(t, E2, label = '$E_2$', color='tomato')
plt.plot(t, E3, label = '$E_3$', color='sandybrown', linewidth=th)
plt.plot(t, E4, label = '$E_4$', color='gold',linewidth=th)
plt.plot(t, E5, label = '$E_5$', color='yellowgreen',linewidth=th)
plt.plot(t, E6, label = '$E_6$', color='seagreen',linewidth=th)
plt.plot(t, E7, label = '$E_7$', color='turquoise',linewidth=th)
plt.plot(t, E8, label = '$E_8$', color='deepskyblue',linewidth=th)
plt.plot(t, E9, label = '$E_9$', color='royalblue',linewidth=th)
plt.plot(t, E10, label = '$E_{10}$', color='mediumblue',linewidth=th)
plt.plot(t, E11, label = '$E_{11}$', color='indigo',linewidth=th)
#plt.plot(t, E12, label = '$E_{12}$', color='palevioletred')
#plt.plot(t, E13, label = '$E_{13}$', color='lightpink')


plt.show
plt.xlabel('Distance / a.u.', size = 20)
plt.ylabel('Energy / a.u.', size = 20)
plt.xticks(size = 14)
plt.yticks(size = 14)
#plt.title('Energies of states allowed in 13x13 bichromatic dressed Hamiltonian')
plt.legend()
plt.close