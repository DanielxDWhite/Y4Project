# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 20:14:08 2019

@author: Daniel White
"""
import matplotlib.pyplot as plt
import numpy as np


n= 300
l_of_t = 100
t = np.linspace(0,l_of_t,n)          #time / x axis

def Or(i):
    '''Rabi Frequency of Red transition formula'''
    o_r = 15*np.sin(i*2*np.pi/n + np.pi/3)
    return float(o_r)
    return o_r
    
def Ob(i):
    '''Rabi Frequency of Blue transition formula'''
    ob = 15*np.sin(i*2*np.pi/n)
    return ob   
    return float(ob)

E, y0 = np.zeros((1,0)) ,np.zeros((1,0))     
    
for i in range(0,n):
    E = np.append(E, Or(i))

plt.plot(t, E)


'''Finding Min Index'''

Ep = np.square(E)
xW = np.where(Ep < 1)
len_of_x = len(xW[0])

Ex = np.zeros((0,1))    
for i in range(0, NoIntR):
    ax1.axvline(x = ((float(xWr[0][i])-n/2)*l_of_t/n), linestyle=':',color='red')

    EpR = np.square(Er)                                                             # Intercepts
EpB = np.square(Eb)        
xWb = np.where(EpR < Int_Thresh)
xWr = np.where(EpB < Int_Thresh)
NoIntR = len(xWr[0])
NoIntB = len(xWb[0])
for i in range(0, NoIntR):
    ax1.axvline(x = ((float(xWr[0][i])-n/2)*l_of_t/n), linestyle=':',color='red')
    ax2.axvline(x = ((float(xWr[0][i])-n/2)*l_of_t/n), linestyle=':',color='red')
for i in range(0, NoIntB):
    ax1.axvline(x = ((float(xWb[0][i])-n/2)*l_of_t/n), linestyle=':',color='b')
    ax2.axvline(x = ((float(xWb[0][i])-n/2)*l_of_t/n), linestyle=':',color='b')

print(Ex)
print(Ep[177])




    
#plt.axvline(x = Ex[0])
'''Mass labelling'''

def j(i):
    j = np.array((0,len(Ex[0])))
    return np.append(j, i)
    





print()

        

    







