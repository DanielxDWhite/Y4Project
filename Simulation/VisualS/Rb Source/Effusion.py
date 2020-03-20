"""
= Effusion Test       https://en.wikipedia.org/wiki/Effusion
"""
import math

import numpy as np
import scipy.constants as sc

Col_Gap = 0.004
C = 100     #Temp in Celcius    
   # 
T = 273 + C
P = 133.322*(3e-7/(273+25))*T
kb = sc.k
M = 87 * sc.proton_mass
A = np.pi*(Col_Gap/2)**2

Q = P*A/((2*np.pi*M*kb*T)**0.5)
def OoM(number):
    return math.floor(math.log(number, 10))

print( 'Flow Rate of Particles = {}  [{}]'.format(round(Q,3), OoM(Q)) )