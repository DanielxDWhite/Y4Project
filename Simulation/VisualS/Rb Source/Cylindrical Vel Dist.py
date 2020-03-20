'''
Cylindrical Velocity Distriution
'''
import numpy as np 
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.integrate as si
import scipy.optimize as so

T=300
m = 10**-20
vx = 300
def Cyl_V_D(vx, T):
    N = (m/(2*np.pi*sc.k*T))**(3/2)
    exp = -m*(vx**2)/(2*sc.k*T)
    MB = 2*N*np.pi*vx*np.exp(exp)
    return MB
    
def Int(x):   
    Int = si.quad(Cyl_V_D, 0, x, args=(T))[0]-0.8
    return Int

vr = so.fsolve(Int, 0.5)

print(vr)

