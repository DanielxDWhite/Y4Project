import numpy as np 
import matplotlib.pyplot as plt
import scipy.constants as sc


Xi = 25
E0 = 1

""" Physical & Atomic Constants """
kb=sc.Boltzmann
mu0 = sc.mu_0
muB = 9.2740099*10**-24
u=sc.proton_mass
hbar=sc.hbar
c=sc.c
pi=np.pi 
e=sc.e
M=87*u
wab=2*pi*384.23e12
G=38.11e6
Z =337
dip= 3.485e-29
##d = G*Xi
''' Variable Orgy '''
Rabi = dip*E0/hbar
IoIs = 2*Rabi**2/G**2
IrE = c*8.85e-12/2*E0**2/10000   # This is W/cm^2
Fplot = [2.41418817e15, 2.41418821e15, 2.41418809e15]
#d = np.subtract(wab, Fplot)
i=0
zs=[]
vs=[]
ts=[]
a=0

d = [10e6,30e6,50e6,70e6,90e6]
w = np.subtract(wab , d)

v = np.linspace(0,100,500)
def dv(v,F,D):
    Lambda=2*pi*c/F
    k = 2*pi/Lambda
    'Incremental Acceleration'
    O = F/(2*pi*c)
    c1 =1+IoIs+4*D**2/G**2
    c2 = O*8/G**2*D
    c3 = 4*O**2/G**2  
    rhoaa = + IoIs**2/(c1-c2*v+c3*v**2)   
    return rhoaa*hbar*k*G/M
print(dv(100,20,2))

L = len(w)
for i in range(L):
    plt.plot(v, dv(v,w[i],d[i]),linewidth=10,label=' {}$MHz$'.format(int(d[i]/1000000)))
S = []
'''
for i in range(L):
    S.append(dv(v,w[i],d[i]))
S = np.reshape(S,(len(w),len(v)))
#print(S)
Ss = []
for i in range(L):
    Sim = np.argmax(S[i])
    Ss = np.append(Ss,Sim)
print(Ss)
Vv = []
for i in range(L):
    Vv.append(v[int(Ss[i])])
print(Vv)
    
for i in range(L):
    plt.axvline(Vv[i])
    '''
plt.legend()
#plt.title('$Force$ $on$ $Rb$ $Atom$ $v$ $Velocity$ $for$ $different$ $detunings$', fontsize=20)
plt.ylabel('$\propto$ $Force$', fontsize=23)
plt.xlabel('$Velocity$ / $ms^{-1}$' , fontsize=23)
plt.show()