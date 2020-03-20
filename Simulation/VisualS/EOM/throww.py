import numpy as np 
import matplotlib.pyplot as plt
import scipy.constants as sc

v = np.linspace(0,1000,50)
Xi = 25
E0 = 130

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
d = np.subtract(wab, Fplot)
w = wab - d


i=0
zs=[]
vs=[]
ts=[]
a=0
print(Fplot[0])
print(v)
#print(mu0)
#print()
#Ir=power/(a**2*pi)
v = np.linspace(0,1000,100)
def dv(v,F,D):
    Lambda=2*pi*c/F
    k = 2*pi/Lambda
    'Incremental Acceleration'
    O = F/(2*pi*c)
    c1 =1+IoIs+4*D**2/G**2
    c2 = O*8/G**2*D
    c3 = 4*O**2/G**2  
    rhoaa =  -IoIs**2/(c1+c2*v+c3*v**2) + IoIs**2/(c1-c2*v+c3*v**2)   
    return rhoaa*hbar*k*G/M
print(dv(100,20,2))

L = len(Fplot)
for i in range(L):
    plt.plot(v, dv(v,Fplot[i],d[i]))
S = []
for i in range(L):
    S.append(dv(v,Fplot[i],d[i]))
S = np.reshape(S,(len(Fplot),len(v)))
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

plt.title('Force on Rb Atom v Velocity for different detunings', fontsize=20)
plt.ylabel('∝ Force', fontsize=20)
plt.xlabel('Velocity m/s', fontsize=20)
plt.show()