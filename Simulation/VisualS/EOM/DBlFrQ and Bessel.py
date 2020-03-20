# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 17:31:36 2020

@author: white
"""

import numpy as np 
import scipy.stats as sts
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.special as scp
import timeit
start = timeit.default_timer()
plt.close('all')

fig = plt.figure()
ax1 = plt.subplot2grid((1,3), (0,0), colspan=2,rowspan=2)
ax2 = plt.subplot2grid((1,3), (0,2) )
ax1.set_axisbelow(True)

# Turn on the minor TICKS, which are required for the minor GRID
ax1.minorticks_on()

# Customize the major grid
ax1.grid(which='major', linestyle='-', linewidth='0.5', color='red')
# Customize the minor grid
ax1.grid(which='minor', linestyle=':', linewidth='0.5', color='black')

E0 = 900
Beta = 0.55

RF = 5e6                       # EOM's RF input 
AOM = 335e6                     # "AOM"'s detuning 

'''
I M P O R T E D    C S V  O F   C A V I T Y 
'''

#  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  #
'''                E O M   &  F r e q  D u b    t i m e                   '''
#  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  #

G = 38.11e6                                  
wab = 2*np.pi*384.23e12    # Freq of transation

beta = np.linspace(0,4,80)


def EOM(RF, Beta): # [][1=Freq  0=E],AOM
    ' [0,1,2,3,4] = Center High1stSideBand Low1stSideBand High2ndSideBand... '
    Ei0 = scp.jv(0,Beta)
    Ei1 = scp.jv(1,Beta)
    Ei2 = scp.jv(2,Beta)
    Ei3 = scp.jv(3,Beta)
    F0 = wab/2 
    Fp1, Fm1 = wab/2+RF, wab/2-RF
    Fp2, Fm2 = wab/2+2*RF, wab/2-2*RF
    Fp3, Fm3 = wab/2+3*RF, wab/2-3*RF
    return [Ei0,F0],[Ei1,Fp1],[-Ei1,Fm1],[Ei2,Fp2],[Ei2,Fm2],[Ei3,Fp3],[-Ei3,Fm3]
#   -- EOM[band][E,F] --

ax2.plot(beta, abs(EOM(RF, beta)[0][0]),linewidth=9,c='darkviolet')
ax2.plot(beta, abs(EOM(RF, beta)[2][0]),linewidth=7,c='mediumorchid')
ax2.plot(beta, abs(EOM(RF, beta)[4][0]),linewidth=5,c='violet')
ax2.plot(beta, abs(EOM(RF, beta)[6][0]),linewidth=3,c='indigo')
ax2.axvline(Beta,c='k', linestyle='dashed')
ax2.axhline(abs(EOM(RF, Beta)[0][0]),linewidth=8,c='darkviolet', linestyle='dotted')
ax2.axhline(abs(EOM(RF, Beta)[2][0]),linewidth=6,c='mediumorchid', linestyle='dotted')
ax2.axhline(abs(EOM(RF, Beta)[4][0]),linewidth=4,c='violet', linestyle='dotted')
ax2.axhline(abs(EOM(RF, Beta)[6][0]),linewidth=3,c='indigo', linestyle='dotted')

# Freq and Efield together
A0  = [EOM(RF, Beta)[0][1],EOM(RF, Beta)[0][0]]
A1p = [EOM(RF, Beta)[1][1],EOM(RF, Beta)[1][0]]
A1m = [EOM(RF, Beta)[2][1],EOM(RF, Beta)[2][0]]
A2p = [EOM(RF, Beta)[3][1],EOM(RF, Beta)[3][0]]
A2m = [EOM(RF, Beta)[4][1],EOM(RF, Beta)[4][0]]
A3p = [EOM(RF, Beta)[5][1],EOM(RF, Beta)[5][0]]
A3m = [EOM(RF, Beta)[6][1],EOM(RF, Beta)[6][0]]


# Adding and multiplying all permutations
# abcde = 1st thing   0,1p.. = 2nd thing
e0  = [A0[0]+A0[0],  A0[1]*A0[1]]
e1p = [A0[0]+A1p[0], A0[1]*A1p[1]]
e1m = [A0[0]+A1m[0], A0[1]*A1m[1]]
e2p = [A0[0]+A2p[0], A0[1]*A2p[1]]
e2m = [A0[0]+A2m[0], A0[1]*A2m[1]]
e3p = [A0[0]+A3p[0], A0[1]*A3p[1]]
e3m = [A0[0]+A3m[0], A0[1]*A3m[1]]
EE = [e0,e1p,e1m,e2p,e2m,e3p,e3m]

a0  = [A1p[0]+A0[0],  A1p[1]*A0[1]]
a1p = [A1p[0]+A1p[0], A1p[1]*A1p[1]]
a1m = [A1p[0]+A1m[0], A1p[1]*A1m[1]]
a2p = [A1p[0]+A2p[0], A1p[1]*A2p[1]]
a2m = [A1p[0]+A2m[0], A1p[1]*A2m[1]]
a3p = [A1p[0]+A3p[0], A1p[1]*A3p[1]]
a3m = [A1p[0]+A3m[0], A1p[1]*A3m[1]]
AA = [a0,a1p,a1m,a2p,a2m,a3p,a3m]

b0  = [A1m[0]+A0[0],  A1m[1]*A0[1]]
b1p = [A1m[0]+A1p[0], A1m[1]*A1p[1]]
b1m = [A1m[0]+A1m[0], A1m[1]*A1m[1]]
b2p = [A1m[0]+A2p[0], A1m[1]*A2p[1]]
b2m = [A1m[0]+A2m[0], A1m[1]*A2m[1]]
b3p = [A1m[0]+A3p[0], A1m[1]*A3p[1]]
b3m = [A1m[0]+A3m[0], A1m[1]*A3m[1]]
BB = [b0,b1p,b1m,b2p,b2m,b3p,b3m]

c0  = [A2p[0]+A0[0],  A2p[1]*A0[1]]
c1p = [A2p[0]+A1p[0], A2p[1]*A1p[1]]
c1m = [A2p[0]+A1m[0], A2p[1]*A1m[1]]
c2p = [A2p[0]+A2p[0], A2p[1]*A2p[1]]
c2m = [A2p[0]+A2m[0], A2p[1]*A2m[1]]
c3p = [A2p[0]+A3p[0], A2p[1]*A3p[1]]
c3m = [A2p[0]+A3m[0], A2p[1]*A3m[1]]
CC = [c0,c1p,c1m,c2p,c2m,c3p,c3m]

d0  = [A2m[0]+A0[0],  A2m[1]*A0[1]]
d1p = [A2m[0]+A1p[0], A2m[1]*A1p[1]]
d1m = [A2m[0]+A1m[0], A2m[1]*A1m[1]]
d2p = [A2m[0]+A2p[0], A2m[1]*A2p[1]]
d2m = [A2m[0]+A2m[0], A2m[1]*A2m[1]]
d3p = [A2m[0]+A3p[0], A2m[1]*A3p[1]]
d3m = [A2m[0]+A3m[0], A2m[1]*A3m[1]]
DD = [d0,d1p,d1m,d2p,d2m,d3p,d3m]

f0  = [A3m[0]+A0[0],  A3p[1]*A0[1]]
f1p = [A3m[0]+A1p[0], A3p[1]*A1p[1]]
f1m = [A3m[0]+A1m[0], A3p[1]*A1m[1]]
f2p = [A3m[0]+A2p[0], A3p[1]*A2p[1]]
f2m = [A3m[0]+A2m[0], A3p[1]*A2m[1]]
f3p = [A3m[0]+A3p[0], A3p[1]*A3p[1]]
f3m = [A3m[0]+A3m[0], A3p[1]*A3m[1]]

FF = [f0,f1p,f1m,f2p,f2m,f3p,f3m]
g0  = [A3m[0]+A0[0],  A3m[1]*A0[1]]
g1p = [A3m[0]+A1p[0], A3m[1]*A1p[1]]
g1m = [A3m[0]+A1m[0], A3m[1]*A1m[1]]
g2p = [A3m[0]+A2p[0], A3m[1]*A2p[1]]
g2m = [A3m[0]+A2m[0], A3m[1]*A2m[1]]
g3p = [A3m[0]+A3p[0], A3m[1]*A3p[1]]
g3m = [A3m[0]+A3m[0], A3m[1]*A3m[1]]
GG = [g0,g1p,g1m,g2p,g2m,g3p,g3m]


for j in range(7): 
    ax1.scatter(EE[j][0]-AOM ,EE[j][1],s=100,alpha=0.3) 
    ax1.scatter(AA[j][0]-AOM ,AA[j][1],s=100,alpha=0.3)
    ax1.scatter(BB[j][0]-AOM ,BB[j][1],s=100,alpha=0.3)
    ax1.scatter(CC[j][0]-AOM ,CC[j][1],s=100,alpha=0.3)
    ax1.scatter(DD[j][0]-AOM ,DD[j][1],s=100,alpha=0.3)
    ax1.scatter(FF[j][0]-AOM ,FF[j][1],s=100,alpha=0.3)
    ax1.scatter(GG[j][0]-AOM ,GG[j][1],s=100,alpha=0.3)

BBig = np.array([a0,a1p,a1m,a2p,a2m,b0,b1p,b1m,b2p,b2m,c0,c1p,c1m,c2p,c2m,d0,d1p,d1m,d2p,d2m,e0,e1p,e1m,e2p,e2m,f0,f1p,f1m,f2p,f2m,f3p,f3m,g0,g1p,g1m,g2p,g2m,g3p,g3m])
Flist = BBig[:,0]
XX = []         # Index of frequency values
for l in range(len(BBig)):
    '''test for dups - get index symm - choose 1st'''
    x = np.where(BBig[:,0] == Flist[l])[0]
    xx = x[0]
    XX.append(x[0])
XX_0 = list(dict.fromkeys(XX)) # Dup Prepended list
XX_ = np.sort(XX_0)
Fplot_ = []
for i in range(len(XX_)):
    '''Takes the 1st index & gives freq'''
    Fplot_.append(Flist[XX_[i]])
    Fplot = np.subtract(Fplot_,AOM)
print('XXXXXXX',XX_)
Jsum = []
X0 = []
X0_ = []
for I in range(len(Flist)):
    x0 = np.where(BBig[:,0] == Flist[I])[0]
    X0.append(x0)
u =[]
for i in range(len(X0)):
    U=[]
    u = (X0[i])
    for j in range(len(X0[i])):
        U.append(u[j])
    X0_.append(U) 
X0_U = np.unique(X0_)  # Index Prepended list

Yy=[]
for j in range(len(Fplot)):
    for l in range(len(X0_U[j])):
        y = BBig[X0_U[j][l]][1]
        Yy.append([y])
    Y = np.sum(Yy)
    Jsum.append(Y)
    Yy=[]

ax1.scatter(Fplot,np.square(Jsum),s=1500,c='y',alpha=0.4)
ax1.scatter(Fplot,Jsum,s=200,c='G',alpha=0.3)

ax1.set_xlabel('MHz',size=10)
ax1.set_ylabel('Electric Field / E0',size=15)
ax2.set_xlabel(r'$ \beta $',size=22)
ax2.set_ylabel('Bessel Value',size=13)
ax1.grid(which='major')
#ax1.plot(Fplot,np.square(Jsum),c='y')

for i in range(len(Fplot)):
    ax1.axvline(Fplot[i],alpha=0.3)
#ax1.axvline(wab, c='pink',linewidth=9)
#ax1.set_xlim(left=min(Fplot)-2*RF, right=wab+RF*2)
ax1.set_xlim(left=2*RF+min(Fplot)-5*RF, right=2*RF+min(Fplot)+RF*12)

####################################################
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
Rabiii = dip*max(Jsum)/hbar
IoIs = 2*Rabiii**2/G**2
IrE = c*8.85e-12/2*E0**2/10000   # This is W/cm^2

d = np.subtract(wab, Fplot)
w = wab - d

i=0
zs=[]
vs=[]
ts=[]
a=0
L = len(Fplot)
#print(mu0)
#print()
#Ir=power/(a**2*pi)
v = np.linspace(0,600,100)
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

for i in range(L):
    ax1.plot(v, dv(v,Fplot[i],d[i]))
S = []
for i in range(L):
    S.append(dv(v,Fplot[i],d[i]))
S = np.reshape(S,(len(Fplot),len(v)))
#print(S)
Ss = []
for i in range(L):
    Sim = np.argmax(S[i])
    Ss = np.append(Ss,Sim)

Vv = []
for i in range(L):
    Vv.append(v[int(Ss[i])])
    
for i in range(L):
    ax1.axvline(Vv[i])
ax1.grid()

plt.title('Force on Rb Atom v Velocity for different detunings', fontsize=12)
plt.ylabel('∝ Force', fontsize=12)
plt.xlabel('Velocity m/s', fontsize=12)
ax1.axhline(0,c='k')

print(Fplot)

plt.show()