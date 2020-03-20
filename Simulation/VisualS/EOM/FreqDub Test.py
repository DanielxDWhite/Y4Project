import matplotlib.pyplot as plt
import numpy as np 
import scipy.special as scp

B=0.6


J0  = scp.jv(0,B )
J1p = scp.jv(1,B )
J1m = -scp.jv(1,B )
J2p = scp.jv(2,B )
J2m = -scp.jv(2,B )
#print(J1m,J1p)
#print(J2m)

RF = 10
AOM = 100
w = 1000

f = w-AOM
F0 = f
F1p = f+RF 
F1m = f-RF
F2p = f+2*RF 
F2m = f-2*RF

# Freq and Efield together
A0 = [F0,J0]
print('A0',A0)

A1p = [F1p,J1p]
print('A1p',A1p)
A1m = [F1m,J1m]
A2p = [F2p,J2p]
A2m = [F2m,J2m]

# Adding and multiplying all permutations
# abcde = 1st thing   0,1p.. = 2nd thing

e0  = [A0[0]+A0[0],  A0[1]*A0[1]]
e1p = [A0[0]+A1p[0], A0[1]*A1p[1]]
e1m = [A0[0]+A1m[0], A0[1]*A1m[1]]
e2p = [A0[0]+A2p[0], A0[1]*A2p[1]]
e2m = [A0[0]+A2m[0], A0[1]*A2m[1]]
EE = [e0,e1p,e1m,e2p,e2m]

a0  = [A1p[0]+A0[0],  A1p[1]*A0[1]]
a1p = [A1p[0]+A1p[0], A1p[1]*A1p[1]]
a1m = [A1p[0]+A1m[0], A1p[1]*A1m[1]]
a2p = [A1p[0]+A2p[0], A1p[1]*A2p[1]]
a2m = [A1p[0]+A2m[0], A1p[1]*A2m[1]]
AA = [a0,a1p,a1m,a2p,a2m]

b0  = [A1m[0]+A0[0],  A1m[1]*A0[1]]
b1p = [A1m[0]+A1p[0], A1m[1]*A1p[1]]
b1m = [A1m[0]+A1m[0], A1m[1]*A1m[1]]
b2p = [A1m[0]+A2p[0], A1m[1]*A2p[1]]
b2m = [A1m[0]+A2m[0], A1m[1]*A2m[1]]
BB = [b0,b1p,b1m,b2p,b2m]

c0  = [A2p[0]+A0[0],  A2p[1]*A0[1]]
c1p = [A2p[0]+A1p[0], A2p[1]*A1p[1]]
c1m = [A2p[0]+A1m[0], A2p[1]*A1m[1]]
c2p = [A2p[0]+A2p[0], A2p[1]*A2p[1]]
c2m = [A2p[0]+A2m[0], A2p[1]*A2m[1]]
CC = [c0,c1p,c1m,c2p,c2m]

d0  = [A2m[0]+A0[0],  A2m[1]*A0[1]]
d1p = [A2m[0]+A1p[0], A2m[1]*A1p[1]]
d1m = [A2m[0]+A1m[0], A2m[1]*A1m[1]]
d2p = [A2m[0]+A2p[0], A2m[1]*A2p[1]]
d2m = [A2m[0]+A2m[0], A2m[1]*A2m[1]]
DD = [d0,d1p,d1m,d2p,d2m]

for j in range(5): 
    plt.scatter(EE[j][0],EE[j][1],s=100,alpha=0.3) 
    plt.scatter(AA[j][0],AA[j][1],s=100,alpha=0.3)
    plt.scatter(BB[j][0],BB[j][1],s=100,alpha=0.3)
    plt.scatter(CC[j][0],CC[j][1],s=100,alpha=0.3)
    plt.scatter(DD[j][0],DD[j][1],s=100,alpha=0.3)

XX = []
Dup = []

BBig = np.array([a0,a1p,a1m,a2p,a2m,b0,b1p,b1m,b2p,b2m,c0,c1p,c1m,c2p,c2m,d0,d1p,d1m,d2p,d2m,e0,e1p,e1m,e2p,e2m])
Flist = BBig[:,0]
Fplot = []

for k in range(len(BBig)):
    '''test for dups - get index symm - choose 1st'''
    x = np.where(BBig[:,0] == Flist[k])[0]
    xx = x[0]
    XX.append(x[0])
XX_ = list(dict.fromkeys(XX)) # Dup Prepended list
print('XX_ =',XX_)


for i in range(len(XX_)):
    '''Takes the 1st index & gives freq'''
    Fplot.append(Flist[XX_[i]])
Jsum = []

X0 = []
X0_ = []
for I in range(len(Flist)):
    x0 = np.where(BBig[:,0] == Flist[I])[0]

    X0.append(x0)
print('X0 ----', X0)
print(np.shape(X0))
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
    for k in range(len(X0_U[j])):
    
        y = BBig[X0_U[j][k]][1]
        Yy.append([y])

    Y = np.sum(Yy)
    Jsum.append(Y)
    Yy=[]
    
print("J's = ", Jsum)
print('Freqs =',Fplot)

plt.scatter(Fplot,Jsum,s=800,c='y', alpha=0.4,)

plt.xlabel('Frequency (AU)',size=18)
plt.ylabel('Electric Field / Bessel Value (AU)',size=18)
plt.title('EOM to Frequency Doubler',size=20)
plt.grid()
plt.show()