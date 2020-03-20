import numpy as np 
import scipy.stats as sts
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.special as scp
import timeit
start = timeit.default_timer()
plt.close('all')

fig = plt.figure()
ax1 = plt.subplot2grid((4,4), (0,1), colspan=3, rowspan=1)
ax2 = plt.subplot2grid((4,4), (1,1), colspan=3,sharex=ax1, rowspan=2)
ax5 = plt.subplot2grid((4,4), (3,1), colspan=3,sharex=ax1, rowspan=1)
ax3 = plt.subplot2grid((4,4), (0,0), rowspan=3)
ax4 = plt.subplot2grid((4,4), (3,0), rowspan=1)

E0 = 300
Beta = 0.3

RF = 5e6                       # EOM's RF input 
AOM = 35e6                     # "AOM"'s detuning 

aa = 0.15                       # Coil Radius
s = 0.11                        # Coil Separation
Curr = 0.6805                   # Current (16G/cm @ 0.6805)
z0 = 0.173                      # Position of MOT centre

h = 0.0001                     # Step Size
Natoms = 100                 # Atoms
Nsteps = 170                 # Steps

Col_Gap = 0.002                 # Pinhole Diameter
ymot = 0.008                    # Diameter of Cold Atoms


#  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  #
'''                E O M   &  F r e q  D u b    t i m e                   '''
#  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  #

G = 38.11e6                                  
wab = 2*np.pi*384.23e12    # Freq of transation

hbar = sc.hbar             # hbar
dip = 3.485e-29            # dipole moment
c = sc.c
u = sc.proton_mass         # Proton Mass
M = 87*u                   # Mass of 87Rb
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

ax4.plot(beta, abs(EOM(RF, beta)[0][0]),linewidth=9,c='darkviolet')
ax4.plot(beta, abs(EOM(RF, beta)[2][0]),linewidth=7,c='mediumorchid')
ax4.plot(beta, abs(EOM(RF, beta)[4][0]),linewidth=5,c='violet')
ax4.plot(beta, abs(EOM(RF, beta)[6][0]),linewidth=3,c='indigo')
ax4.axvline(Beta,c='k', linestyle='dashed')
ax4.axhline(abs(EOM(RF, Beta)[0][0]),linewidth=8,c='darkviolet', linestyle='dotted')
ax4.axhline(abs(EOM(RF, Beta)[2][0]),linewidth=6,c='mediumorchid', linestyle='dotted')
ax4.axhline(abs(EOM(RF, Beta)[4][0]),linewidth=4,c='violet', linestyle='dotted')
ax4.axhline(abs(EOM(RF, Beta)[6][0]),linewidth=3,c='indigo', linestyle='dotted')

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
    ax3.scatter(EE[j][0]-AOM ,EE[j][1],s=100,alpha=0.3) 
    ax3.scatter(AA[j][0]-AOM ,AA[j][1],s=100,alpha=0.3)
    ax3.scatter(BB[j][0]-AOM ,BB[j][1],s=100,alpha=0.3)
    ax3.scatter(CC[j][0]-AOM ,CC[j][1],s=100,alpha=0.3)
    ax3.scatter(DD[j][0]-AOM ,DD[j][1],s=100,alpha=0.3)
    ax3.scatter(FF[j][0]-AOM ,FF[j][1],s=100,alpha=0.3)
    ax3.scatter(GG[j][0]-AOM ,GG[j][1],s=100,alpha=0.3)

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

ax3.scatter(Fplot,np.square(Jsum),s=1000,c='y',alpha=0.4)
ax3.plot(Fplot,np.square(Jsum),c='y')
ax3.scatter(Fplot,Jsum,s=200,c='G',alpha=0.3)
ax3.set_xlabel('MHz',size=10)
ax3.set_ylabel('Electric Field / E0',size=15)
ax4.set_xlabel(r'$ \beta $',size=22)
ax4.set_ylabel('Bessel Value',size=13)
ax3.grid(which='major')
ax3.set_xlim(left=min(Fplot)-2*RF, right=wab+RF*2)
for i in range(len(Fplot)):
    ax3.axvline(Fplot[i],alpha=0.3)
ax3.axvline(wab, c='pink',linewidth=9)



 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
'''                     F O R C E    L O O P                         '''
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
Forder = sorted(range(len(Jsum)), key=lambda kk: Fplot[kk])
# The Index list that gives Increasing FREQUENCIES ^
Dorder = Forder[::-1] # reversed order
D = []                # Detuning
Rabi = []             # Rabi Freq
II = []               # II is the new IoIs = Intensity / Saturation Intensity 

L = len(Forder)
for i in range(L):
    dd = wab - Fplot[Dorder[i]]
    RR = dip*Jsum[Forder[i]]*E0/hbar
    D.append(dd)
    Rabi.append(RR)

for i in range(L):
    ii = 2*Rabi[i]**2/G**2 
    II.append(ii)

w_=[]
c1 = []
for i in range(L):
    ww = Fplot[Forder[i]]
    ccaa = 1+II[i]+4*D[i]**2/G**2
    w_.append(ww)
    c1.append(ccaa)

w = wab - AOM                   # Average Freq of colliding photon
Lambda = 2*np.pi*c/w            # Avg Wavelength
k = 2*np.pi/Lambda              # Average wavenumber of a momentum transfering photon

def dv(t,z,v):        
    "  The 'complete' Force Equation for a 7 freq 1 dimensional slower inc. magnetic field  "#    F O R C E   F O R C E   F O R C E   F O R C E   F O R C E   F O R C E
    O = []
    for i in range(L):
        ooo = w_[i]/(2*np.pi*c)#-muB*MagLeak(z, z0, Curr)/hbar
        O.append(ooo)
    c2,c3 = [],[] 
    for i in range(L):
        cc22 = O[i]*8*D[i]/G**2
        cc33 = 4*O[i]**2/G**2 
        c2.append(cc22)
        c3.append(cc33)
    Rho = []
    for i in range(L):
        rrr = II[i]/(c1[i]-c2[i]*v+c3[i]*v**2)
        Rho.append(rrr)
    rhoaa = np.sum(Rho)
    return -rhoaa*hbar*k*G/M/2

def dv4y(t,z,v):
    return 0

def dz(t,z,v):
    return v  
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
'''                I N I T I A L   C O N D I T I O N S                  '''
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def y_vGauss(n):
#    T = 273+20
#    g = np.random.normal(loc=0.0, scale=sc.k*T/M, size=n)
    g = np.random.normal(loc=0.0, scale=0.05, size=n)
#   g = np.linspace(-1,1, n)
    return g
def MBrand(n):
    #x = np.random.rand(n)
    x = np.linspace(0.0001,1,n)
    kb = sc.Boltzmann          # Boltzmann Constant
    u = sc.proton_mass         # Proton Mass
    M = 87*u                   # Mass of 87Rb
    a=abs(np.sqrt((kb*370)/(M)))
    X = sts.maxwell.isf(x, scale=a)
    return X


def Jgen(T,n,M):
    ''' Maxwell-Boltzmann Distrb. '''
    T = 273+20
    MB=sts.maxwell
    kb=1.38e-23
    a=abs(np.sqrt((kb*T)/(M)))
    vt=MB.rvs(loc=0, scale=a, size=n, random_state=None)    
    return vt

Vn = 10
Vx = 260
def vgen(n):
    ''' Linear Velocity Distrb. '''
    lin = np.linspace(Vn,Vx,n)
    return lin
    
def zgen(n):
    ''' Linear Coordinate Distrb. '''
    lin = np.linspace(0,0,n)
    return lin


def yrand(n):
    ran = np.random.random(n)
    return (ran-0.5)*2*Col_Gap

v_ = np.mean(vgen(Natoms))
Ry = ((ymot-Col_Gap/2)*v_)/(z0**2+(ymot-Col_Gap/2)**2)**0.5
def ygen(n):
    ''' Linear Coordinate Distrb. '''
    
    #lin = np.linspace(-Col_Gap/2,Col_Gap/2,n)
    lin = np.linspace(0,0,n)
    return lin

def vyrand(n):
    ran = np.random.random(n)
    it = (ran-0.5)*2*Ry
    return it
'''
def MagLeak(z, z0, Curr): 
    #Mag Field from AntiHlmHltz coils (of center z0 [ >0 ]) that leaks into our slower
    x  = s/2
    ZZ = -z+z0
    zz = -ZZ
    A,B = ZZ/aa, x/aa
    Q = B**2+(1+A)**2
    k = (4*A/Q)**0.5
    B0 = Curr*sc.mu_0/(2*aa)
    K = scp.ellipk(k**2)
    E = scp.ellipe(k**2)    
    Br = 2*B0*(x/ZZ)/(np.pi*Q**0.5)*(E*(1+A**2+B**2)/(Q-4*A)-K)
    Bro = np.nan_to_num(Br)
    #
    A_ = zz/aa
    Q_ = B**2+(1+A_)**2
    k_ = (4*A_/Q_)**0.5
    K_ = scp.ellipk(k_**2)
    E_ = scp.ellipe(k_**2)    
    Br_ = -2*B0*(x/zz)/(np.pi*Q_**0.5)*(E_*(1+A_**2+B**2)/(Q_-4*A_)-K_)
    Br_o = np.nan_to_num(Br_)
    return  Br_o + Bro
'''
def RK4step(ti,zi,vi,h,dv,dz):            
    k11=dz(ti,zi,vi)
    k21=dv(ti,zi,vi)    
    k12=dz(ti+h/2,zi +(h/2)*k11,vi +(h/2)*k21)
    k22=dv(ti+h/2,zi +(h/2)*k11,vi +(h/2)*k21)    
    k13=dz(ti+h/2,zi +(h/2)*k12,vi +(h/2)*k22)
    k23=dv(ti+h/2,zi +(h/2)*k12,vi +(h/2)*k22)    
    k14=dz(ti+h,zi +(h)*k13,vi +(h)*k23)
    k24=dv(ti+h,zi +(h)*k13,vi +(h)*k23)    
    z1=zi+(h/6.0)*(k11+2.0*k12+2.0*k13+k14)  
    v1=vi+(h/6.0)*(k21+2.0*k22+2.0*k23+k24)               
    zi = z1
    vi = v1    
    return zi,vi

""" Physical & Atomic Constants """
kb = sc.Boltzmann          # Boltzmann Constant
mu0 = sc.mu_0              # Vacc Permtivity
muB = 9.2740099*10**-24    # Borh Magnetron
                           # Proton Mass
                           # hbar
dip = 3.485e-29            # dipole moment
                           # speed of light
pi = np.pi                 # pi
                           # Mass of 87Rb
                           # Freq of transation
#G = 38.11e6               # Gamma / Rate of SpE 
                           # dipole moment

"""creation of our array of velocities"""
#vlin=vgen(Natoms)
zlin=zgen(Natoms)
y_vlin=vyrand(Natoms)
#ylin=yrand(Natoms)
vlin=MBrand(Natoms)
ylin=zgen(Natoms)
#y_vlin=y_vGauss(Natoms)
zs,vs,ts=[],[],[]
ys,yvs=[],[]
print('Trajectorizing...')
"""this loop goes through all the atoms we've got and applies the force dv to them for a number of steps, Nsteps"""
for j in range(Natoms):
    vi = vlin[j] 
    zi = zlin[j]
    yvi= y_vlin[j]
    yi = ylin[j]
    for i in range(Nsteps):
        ti=h*i 
        zs.append(zi)
        vs.append(vi)   
        ts.append(ti)     
        ys.append(yi)
        yvs.append(yvi)
        z1,v1=RK4step(ti,zi,vi,h,dv,dz)
        y1,yv1=RK4step(ti,yi,yvi,h,dv4y,dz)
        yvi = yv1
        yi = y1
        zi = z1
        vi = v1
        
Y = np.reshape(ys, (Natoms,Nsteps))
V = np.reshape(vs, (Natoms,Nsteps))
Z = np.reshape(zs, (Natoms,Nsteps))
tt = np.array(ts)
thet = np.split(tt, Natoms)[1]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
'''        C A P T U R E    D E T E C T I O N          '''
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
z_ , z__ = z0 - 0.01, z0 + 0.01
y_ = 0.01
capV, capv = 50,15
n_ = []
print('Accumulating...')
for j in range(Natoms):
    for i in range(Nsteps):
        if (z_ < Z[j][i] < z__ and abs(Y[j][i]) < y_ and abs(V[j][i]) < capv):
            #Y[j],Z[j] = np.linspace(0,0,Nsteps), np.linspace(0,0.0001,Nsteps)
            nnn = 2
            n_ = np.append(n_, nnn)
        else:
            nnn = 0
            n_ = np.append(n_, nnn)
        
N = np.reshape(n_, (Natoms, Nsteps))
N0 = 0
for j in range(Natoms):
    for i in range(Nsteps):
        if N[j][i] == 2:
            N0 += 1
            break
print('Doing the rest...')
print('#Captured = ',N0)

V_ = np.mean(V[:,0])
Vmp,Vrms = V_*(3*pi/8)**0.5, V_*(2*pi/8)**0.5
print(V_)


Cslow = 'blue'
Cslowm = 'green'
Cfastm = 'orange'
Cfast = 'red'
vVv1 = min(V[:,0], key=lambda g:abs(g-Vmp))
vVv1v = int(np.where(vVv1 == V[:,0])[0])
vVv2 = min(V[:,0], key=lambda g:abs(g-V_))
vVv2v = int(np.where(vVv2 == V[:,0])[0])
vVv3 = min(V[:,0], key=lambda g:abs(g-Vrms))
vVv3v = int(np.where(vVv3 == V[:,0])[0])

#print(int(vVv1v),vVv2v,vVv3v)
col=[]


for i in range(vVv1v):
    col.append(Cfast)
for j in range(vVv1v,vVv2v):
    col.append(Cfastm)
for k in range(vVv2v,vVv3v):
    col.append(Cslowm)
for l in range(vVv3v,Natoms):
    col.append(Cslow)
print(len(col))
for i in range(Natoms):
    'A plot for each of the Natoms particles'
    th = 0.4
    ax1.plot(Z[i],Y[i],linewidth=th, color = col[i])
    ax2.plot(Z[i],V[i],linewidth=th, color = col[i])
print(len(V[4]))

ax1.axhspan(-0.01,0.01, alpha=0.05, color='green')
ax1.axvspan(z0-0.01,z0+0.01, alpha=0.05, color='purple')
ax1.axvline(x = z0 - aa, color = 'k', linestyle='dotted')
ax1.axvline(x = z0, color = 'k', linestyle='dashed')
ax1.axvline(x = z0-0.01, color = 'k',linewidth=1.3)
ax1.axvline(x = z0+0.01, color = 'k',linewidth=1.3)
ax1.axhline(y = Col_Gap/2, color = 'k',linewidth=1.3)
ax1.axhline(y = -Col_Gap/2, color = 'k',linewidth=1.3)
ax2.set_ylim(top=1.5*ymot, bottom=-1.5*ymot)

ax2.axvspan(z0-0.01,z0+0.01, alpha=0.05, color='purple')
#ax2.axhspan(-capV,capV, alpha=0.05, color='b')
ax2.axhspan(-capv,capv, alpha=0.05, color='red')
ax2.axvline(x = z0 - aa, color = 'k', linestyle='dotted')
ax2.axvline(x = z0, color = 'k', linestyle='dashed')
ax2.axvline(x = z0-0.01, color = 'k',linewidth=1.3)
ax2.axvline(x = z0+0.01, color = 'k',linewidth=1.3)
ax2.axhline(y = capv, color = 'k',linewidth=1.3)
ax2.axhline(y = -capv, color = 'k',linewidth=1.3)

ax2.set_ylim(top=1.7*Vmp, bottom=-20)

#ax1.subplots_adjust(hspace=0) # Makes the plots that share the 
#                             # same x axis on top of each other
ax1.set_ylabel("y coordinate / m", size = 15)
ax2.set_ylabel("Speed / ms`'", size = 11)
ax1.set_title('Multi-Frequency: Sim-ya-lator Simulator', size=31)
ax5.set_xlabel('Distance / m (0 = atom source)', size = 15)
ax5.set_ylabel('Spatial Desity @ v=0',size=11)
#ax1.set_yticks(np.arange(-0.002, 0.002, step=0.0005))
from datetime import date
today = date.today()
d4 = today.strftime("%d-%b-%Y")

stop = timeit.default_timer()
IrE = c*8.85e-12/2*E0**2/10000   # Intensity (This /10000 makes it W/cm^2)

print('Run Time =',round(stop - start, 3),'sec')
print(np.sum(Jsum))
STR = '{}mW/cm2 [{}] \n{}MHz =RF \n{}MHz =AOM \n{} =Beta [{}]\n{}% Captured\n[{}/{}]\n    {} {}s' 
iSTR = STR.format(round(IrE,3)*1000, E0, RF/1000000, (AOM)/1000000, Beta,  round(np.sum(Jsum),4), round(N0/Natoms,3)*100, N0,Natoms,d4,round(stop - start, 3))
ax2.text(z0*0.4, 230, iSTR,fontsize=18, bbox = dict(boxstyle='round', fc=(0.0,0.0,0.0), alpha=0.2))
 # # # # # # # # # # # # # # # # # # # 
'''    Velocity class tracking      '''
v_ = np.linspace(0,1000,1000)
def dv_(v_,F,D):
    Lambda=2*pi*c/F
    k = 2*pi/Lambda
    'Incremental Acceleration'
    O = F/(2*pi*c)
    c1 =1+II[i]+4*D**2/G**2
    c2 = O*8/G**2*D
    c3 = 4*O**2/G**2  
    rhoaa =  -II[i]**2/(c1+c2*v_+c3*v_**2) + II[i]**2/(c1-c2*v_+c3*v_**2)   
    return rhoaa*hbar*k*G/M


S,Ss,Vv=[],[],[] 
for i in range(L):#reshaped Force(v) 
    S.append(dv_(v_,Fplot[i],D[i]))
S = np.reshape(S,(len(Fplot),len(v_)))
for i in range(L):#max indexes
    Sim = np.argmax(S[i])
    Ss = np.append(Ss,Sim)
for i in range(L):#indexing velocity
    Vv.append(v_[int(Ss[i])])
for i in range(L):
    ax2.axhline(Vv[i],c='brown', linewidth=1.7, alpha=0.34)
stop = timeit.default_timer()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
'''                V=0 density curve                  '''
 # # # # # # # # # # # # # # # # # # # # # # # # # # # #
'''
BinV=[]
V0s = []
for j in range(Natoms):    
        inx = (np.abs(V[j]-0)).argmin()
       # print(inx,Z[j][inx])
        V0s.append(Z[j][inx])
#print(V0s)
        
BinFactor = 4
Bins = np.linspace(0,2.5*z0,int(Natoms/BinFactor))  # Bins # Bins # Bins # Bins # Bins # Bins # Bins # Bins # Bins #
print(Bins)
His =np.histogram(V0s,bins=Bins)[0]
print('His',His,len(His))

Bplot = []
for i in range(len(Bins)-1):
    Bplot.append((Bins[i]+Bins[1+i])/2)
#print('b',Bplot)
ax5.plot(Bplot,His,linewidth=4,color='turquoise')
ax5.fill_between(Bplot,His,0,alpha=0.5,color='turquoise')
        
#  Hpeak[0,1]=value, index
Hpeak = [max(His), int(round( np.median(np.where(His==max(His))) ))]

Lmin = max(np.where(His[:Hpeak[1]] == min(His[:Hpeak[1]] ) )[0])

Lmax = max(np.where(His[Hpeak[1]:] == min(His[Hpeak[1]:] ) )[0])

#vLmin,vLmax = BinFactor*  IS IT POSSIBLE TO CONVERT INTO ORGINAL VELOCITY = not needed right now
FWi = Lmax-Lmin

Bot = max(His[Lmax],His[Lmin])

#print(Bot)
HM = Bot + (Hpeak[0]+Bot)/2


lHM = np.abs(His[:Hpeak[1]]-HM).argmin()
rHM = np.abs(His[Hpeak[1]:]-HM).argmin()+Hpeak[1]
print(lHM,rHM)
#print(lHM,rHM)
Skew = -1*(Bplot[Hpeak[1]]-Bplot[lHM]-Bplot[rHM]+Bplot[Hpeak[1]])
print('Skew =',Skew,'  +=MB')
FWHM = Bplot[rHM]-Bplot[lHM]
print('FWHM =', FWHM)

ax5.axhline(y=Hpeak[0],c='k',linewidth=0.4)
ax5.axhline(y=Bot,c='dimgrey',linewidth=0.4)
ax5.axvline(x=Bplot[lHM],c='k',linewidth=0.4)
ax5.axvline(x=Bplot[rHM],c='k',linewidth=0.4)
Strd = 'FWHM = {}mm\n÷Cloud Diam = {}\n[Skew = {}]'
StrD = Strd.format(round(FWHM,2)*1000,round(FWHM/ymot,4),round(Skew,4)) 
ax5.text(0.01,2,StrD,fontsize=16)
ax2.set_xlim(left=0, right=3*z0)
#print('# Particles = {}'.format(Natoms))
#print('Beam Intensity = {}W/cm^2'.format(round(IrE, 3)))
#print('Run Time =',round(stop - start, 3),'sec')
print('(',E0,AOM,RF,Beta,')',d4)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#    APPENDING DISTDENS TO A FILE TO BE USED LATER    B First    

import winsound
frequency = 2300  # Set Frequency To 2500 Hertz
duration = 1300  # Set Duration To 1000 ms == 1 second
winsound.Beep(frequency, duration)

import csv

with open('EOM_Beta_data{}.csv'.format(int(Beta*10)),'w',newline='') as f:#
    Wri = csv.writer(f)
    Wri.writerow(['Bplot,His,w/ Beta=','{}'.format(Beta)])
    for i in range(len(Bplot)):
        Wri.writerow(['{}'.format(Bplot[i]),'{}'.format(His[i]) ])
'''
plt.show()