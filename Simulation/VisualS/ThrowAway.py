import numpy as np 
import scipy.stats as sts
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import scipy.constants as sc
import scipy.special as scp
import timeit
start = timeit.default_timer()
plt.close('all')
fig = plt.figure()
ax1 = plt.subplot2grid((4,4), (0,1), colspan=3, rowspan=2)
ax2 = plt.subplot2grid((4,4), (2,1), colspan=3,sharex=ax1, rowspan=2)
ax3 = plt.subplot2grid((4,4), (0,0), rowspan=3)#, rowspan=2)
ax4 = plt.subplot2grid((4,4), (3,0), rowspan=1)

E0 = 100
Beta = 3.0 # currently not sure if this is a fn Beta(RF) only & how it is connected to it ;  if it is then need only change EOM(args*)

RF = 0.5e6
AOM = (100e6) #    AOM s detuning (effective detuning, will be doubled)

aa = 0.15       # Coil Radius
s = 0.11        # Coil Separation
Curr = 0.6805   # Current (16G/cm @ 0.6805)
z0 = 0.24        # Position of MOT centre

h = 0.0001     #Step Size
Natoms = 50
Nsteps = 80

Col_zPos = 0.008
Col_Gap = 0.006
#  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  #
'''                E O M   &  F r e q  D u b    t i m e                   '''
#  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  #
# E0 is central, E1p E1m is +- 1st Sideband Pair    &    same for F

G = 38.11e6                                  
#wab = 2*np.pi*384.23e12    # Freq of transation
wab = 400e12    # Freq of transation
hbar = sc.hbar             # hbar
dip = 3.485e-29            # dipole moment
c = sc.c
u = sc.proton_mass         # Proton Mass
M = 87*u                   # Mass of 87Rb
beta = np.linspace(0,4,80)

#  Params ^^

def EOM(RF, Beta): # [][1=Freq  0=E],AOM
    'Convention: [0,1,2,3,4] = Center High1stSideBand Low1stSideBand High2ndSideBand... '
    Ei0 = scp.jv(0,Beta)*E0
    Ei1 = scp.jv(1,Beta)*E0
    Ei2 = scp.jv(2,Beta)*E0
    F0 = wab/2 
    Fp1, Fm1 = wab/2/+RF, wab/2-RF
    Fp2, Fm2 = wab/2+2*RF, wab/2-2*RF
    return [Ei0,F0],[Ei1,Fp1],[-Ei1,Fm1],[Ei2,Fp2],[-Ei2,Fm2]
    #So EOM[band][E,F]
print('+++++++++', EOM(RF,Beta))
ax4.plot(beta, abs(EOM(RF, beta)[0][0]/E0),linewidth=6,c='darkviolet')
ax4.plot(beta, abs(EOM(RF, beta)[2][0]/E0),linewidth=4,c='mediumorchid')
ax4.plot(beta, abs(EOM(RF, beta)[4][0]/E0),linewidth=2,c='violet')
ax4.axvline(Beta,c='k', linestyle='dashed')
ax4.axhline(abs(EOM(RF, Beta)[0][0]/E0),linewidth=6,c='darkviolet', linestyle='dotted')
ax4.axhline(abs(EOM(RF, Beta)[2][0]/E0),linewidth=4,c='mediumorchid', linestyle='dotted')
ax4.axhline(abs(EOM(RF, Beta)[4][0]/E0),linewidth=2,c='violet', linestyle='dotted')

# Freq and Efield together
A0  = [EOM(RF, Beta)[0][1],EOM(RF, Beta)[0][0]]
A1p = [EOM(RF, Beta)[1][1],EOM(RF, Beta)[1][0]]
A1m = [EOM(RF, Beta)[2][1],EOM(RF, Beta)[2][0]]
A2p = [EOM(RF, Beta)[3][1],EOM(RF, Beta)[3][0]]
A2m = [EOM(RF, Beta)[4][1],EOM(RF, Beta)[4][0]]
#print('A0',A0)
#print('A1p',A1p)
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
    ax3.scatter(EE[j][0]-AOM,EE[j][1],s=100,alpha=0.3) 
    ax3.scatter(AA[j][0]-AOM,AA[j][1],s=100,alpha=0.3)
    ax3.scatter(BB[j][0]-AOM,BB[j][1],s=100,alpha=0.3)
    ax3.scatter(CC[j][0]-AOM,CC[j][1],s=100,alpha=0.3)
    ax3.scatter(DD[j][0]-AOM,DD[j][1],s=100,alpha=0.3)
#print('-_-_-_-',DD[j][0],DD[j][1])

XX = []
Dup = []

BBig = np.array([a0,a1p,a1m,a2p,a2m,b0,b1p,b1m,b2p,b2m,c0,c1p,c1m,c2p,c2m,d0,d1p,d1m,d2p,d2m,e0,e1p,e1m,e2p,e2m])
Flist = BBig[:,0]
#print(' s s s -',len(Flist),Flist)
Fplot_ = []
print('@@@', BBig)
for k in range(len(BBig)):
    '''test for dups - get index symm - choose 1st'''
    x = np.where(BBig[:,0] == Flist[k])[0]
    xx = x[0]
    XX.append(x[0])
XX_ = list(dict.fromkeys(XX)) # Dup Prepended list
print('@b',XX_)
for i in range(len(XX_)):
    '''Takes the 1st index & gives freq'''
    Fplot_.append(Flist[XX_[i]])
    Fplot = np.subtract(Fplot_,AOM)
print('@2b',Fplot)
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
print('@3b',X0_U)
Yy=[]
for j in range(len(Fplot)):
    for k in range(len(X0_U[j])):
        y = BBig[X0_U[j][k]][1]
      #  print('-',y)
        Yy.append([y])
    Y = np.sum(Yy)
 #   print('Y',Y)
    Jsum.append(Y)
    Yy=[]
    
print("@4b",len(Jsum), Jsum)
print('bFreqs =',len(Fplot),Fplot)

ax3.scatter(Fplot,Jsum,s=800,c='y', alpha=0.4)
ax3.set_xlabel('MHz',size=10)
ax3.set_ylabel('Electric Field / E0',size=15)
ax4.set_xlabel(r'$ \beta $',size=22)
ax4.set_ylabel('Bessel Value',size=13)
ax3.grid(which='major')
#ax3.set_xlim(left=wab-2*AOM-RF, right=wab-2*AOM+RF) fucked


'''         E N D   of   E O M   &   F D         '''
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#Db,Dc,Dd,De,Df,Dg = Da+d, Da+2*d, Da+3*d, Da+4*d, Da+5*d, Da+6*d
xxxxx = np.linspace(0.000001, 1, 100)

#    I N T E N S I T Y S     --    we will do in order of increasing Intensities! J will follow 
Forder = sorted(range(len(Jsum)), key=lambda k: Fplot[k])
print(Forder)
# Da b c d e = mid right left right left 



# # # # # # # # # # # # # # # # # # # # # # # # # # # #
'''  N E W   L O O P  T I M E '''
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
L = len(Forder)

D = []
for i in range(L):
    dd = wab - Fplot[Forder[i]]
    D.append(dd)
print('D',D)
#Da = wab - Fplot[Forder[0]]

Rabi = []
for i in range(L):
    RR = dip*Jsum[Forder[i]]/hbar
    Rabi.append(RR)
print('RAbi',Rabi)
#Rabi7 = dip*Jsum[Forder[7]]/hbar   

                         #      II is the new IoIs = Intensity / Saturation Intensity 
II = []
for i in range(L):
    ii = 2*Rabi[i]**2/G**2 
    II.append(ii)
print('II',II)
print('IoIs-------::',II)
print('Dd--------::',D)
#II7 = 2*Rabi7**2/G**2 
# Intensity (This /10000 makes it W/cm^2)
#IrE = c*8.85e-12/2*II7 = 2*Rabi7**2/G**2**2/10000
w = (wab) - AOM                   # Average Freq of colliding photon
Lambda = 2*np.pi*c/w                  # Avg Wavelength
k = 2*np.pi/Lambda                  # Average wavenumber of a momentum transfering photon

w_=[]
for i in range(L):
    ww = wab - D[i]
    w_.append(ww)
print(w_)

c1 = []
for i in range(L):
    ccaa = 1+II[i]+4*D[i]**2/G**2
    c1.append(ccaa)

def dv(t,z,v):        
    "  The 'complete' Force Equation for a 7 freq 1 dimensional slower inc. magnetic field  "#    F O R C E   F O R C E   F O R C E   F O R C E   F O R C E   F O R C E

    #w_e = wab - De 

    O = []
    for i in range(L):
        ooo = w_[i]/(2*np.pi*c)#-muB*MagLeak(z, z0, Curr)/hbar
        O.append(ooo)
    #Og = w_g/(2*np.pi*c)      #-muB*MagLeak(z, z0, Curr)/hbar          #         #         #        #       #   #  #  # # ### ###### ########################

    c2 = []
    for i in range(L):
        cc22 = O[i]*8*D[i]**2/G**2
        c2.append(cc22)
    #c1g = 1+II6+4*Dg**2/G**2 
    #c2g = Og*8*Dg/G**2 

    c3 = [] 
    for i in range(L):
        cc33 = 4*O[i]**2/G**2 
        c3.append(cc33)
    #c3g = 4*Og**2/G**2 

    Rho = []
    for i in range(L):
        rrr = -II[i]/(2*(c1[i]-c2[i]*v+c3[i]*v**2))
        Rho.append(rrr)
    rhoaa = np.sum(Rho)
    #rhoaa =  -(II0/2)*(1/(c1a-c2a*v+c3a*v**2)+1/(c1b-c2b*v+c3b*v**2)+1/(c1c-c2c*v+c3c*v**2)+1/(c1d-c2d*v+c3d*v**2)+1/(c1e-c2e*v+c3e*v**2)+1/(c1f-c2f*v+c3f*v**2)+1/(c1g-c2g*v+c3g*v**2))   
    return rhoaa*hbar*k*G/M
print('F', dv(2,2,200))

def y_vgen(n):
    ''' Linear Velocity distribution in Y / spread of particles '''
    lin = np.linspace(0.5,35,n)
    return lin

def y_vGauss(n):
#    T = 273+20
#    g = np.random.normal(loc=0.0, scale=sc.k*T/M, size=n)
    g = np.random.normal(loc=0.0, scale=1.5, size=n)
#   g = np.linspace(-1,1, n)
    return g

def ygen(n):
    ''' Linear Coordinate Distrb. '''
    lin = np.linspace(0,0,n)
    return lin

def Jgen(T,n,M):
    ''' Maxwell-Boltzmann Distrb. '''
    T = 273+20
    MB=sts.maxwell
    kb=1.38e-23
    a=abs(np.sqrt((kb*T)/(M)))
    vt=MB.rvs(loc=0, scale=a, size=n, random_state=None)    
    return vt

Vn = 70
Vx = 260
def vgen(n):
    ''' Linear Velocity Distrb. '''
    lin = np.linspace(Vn,Vx,n)
    return lin
    
def zgen(n):
    ''' Linear Coordinate Distrb. '''
    lin = np.linspace(0,0,n)
    return lin
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

''' Variable Dance '''

def dv4y(t,z,v):
    return 0

def dz(t,z,v):
    return v  
    


"""creation of our array of velocities"""
#vlin=Jgen(T,Natoms,M)
vlin=vgen(Natoms)
zlin=zgen(Natoms)
#y_vlin=y_vgen(Natoms)
y_vlin=y_vGauss(Natoms)
ylin=ygen(Natoms)
zs,vs,ts=[],[],[]
ys,yvs=[],[]
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
        
        z1=RK4step(ti,zi,vi,h,dv,dz)[0]
        v1=RK4step(ti,zi,vi,h,dv,dz)[1]
        y1=RK4step(ti,yi,yvi,h,dv4y,dz)[0]
        yv1=RK4step(ti,yi,yvi,h,dv4y,dz)[1]
        
        yvi = yv1
        yi = y1
        zi = z1
        vi = v1
        
Y = np.reshape(ys, (Natoms,Nsteps))
V = np.reshape(vs, (Natoms,Nsteps))
Z = np.reshape(zs, (Natoms,Nsteps))
tt = np.array(ts)
thet = np.split(tt, Natoms)[1]

Top, Thicc = 0.002, 0.003
ax1.bar(Col_zPos, Top, Thicc, bottom= Col_Gap/2, color='k')
ax1.bar(Col_zPos,-Top, Thicc, bottom=-Col_Gap/2, color='k')
#print(Y, Y.shape)
    
'''    Collimation Collision Detection     '''

nn=0

    
'''    Capture Detection     '''

z_ , z__ = z0 - 0.01, z0 + 0.01
y_ = 0.01
capV, capv = 50,15
n_ = []

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

#print(n_)
#rint(N)
N0 = 0
for j in range(Natoms):
    NN = False
    for i in range(Nsteps):
        if N[j][i] == 2:
            NN = True
    if NN == True:
        N0 += 1
print(N0)


for i in range(Natoms):
    'A plot for each of the Natoms particles'
    th = 0.5
    col = (0.1, float(i/(Natoms+1)+0.0001), 1-float(i/(Natoms+5)+0.0001))
    ax1.plot(Z[i],Y[i],linewidth=th, color = col)
    ax2.plot(Z[i],V[i],linewidth=th, color = col)

ax1.axhspan(-0.01,0.01, alpha=0.05, color='green')
ax1.axvspan(z0-0.01,z0+0.01, alpha=0.05, color='purple')
ax1.axvline(x = z0 - aa, color = 'k', linestyle='dotted')
ax1.axvline(x = z0, color = 'k', linestyle='dashed')
ax1.axvline(x = z0-0.01, color = 'r')
ax1.axvline(x = z0+0.01, color = 'r')
ax1.axhline(y = Col_Gap/2, color = 'r')
ax1.axhline(y = -Col_Gap/2, color = 'r')

ax2.axvspan(z0-0.01,z0+0.01, alpha=0.05, color='purple')
ax2.axhspan(-capV,capV, alpha=0.05, color='b')
ax2.axhspan(-capv,capv, alpha=0.05, color='red')
ax2.axvline(x = z0 - aa, color = 'k', linestyle='dotted')
ax2.axvline(x = z0, color = 'k', linestyle='dashed')
ax2.axvline(x = z0-0.01, color = 'r')
ax2.axvline(x = z0+0.01, color = 'r')
ax2.axhline(y = capv, color = 'r')
ax2.axhline(y = -capv, color = 'r')
#ax2.set_xlim(left=-0.009, right=z0+0.4*aa)
#ax2.set_ylim(top=Vx+5, bottom=-40)



#ax1.subplots_adjust(hspace=0) # Makes the plots that share the 
#                             # same x axis on top of each other
ax1.set_ylabel("y coordinate / m", size = 15)
ax2.set_ylabel("Speed / ms`'", size = 11)
#ax3.set_title('Multi-Frequency Slowing Simulation: $\it{7}$ $\it{Frequencies}$, $\it{MOT}$ $\it{Magnetic}$ $\it{Field}$', size=17)
ax1.set_title('Multi-Frequency: Sim-ya-lator Simulator', size=31)
ax2.set_xlabel('Distance / m (0 = atom source)', size = 15)
#ax1.set_yticks(np.arange(-0.002, 0.002, step=0.0005))
from datetime import date
today = date.today()
d4 = today.strftime("%d-%b-%Y")

STR = '{}mW/cm2 =Intns \n{}MHz = RF \n{}= AOM/2 \n{}=Beta \n{}' 
iSTR = STR.format('-', RF/1000000, (AOM*2)/1000000, Beta,d4)
ax2.text(z0+1, 150, iSTR,fontsize=18, bbox = dict(boxstyle='round', fc=(0.19,0.68,0.5), alpha=0.4))


stop = timeit.default_timer()
#print('# Particles = {}'.format(Natoms))
#print('Beam Intensity = {}W/cm^2'.format(round(IrE, 3)))
#print('Run Time =',round(stop - start, 3),'sec')

plt.show()