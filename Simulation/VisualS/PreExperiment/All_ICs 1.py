import numpy as np 
import scipy.stats as sts
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import scipy.constants as sc
import scipy.special as scp
import timeit
start = timeit.default_timer()

G = 38.11e6     # See ln 96

E0 = 130
Da = 4.8*G    # Original detuning
d = 0.3*G    # Increment of Detuning

aa = 0.15       # Coil Radius
s = 0.11        # Coil Separation
Curr = 0.6805   # Current (16G/cm @ 0.6805)
z0 = 0.24        # Position of MOT centre

h = 0.0001     #Step Size
Natoms = 4
Nyvel = 7
Nycor = 10
Nsteps = 50

Col_Gap = 0.02
Col_zPos = 0
Db,Dc,Dd,De,Df,Dg = Da+d, Da+2*d, Da+3*d, Da+4*d, Da+5*d, Da+6*d
xx = np.linspace(0.000001, 1, 100)

def y_vgen(n):
    ''' Linear Velocity distribution in Y / spread of particles '''
    lin = np.linspace(-30,30,n)
    return lin
'''
def y_vGauss(n):
#    T = 273+20
#    g = np.random.normal(loc=0.0, scale=sc.k*T/M, size=n)
    g = np.random.normal(loc=0.0, scale=1.0, size=n)
#   g = np.linspace(-1,1, n)
    return g
'''
def ygen(n):
    ''' Linear Coordinate Distrb. '''
    #lin = np.linspace(-Col_Gap/2,Col_Gap/2,n)
    lin = np.linspace(-0.01,0.01,n)
    return lin
    
Vn = 50
Vx = 250
def vgen(n):
    ''' Linear Velocity Distrb. '''
    lin = np.linspace(Vn,Vx,n)
    return lin
    
def Jgen(T,n,M):
    ''' Maxwell-Boltzmann Distrb. '''
    T = 273+20
    MB=sts.maxwell
    kb=1.38e-23
    a=abs(np.sqrt((kb*T)/(M)))
    vt=MB.rvs(loc=0, scale=a, size=n, random_state=None)    
    return vt


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
    k12=dz(ti+h/2, zi + (h/2)*k11, vi + (h/2)*k21)
    k22=dv(ti+h/2, zi + (h/2)*k11, vi + (h/2)*k21)    
    k13=dz(ti+h/2, zi + (h/2)*k12, vi + (h/2)*k22)
    k23=dv(ti+h/2, zi + (h/2)*k12, vi + (h/2)*k22)    
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
u = sc.proton_mass         # Proton Mass
hbar = sc.hbar             # hbar
c = sc.c                   # speed of light
pi = np.pi                 # pi
M = 87*u                   # Mass of 87Rb
wab = 2*pi*384.23e12       # Freq of transation
#G = 38.11e6               # Gamma / Rate of SpE 
dip = 3.485e-29            # dipole moment

''' Variable Dance '''
Rabi = dip*E0/hbar               # Rabi Frequency 
IoIs = 2*Rabi**2/G**2            # Intensity / Saturation Intensity
IrE = c*8.85e-12/2*E0**2/10000   # Intensity (This /10000 makes it W/cm^2)
w = wab - Dd                     # Average Freq of colliding photon
Lambda = 2*pi*c/w                  # Avg Wavelength
k = 2*pi/Lambda                  # Average wavenumber of a momentum transfering photon

def dv4y(t,z,v):
    return 0

def dv(t,z,v):        
    " The 'complete' Force Equation for a 7 freq 1 dimensional slower inc. magnetic field "
    w_a = wab - Da
    w_b = wab - Db 
    w_c = wab - Dc 
    w_d = wab - Dd 
    w_e = wab - De 
    w_f = wab - Df 
    w_g = wab - Dg 
    Oa = w_a/(2*pi*c)#-muB*MagLeak(z, z0, Curr)/hbar
    Ob = w_b/(2*pi*c)#-muB*MagLeak(z, z0, Curr)/hbar
    Oc = w_c/(2*pi*c)#-muB*MagLeak(z, z0, Curr)/hbar
    Od = w_d/(2*pi*c)#-muB*MagLeak(z, z0, Curr)/hbar
    Oe = w_e/(2*pi*c)#-muB*MagLeak(z, z0, Curr)/hbar
    Of = w_f/(2*pi*c)#-muB*MagLeak(z, z0, Curr)/hbar
    Og = w_g/(2*pi*c)#-muB*MagLeak(z, z0, Curr)/hbar
    c1a = 1+IoIs+4*Da**2/G**2
    c2a = Oa*8*Da/G**2
    c1b = 1+IoIs+4*Db**2/G**2
    c2b = Ob*8*Db/G**2
    c1c = 1+IoIs+4*Dc**2/G**2
    c2c = Oc*8*Dc/G**2
    c1d = 1+IoIs+4*Dd**2/G**2
    c2d = Od*8*Dd/G**2
    c1e = 1+IoIs+4*De**2/G**2
    c2e = Oe*8*De/G**2
    c1f = 1+IoIs+4*Df**2/G**2
    c2f = Of*8*Df/G**2
    c1g = 1+IoIs+4*Dg**2/G**2
    c2g = Og*8*Dg/G**2
    c3a = 4*Oa**2/G**2 
    c3b = 4*Ob**2/G**2 
    c3c = 4*Oc**2/G**2 
    c3d = 4*Od**2/G**2 
    c3e = 4*Oe**2/G**2 
    c3f = 4*Of**2/G**2 
    c3g = 4*Og**2/G**2 
    rhoaa =  -(IoIs/2)*(1/(c1a-c2a*v+c3a*v**2)+1/(c1b-c2b*v+c3b*v**2)+1/(c1c-c2c*v+c3c*v**2)+1/(c1d-c2d*v+c3d*v**2)+1/(c1e-c2e*v+c3e*v**2)+1/(c1f-c2f*v+c3f*v**2)+1/(c1g-c2g*v+c3g*v**2))   
    return rhoaa*hbar*k*G/M

def dz(t,z,v):
    return v  
    
plt.close('all')
fig = plt.figure()
ax1 = plt.subplot2grid((2,1), (0,0))#, rowspan=2)
ax2 = plt.subplot2grid((2,1), (1,0), sharex=ax1)



"""creation of our array of velocities"""
#vlin=Jgen(T,Natoms,M)
vlin=vgen(Natoms)
zlin=zgen(Natoms)
#y_vlin=y_vgen(Natoms)
y_vlin=y_vgen(Nyvel)       # 20 velocities at each atom site
ylin=ygen(Nycor)
zs,vs,ts=[],[],[]
ys,yvs=[],[]

#YVy = np.column_stack((ylin,y_vlin))

def YVy(y_vlin, ylin, i, j):
    '''   i = v  j = y  ===>    [0] = v   [1] = y  '''
    x = y_vlin[i]
    x_ = ylin[j]
    return x,x_

def All_IC(y_vlin, ylin, vlin, j, k, l):
    x = y_vlin[j]
    x_ = ylin[k]
    x__ = vlin[l]
    return x,x_,x__
#yy = np.transpose(np.column_stack((ylin,ylin,ylin,ylin,ylin,ylin,ylin,ylin,ylin,ylin)))

#for i in range(Nvel):
#    print(YVy(y_vlin, ylin, i, 2)[0])

"""this loop goes through all the atoms we've got and applies the force dv to them for a number of steps, Nsteps"""
for l in range(Natoms):    
    for j in range(1,Nyvel):
        for k in range(Nycor):
            zi = zlin[l]
            yvi = All_IC(y_vlin, ylin, vlin, j, k, l)[0]
            yi = All_IC(y_vlin, ylin, vlin, j, k, l)[1]
            #yi = 0    NOTHING WRONG WITH THIS^
            vi = All_IC(y_vlin, ylin, vlin, j, k, l)[2]
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

'''
       Z   V e l o c i t y    O N L Y
'''
zs_,vs_,ts_=[],[],[]
for j in range(Natoms):
    vi_ = vlin[j] 
    zi_ = zlin[j]
    
    for i in range(Nsteps):
        ti_=h*i 
        zs_.append(zi_)
        vs_.append(vi_)   
        ts_.append(ti_)     
        z1_=RK4step(ti_,zi_,vi_,h,dv,dz)[0]
        v1_=RK4step(ti_,zi_,vi_,h,dv,dz)[1]
        zi_ = z1_
        vi_ = v1_
        print(vi_, zi_)
        
#Y = np.reshape(ys, (Natoms,Nsteps))
V_ = np.reshape(vs_, (Natoms,Nsteps))
Z_ = np.reshape(zs_, (Natoms,Nsteps))

print(len(ys))
print(All_IC(y_vlin, ylin, vlin, 2, 2, 2)[2])

Y = np.reshape(ys, (Natoms*Nycor*(Nyvel-1),Nsteps))
V = np.reshape(vs, (Natoms*(Nyvel-1)*Nycor,Nsteps))
Z = np.reshape(zs, (Natoms*(Nyvel-1)*Nycor,Nsteps))
#tt = np.array(ts)
#thet = np.split(tt, Natoms)[1]
print('  -  ',len(V[4]))

'''
Top, Thicc = 0.002, 0.003
ax1.bar(Col_zPos, Top, Thicc, bottom= Col_Gap/2, color='k')
ax1.bar(Col_zPos,-Top, Thicc, bottom=-Col_Gap/2, color='k')
'''
'''     Target Region  TOO BE     '''
'''
nn=0
for j in range(Natoms*10):
    for i in range(Nsteps):
        if (Z[j][i] < Col_zPos and abs(Y[j][i]) > Col_Gap/2):
            Y[j],Z[j] = np.linspace(0,0,Nsteps), np.linspace(0,0.0001,Nsteps)
            nn += 1

Z_o = Z_.append(Y.append(Y))
for i in range(Natoms*(Nyvel-1)*Nycor):
    'A plot for each of the Natoms particles'
    th = 0.2
    col = (0.3, 1-i/(Natoms*Nyvel*Nycor+2), 0.8)
    ax1.plot(Z_[i],Y[i], linewidth=th, color = col)
    #ax2.plot(Z_[i],V_[i], linewidth=th, color = col)
'''
for i in range(Natoms):
    th = 2
    col = (0.8, 0, 0.2)#-i/(Natoms+2)
    #ax1.plot(Z[i],Y[i], linewidth=th, color = col)
    ax1.plot(Z_[i],V_[i], linewidth=th, color = col) 
    ax2.plot(Z_[i],V_[i], linewidth=th, color = col) 

capV, capv = 50,15
ax2.axhspan(-capV,capV, alpha=0.05, color='b')
ax2.axhspan(-capv,capv, alpha=0.05, color='red')
ax2.axvline(x = z0 - aa, color = 'k', linestyle='dotted')
ax2.axvline(x = z0, color = 'k', linestyle='dashed')
ax2.set_xlim(left=-0.009, right=z0+0.4*aa)
ax2.set_ylim(top=Vx+5, bottom=-40)



fig.subplots_adjust(hspace=0) # Makes the plots that share the 
#                             # same x axis on top of each other
ax1.set_ylabel("y coordinate / m", size = 15)
ax2.set_ylabel("Speed / ms`'", size = 11)
#ax3.set_title('Multi-Frequency Slowing Simulation: $\it{7}$ $\it{Frequencies}$, $\it{MOT}$ $\it{Magnetic}$ $\it{Field}$', size=17)
ax1.set_title('Multi-Frequency: 7 Frequencies and COLLIMATION', size=17)
ax2.set_xlabel('Distance / m (0 = atom source)', size = 15)
#ax1.set_yticks(np.arange(-0.002, 0.002, step=0.0005))


'''Slope Finding'''  
'''  
LenLin,dd = 1000, 0.03
LinGrad = np.linspace(z0-dd, z0+dd, LenLin)
Bp = MagLeak(LinGrad,z0,Curr)[0]
Bm = MagLeak(LinGrad,z0,Curr)[LenLin-1]
ax4.plot(LinGrad, MagLeak(LinGrad, z0, Curr)*1000, color='cyan', lw=0.5)
Grad = abs(Bp-Bm)/(2*dd)
'''
from datetime import date
today = date.today()
d4 = today.strftime("%d-%b-%Y")
#print("d4 =", d4)
#ax3.legend(title='     {}\nIntensity = {}W/cm2\nDetuning = {} w/ Increment {}MHz\nE0 = {}    no. atoms = {} \nLength of Tube = {}cm\nMag Field Gradient = {}G/cm'.format(d4,round(IrE, 3),Da/1000000,d/1000000,E0,nj, round((z0-aa)*100,3),round(Grad*1000000,2), loc=2, prop={'size': 18}))
nn=0
textstr = '     {}\nIntensity = {}mW/cm2\nDetuning = {}\n w/ Increment {}MHz\nE0 = {}    Total Atoms = {} \n #/Vz = {} #/Vy = {} #/y = {}\nLength of Tube = {}cm\nMag Field Gradient = {}G/cm\n\n    Hole Diameter = {}mm\n    Hole Distance = {}cm\n    Pariticles lost = {}%'.format(d4,round(IrE*1000, 3),Da/1000000,d/1000000,E0,Natoms*Nyvel*Nycor,Natoms,Nyvel,Nycor, round((z0-aa)*100,3),'Na', Col_Gap*1000, Col_zPos*100, (round((nn/Natoms)*100,3) ))#round(Grad*1000000,2))

ax2.text(z0-0.06, 150, textstr, fontsize=10)

stop = timeit.default_timer()
print('# Particles = {}'.format(Natoms))
print('Beam Intensity = {}W/cm^2'.format(round(IrE, 3)))
print('Run Time =',round(stop - start, 3),'sec')

plt.show()