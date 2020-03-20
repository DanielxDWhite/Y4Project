"""
             MFreq   w/         Atom Counter   &   f(Col_)
"""
import numpy as np 
import scipy.stats as sts
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.special as scp
import timeit
start = timeit.default_timer()

ymot = 0.004      # Radius of Col Atoms

G = 38.11e6     # See ln 96

Xi_D = 8.6
Xi_d = 0.15

E0 = 280
Da = Xi_D*G    # Original detuning
d = Xi_d*G    # Increment of Detuning
print(Da,d)

aa = 0.15       # Coil Radius
s = 0.11        # Coil Separation
Curr = 0.6805   # Current (16G/cm @ 0.6805)
z0 = 0.24        # Position of MOT centre

h = 0.00002     #Step Size
Natoms = 150
Nsteps = 410

Col_zPos = 0
Col_Gap = 0.002
'''If ^ > ymot then simulation is unphysical'''

J = 10e20
Q = J*(Col_Gap/2)**2

print('E0, Da, d is {}, {}, {}'.format(E0,Xi_D,Xi_d))
Db,Dc,Dd,De,Df,Dg = Da+d, Da+2*d, Da+3*d, Da+4*d, Da+5*d, Da+6*d
xx = np.linspace(0.000001, 1, 100)




def y_vgen(n):
    ''' Linear Velocity distribution in Y / spread of particles '''
    lin = np.linspace(-vy,vy,n)
    return lin

def ygen(n):
    ''' Linear Coordinate Distrb. '''
    lin = np.linspace(-Col_Gap/2,Col_Gap/2,n)
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

def yrand(n):
    ran = np.random.random(n)
    return (ran-0.5)*Col_Gap


def MBrand(n):
    x = np.random.rand(n)
    kb = sc.Boltzmann          # Boltzmann Constant
    u = sc.proton_mass         # Proton Mass
    M = 87*u                   # Mass of 87Rb
    a=abs(np.sqrt((kb*300)/(M)))
    X = sts.maxwell.isf(x, scale=a)
    return X

#v_ = 0.5*(max(MBrand(Natoms))+min(MBrand(Natoms))) #Mean Velocity
v_ = np.mean(MBrand(Natoms))
vy = ((ymot-Col_Gap/2)*v_)/(z0**2+(ymot-Col_Gap/2)**2)**0.5

def Vyrand(n):
    ran = np.random.random(n)
    return (ran-0.5)*vy*2

#q = Q*(1-np.cos(np.arctan((vy/v_))))
#print('q , Q =',q,Q)


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
hbar = sc.hbar             # hbar
c = sc.c                   # speed of light
pi = np.pi                 # pi
u = sc.proton_mass         # Proton Mass
M = 87*u                   # Mass of 87Rb
wab = 2*pi*384.23e12       # Freq of transation
#G = 38.11e6               # Gamma / Rate of SpE 
dip = 3.485e-29            # dipole moment
q =1
''' Variable Dance '''
Rabi = dip*E0/hbar               # Rabi Frequency 
IoIs = 2*Rabi**2/G**2 
print('IoIs',IoIs)           # Intensity / Saturation Intensity
IrE = c*8.85e-12/2*E0**2/10000   # Intensity (This /10000 makes it W/cm^2)
w = wab - Dd                     # Average Freq of colliding photon
Lambda = 2*pi*c/w                  # Avg Wavelength
k = 2*pi/Lambda                  # Average wavenumber of a momentum transfering photon
print(k,'k')
print(Rabi)
print((wab-Da)**2/G**2)
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
    #   print(c1g,c2g,c3a)
    c3b = 4*Ob**2/G**2 
    c3c = 4*Oc**2/G**2 
    c3d = 4*Od**2/G**2 
    c3e = 4*Oe**2/G**2 
    c3f = 4*Of**2/G**2 
    c3g = 4*Og**2/G**2 
    rhoaa =  -(IoIs/2)*(1/(c1a-c2a*v+c3a*v**2)+1/(c1b-c2b*v+c3b*v**2)+1/(c1c-c2c*v+c3c*v**2)+1/(c1d-c2d*v+c3d*v**2)+1/(c1e-c2e*v+c3e*v**2)+1/(c1f-c2f*v+c3f*v**2)+1/(c1g-c2g*v+c3g*v**2))   
    return rhoaa*hbar*k*G/M
print('F', dv(2,2,100))
print('F', dv(2,2,300))
print('F', dv(2,2,500))
def dz(t,z,v):
    return v  
    
plt.close('all')
fig = plt.figure()
ax1 = plt.subplot2grid((2,1), (0,0))#, rowspan=2)
ax2 = plt.subplot2grid((2,1), (1,0),sharex=ax1)


"""   s"""

zlin=zgen(Natoms)
yran = yrand(Natoms)
Vyran = Vyrand(Natoms)
vran = MBrand(Natoms)
print('begining Loop 1')
zs,vs,ts=[],[],[]
ys,yvs=[],[]
"""this loop goes through all the atoms we've got and applies the force dv to them for a number of steps, Nsteps"""
for j in range(Natoms):
    vi = vran[j] 
    zi = zlin[j]
    yvi= Vyran[j]
    yi = yran[j]
    for i in range(Nsteps):
        ti=h*i 
        zs.append(zi)
        vs.append(vi)   
        ts.append(ti)     
        ys.append(yi)
        yvs.append(yvi)
        z1,v1 = RK4step(ti,zi,vi,h,dv,dz)
        y1,yv1 = RK4step(ti,yi,yvi,h,dv4y,dz)
        yvi = yv1
        yi = y1
        zi = z1
        vi = v1
print('done Loop 1')
Y = np.reshape(ys, (Natoms,Nsteps))
V = np.reshape(vs, (Natoms,Nsteps))
Z = np.reshape(zs, (Natoms,Nsteps))
tt = np.array(ts)
thet = np.split(tt, Natoms)[1]

#Top, Thicc = 0.002, 0.003
#ax1.bar(Col_zPos, Top, Thicc, bottom= Col_Gap/2, color='k')
#ax1.bar(Col_zPos,-Top, Thicc, bottom=-Col_Gap/2, color='k')
#print(Y, Y.shape)
#ax2.subplots_adjust(hspace=0)
z_ , z__ = z0 - ymot, z0 + ymot
y_ = 0.01
capV, capv = 50,15
n_ = []
print('begin Loop 2')
for j in range(Natoms):
    for i in range(Nsteps):
        if (z_ < Z[j][i] < z__ and abs(Y[j][i]) < ymot and abs(V[j][i]) < capv):
            #Y[j],Z[j] = np.linspace(0,0,Nsteps), np.linspace(0,0.0001,Nsteps)
            nnn = 2
            n_ = np.append(n_, nnn)
        else:
            nnn = 0
            n_ = np.append(n_, nnn)
print('done Loop 2')
N = np.reshape(n_, (Natoms, Nsteps))

#print(n_)
#rint(N)
print('begin loop 3')
N0 = 0
for j in range(Natoms):
    NN = False
    for i in range(Nsteps):
        if N[j][i] == 2:
            NN = True
    if NN == True:
        N0 += 1

print('Number Captured =', N0)
stop = timeit.default_timer() 
print('Run Time =',round(stop - start, 3),'sec')   

for i in range(Natoms):
    'A plot for each of the Natoms particles'
    th = 0.5
    col = (0.1, float(i/(Natoms+1)+0.0001), 1-float(i/(Natoms+5)+0.0001))
    ax1.plot(Z[i],Y[i],linewidth=th, color = col)
    ax2.plot(Z[i],V[i],linewidth=th, color = col)

ax1.axhspan(-ymot/2,ymot/2, alpha=0.05, color='green')
ax1.axvspan(z0-0.01,z0+0.01, alpha=0.05, color='purple')
ax1.axvline(x = z0 - aa, color = 'k', linestyle='dotted')
ax1.axvline(x = z0, color = 'k', linestyle='dashed')
ax1.axvline(x = z0-0.01, color = 'r',linewidth=3)
ax1.axvline(x = z0+0.01, color = 'r',linewidth=3)
ax1.axhline(y = ymot/2, color = 'r',linewidth=3)
ax1.axhline(y = -ymot/2, color = 'r',linewidth=3)
ax1.set_ylim(top = 2*ymot, bottom = -2*ymot)

# ax1.axvline(x = wab-Da, color = 'r')
# ax1.axvline(x = wab-Db, color = 'r')
# ax1.axvline(x = wab-Dc, color = 'r')
# ax1.axvline(x = wab-Dd, color = 'r')
# ax1.axvline(x = wab-De, color = 'r')
# ax1.axvline(x = wab-Df, color = 'r')
# ax1.axvline(x = wab-Dg, color = 'r')
# ax1.axvline(x = wab, color = 'b')


ax2.axvspan(z0-0.01,z0+0.01, alpha=0.05, color='purple')
ax2.axhspan(-capV,capV, alpha=0.05, color='b')
ax2.axhspan(-capv,capv, alpha=0.05, color='red')
ax2.axvline(x = z0 - aa, color = 'k', linestyle='dotted')
ax2.axvline(x = z0, color = 'k', linestyle='dashed')
ax2.axvline(x = z0-0.01, color = 'r',linewidth=3)
ax2.axvline(x = z0+0.01, color = 'r',linewidth=3)
ax2.axhline(y = capv, color = 'r',linewidth=3)
ax2.axhline(y = -capv, color = 'r',linewidth=3)
ax2.set_xlim(left=0, right=z0+1*aa)
ax2.set_ylim(top=500, bottom=-20)



fig.subplots_adjust(hspace=0) # Makes the plots that share the 
#                             # same x axis on top of each other
ax1.set_ylabel("$y$ $coordinate$ / $m$", size = 19)
ax2.set_ylabel("$z$ $coordinate$ / $ms$`'", size = 19)
#ax3.set_title('Multi-Frequency Slowing Simulation: $\it{7}$ $\it{Frequencies}$, $\it{MOT}$ $\it{Magnetic}$ $\it{Field}$', size=17)
#ax1.set_title('Multi-Frequency: 7 Frequencies & Capture Number ', size=17)
ax2.set_xlabel('$z$ $coordinate$ / $m$ (0 = atom source)', size = 19)
#ax1.set_yticks(np.arange(-0.002, 0.002, step=0.0005))

q_pc = q*N0/Natoms
#Q_pc = q_pc/(1-np.cos(np.arctan(vy/v_)))
print('Total Flux % =', q_pc/Q * 100)
  
from datetime import date
today = date.today()
d4 = today.strftime("%d-%b-%Y")
#print("d4 =", d4)
#ax3.legend(title='     {}nIntensity = {}W/cm2nDetuning = {} w/ Increment {}MHznE0 = {}    no. atoms = {} nLength of Tube = {}cmnMag Field Gradient = {}G/cm'.format(d4,round(IrE, 3),Da/1000000,d/1000000,E0,nj, round((z0-aa)*100,3),round(Grad*1000000,2), loc=2, prop={'size': 18}))
textstr = '     {}\nIntensity = {}mW/cm2\nDetuning = {}\n w/ Increment {}MHz\nE0 = {}    no. atoms = {} \nLength of Tube = {}cm\nMag Field Gradient = {}G/cm\nLoading Rate* = {}'.format(d4,round(IrE*1000, 3),Da/1000000,d/1000000,E0,Natoms, round((z0-aa)*100,3),'Na', Q*N0/Natoms)#round(Grad*1000000,2))
bigstring = 'Hole Diameter = {} mm\nS/Atoms Captured = {} %\n Total Flux = {}%%%'
ax2.text(z0+0.05, 100, textstr, fontsize=14)
ax2.text(z0-0.18, 150,bigstring.format(Col_Gap*1000,round((N0/Natoms)*100,3), round(q_pc/Q*1000000,3)),fontweight='bold',fontsize=18)

stop = timeit.default_timer()
print('Velocity Range = [{},{}]'.format(round(min(vran),1),round(max(vran),1)))
print('# Particles = {}'.format(Natoms))
print('Beam Intensity = {}W/cm^2'.format(round(IrE, 3)))
print('Run Time =',round(stop - start, 3),'sec')
print('vy = {}m/s'.format(vy))
print('Flux ')

plt.show()
