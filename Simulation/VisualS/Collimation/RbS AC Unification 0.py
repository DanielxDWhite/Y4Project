'''
       % captured as fn of Col_Gap
'''
import numpy as np 
import scipy.stats as sts
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.special as scp
import timeit
start = timeit.default_timer()

ymot = 0.008      # Diameter of Col Atoms

G = 38.11e6     # See ln 96

Xi_D = 6.6
Xi_d = 0.15
E0 = 170

aa = 0.15       # Coil Radius
s = 0.11        # Coil Separation
Curr = 0.6805   # Current (16G/cm @ 0.6805)
z0 = 0.24        # Position of MOT centre
z_0 = z0-aa

Satoms = 15000
Ssteps = 140
Sh = 0.000001
h = 0.000012     #Step Size
Natoms = 500
Nsteps = 500

T = 40
r_max = 0.04
Col_zPos = 0
Col_Gap = 0.001
tcc = 0.003
'''If ^ > ymot then simulation is unphysical'''

#J = 10e20
#Q = J*(Col_Gap/2)**2

print('E0, Da, d is {}, {}, {}'.format(E0,Xi_D,Xi_d))
Da = Xi_D*G    # Original detuning
d = Xi_d*G    # Increment of Detuning
Db,Dc,Dd,De,Df,Dg = Da+d, Da+2*d, Da+3*d, Da+4*d, Da+5*d, Da+6*d
xx = np.linspace(0.000001, 1, 100)

def Sz_initial(Satoms):
    '''Initial conditions for z position'''
    z0 = np.random.rand(Satoms)
    return Col_zPos*z0

def Sy_initial(Satoms): 
    '''Initial conditions for y position'''
    r0 = np.random.rand(Satoms)-0.5
    return r_max*r0*2

pm=[]
kb = sc.Boltzmann
MMM= 87*sc.proton_mass

a=abs(np.sqrt((kb*(273+T))/(MMM)))
print('scale, a', a)
#pm_=-1*np.sign(y_initial(Satoms))


def Szy_vel(Satoms,PM):
    z0 = np.random.rand(Satoms)
    y0 = 1-z0                                           # vy's share of the random number
    Z = sts.maxwell.isf(z0, scale=a)
    Y = sts.maxwell.isf(y0, scale=a*np.log(np.pi*abs(y0-0.5)))   # *1 no 2 from the e^(sqrt(av**2))
    return Z,np.multiply(PM,Y)


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

def Sdv(t,z,v):
    return 0

def Sdz(t,z,v):
    return v
########################### 
def zgen(n):
    lin = np.linspace(0,0,n)
    return lin

def yrand(n):
    ran = np.random.random(n)
    return (ran-0.5)*Col_Gap

def MBrand(n):
    x = np.random.rand(n)   
    X = sts.maxwell.isf(x, scale=a)
    return X
pm = np.random.rand(Satoms).round(0)
PM = np.power(np.linspace(-1,-1,Satoms),pm)  # This gives us a 1d array of +-1 randomly
#v_ = 0.5*(max(MBrand(Natoms))+min(MBrand(Natoms))) #Mean Velocity
sv_ = np.mean(Szy_vel(Satoms, PM)[0])
svy = ((ymot-Col_Gap/2)*sv_)/(z_0**2+(ymot-Col_Gap/2)**2)**0.5
v_ = np.mean(MBrand(Natoms))
vy = ((ymot-Col_Gap/2)*v_)/(z0**2+(ymot-Col_Gap/2)**2)**0.5
def Vyrand(n):
    ran = np.random.random(n)
    return (ran-0.5)*vy*2

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
''' Variable Dance '''
Rabi = dip*E0/hbar               # Rabi Frequency 
IoIs = 2*Rabi**2/G**2            # Intensity / Saturation Intensity
IrE = c*8.85e-12/2*E0**2/10000   # Intensity (This /10000 makes it W/cm^2)
w = wab - Dd                     # Average Freq of colliding photon
Lambda = 2*pi*c/w                  # Avg Wavelength
k = 2*pi/Lambda                  # Average wavenumber of a momentum transfering photon
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
"""    S O U R C E    L O O P    """
sY  = Sy_initial(Satoms)
sZ  = Sz_initial(Satoms)

#print(PM)
sVz, sVy = Szy_vel(Satoms,PM)
zs,vs,ts=[],[],[]
ys,yvs=[],[]
for j in range(Satoms):
    vi = sVz[j] 
    zi = sZ[j]
    yvi= sVy[j]
    yi = sY[j]
    for i in range(Ssteps):
        ti=Sh*i 
        zs.append(zi)
        vs.append(vi)
        ts.append(ti)
        ys.append(yi)
        yvs.append(yvi)
        z1,v1=RK4step(ti,zi,vi,Sh,Sdv,Sdz)
        y1,yv1=RK4step(ti,yi,yvi,Sh,Sdv,Sdz)
        yvi,yi,zi,vi = yv1,y1,z1,v1
        
Y_data = np.reshape(ys, (Satoms,Ssteps))
Vy_data = np.reshape(yvs, (Satoms,Ssteps))
Z_data = np.reshape(zs, (Satoms,Ssteps))

col = []
th = []
n = []
m = []

for j in range(Satoms):
    col.append('red')
    th.append(0.1)
    n.append(0)
    m.append(0)
    for i in range(Ssteps):
        nnn=False
        if ( Col_zPos+tcc > Z_data[j][i] > Col_zPos and abs(Y_data[j][i]) < Col_Gap/2 ):
            col[j] = 'green'
            th[j] = 2.0
            n[j] = 1
        if ( abs(Vy_data[j][i]) < svy and Col_zPos+tcc > Z_data[j][i] > Col_zPos and abs(Y_data[j][i]) < Col_Gap/2 ):
            col[j] = 'cyan'
            th[j] = 4.0
            m[j] = 1
        else:
            pass
Nn = np.sum(n)
Mm = np.sum(m)
#leString = " nAtoms = {}\nnSteps/Size={}/{}\n Escaped = {} %\n <vy' = {}%\n Runtime = {}s".format(Satoms, Ssteps,h,round(Nn*100/Satoms,2), round(Mm*100/Satoms,3), round(stop - start, 3))
#BIGString = " Pin Hole Diameter = {}mm \n Ratio of Simulated / Escaped = {}%".format(Col_Gap*1000, round(Mm*100/Nn),3)
print('sRed={}, sGreen={}, sCyan = {}'.format(Satoms-Nn-Mm,Nn,Mm))

#                 A T O M    C O U N T E R          i know what I need to do here # # # # # # # # # # # # # # # # # # # # # # # ## # # # 
J = 10e20
Q = J*(Col_Gap/2)**2
q = Q*(Mm/Satoms)
print('q , Q =',q,Q)
    
plt.close('all')
fig = plt.figure()
ax1 = plt.subplot2grid((2,1), (0,0))#, rowspan=2)
ax2 = plt.subplot2grid((2,1), (1,0), sharex=ax1)

"""   P L O T    &    C A P T U R E     L O O P  """

zlin=zgen(Natoms)
yran = yrand(Natoms)
Vyran = Vyrand(Natoms)
vran = MBrand(Natoms)

zsC,vsC,tsC=[],[],[]
ysC,yvsC=[],[]
"""this loop goes through all the atoms we've got and applies the force dv to them for a number of steps, Nsteps"""
for j in range(Natoms):
    viC = vran[j] 
    ziC = zlin[j]
    yviC= Vyran[j]
    yiC = yran[j]
    for i in range(Nsteps):
        tiC=h*i 
        zsC.append(ziC)
        vsC.append(viC)   
        tsC.append(tiC)     
        
        ysC.append(yiC)
        yvsC.append(yviC)
        
        z1C=RK4step(tiC,ziC,viC,h,dv,dz)[0]
        v1C=RK4step(tiC,ziC,viC,h,dv,dz)[1]
        y1C=RK4step(tiC,yiC,yviC,h,Sdv,Sdz)[0]
        yv1C=RK4step(tiC,yiC,yviC,h,Sdv,Sdz)[1]
        
        yviC = yv1C
        yiC = y1C
        ziC = z1C
        viC = v1C
        
Y = np.reshape(ysC, (Natoms,Nsteps))
V = np.reshape(vsC, (Natoms,Nsteps))
Z = np.reshape(zsC, (Natoms,Nsteps))
tt = np.array(tsC)
thet = np.split(tt, Natoms)[1]

Top, Thicc = 0.002, 0.003
ax1.bar(Col_zPos, Top, Thicc, bottom= Col_Gap/2, color='k')
ax1.bar(Col_zPos,-Top, Thicc, bottom=-Col_Gap/2, color='k')
#print(Y, Y.shape)
    
'''    Collimation Collision Detection     '''
'''nn=0
for j in range(Natoms):
    for i in range(Nsteps):
        if (Z[j][i] < Col_zPos and abs(Y[j][i]) > Col_Gap/2):
            Y[j],Z[j] = np.linspace(0,0,Nsteps), np.linspace(0,-0.01,Nsteps)
            nn += 1
'''    
'''    Capture Detection     '''

z_ , z__ = z0 - ymot, z0 + ymot
y_ = 0.01
capV, capv = 50,15
n_ = []

for j in range(Natoms):
    for i in range(Nsteps):
        if (z_ < Z[j][i] < z__ and abs(Y[j][i]) < ymot and abs(V[j][i]) < capV):
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

ax1.axhspan(-ymot/2,ymot/2, alpha=0.05, color='green')
ax1.axvspan(z0-0.01,z0+0.01, alpha=0.05, color='purple')
ax1.axvline(x = z0 - aa, color = 'k', linestyle='dotted')
ax1.axvline(x = z0, color = 'k', linestyle='dashed')
ax1.axvline(x = z0-0.01, color = 'r')
ax1.axvline(x = z0+0.01, color = 'r')
ax1.axhline(y = ymot/2, color = 'r')
ax1.axhline(y = -ymot/2, color = 'r')
ax1.set_ylim(top = 2*ymot, bottom = -2*ymot)

ax2.axvspan(z0-0.01,z0+0.01, alpha=0.05, color='purple')
ax2.axhspan(-capV,capV, alpha=0.05, color='b')
ax2.axhspan(-capv,capv, alpha=0.05, color='red')
ax2.axvline(x = z0 - aa, color = 'k', linestyle='dotted')
ax2.axvline(x = z0, color = 'k', linestyle='dashed')
ax2.axvline(x = z0-0.01, color = 'r')
ax2.axvline(x = z0+0.01, color = 'r')
ax2.axhline(y = capv, color = 'r')
ax2.axhline(y = -capv, color = 'r')
ax2.set_xlim(left=-0.009, right=z0+1*aa)
ax2.set_ylim(top=350, bottom=-20)



fig.subplots_adjust(hspace=0) # Makes the plots that share the 
#                             # same x axis on top of each other
ax1.set_ylabel("y coordinate / m", size = 17)
ax2.set_ylabel("Velocity / ms`'", size = 17)
#ax3.set_title('Multi-Frequency Slowing Simulation: $\it{7}$ $\it{Frequencies}$, $\it{MOT}$ $\it{Magnetic}$ $\it{Field}$', size=17)
ax1.set_title('7-Frequency: Total Loading Fraction', size=20)
ax2.set_xlabel('Distance / m', size = 18)
#ax1.set_yticks(np.arange(-0.002, 0.002, step=0.0005))

q_pc = q*N0/Natoms
#Q_pc = q_pc/(1-np.cos(np.arctan(vy/v_)))
print('Total Flux % =', q_pc/Q * 100)
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
stop = timeit.default_timer()
#print("d4 =", d4)
#ax3.legend(title='     {}\nIntensity = {}W/cm2\nDetuning = {} w/ Increment {}MHz\nE0 = {}    no. atoms = {} \nLength of Tube = {}cm\nMag Field Gradient = {}G/cm'.format(d4,round(IrE, 3),Da/1000000,d/1000000,E0,nj, round((z0-aa)*100,3),round(Grad*1000000,2), loc=2, prop={'size': 18}))

stop = timeit.default_timer()
Ustr = 'Pinhole Diameter = {}mm\n Captured / Released = {}%\n Loading Rate = N/A'.format(Col_Gap*1000,round(N0*Mm/(Ssteps*Nn)*100,4))
ustr = 'E0 = {}  D = {}  d = {}\n\n#Rb = {}\nEscaped = {}\nReasonable = {}\n\n#Simmed = {}\nCaptured = {}\n\nTube = {}cm\n Mag Field Grad = N/A G/m\n\nRuntime = {}\n  {}'.format(E0,Xi_D,Xi_d, Satoms,Nn,Mm, Natoms,N0, z_0*100,     round(stop - start, 3),d4)
ax2.text(z0+0.04,20,ustr, fontsize=19, bbox = dict(boxstyle='round', fc=(0.79,0.98,0.6), alpha=0.4))
ax1.text(z0+0.03,-1.6*ymot,Ustr,fontweight='bold',fontsize=23,  bbox = dict(boxstyle='round', fc=(0.99,0.87,0.12), alpha=0.5))

print('Velocity Range = [{},{}]'.format(round(min(vran),1),round(max(vran),1)))
print('# Particles = {}'.format(Natoms))
print('Beam Intensity = {}W/cm^2'.format(round(IrE, 3)))
print('Run Time =',round(stop - start, 3),'sec')
print('vy = {}m/s'.format(vy))
#print('Flux ')

plt.show()
