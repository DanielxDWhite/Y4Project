import numpy as np 
import matplotlib.pyplot as plt
import scipy.constants as sc

""" Variables """
E0 = 120
d = 1.5e7
B_0=0.16

def vgen(n):
    '''Initial Velocities'''
    lin = np.linspace(-5,5,n)
    #lin = 3*np.ones((20,1))
    return lin

def zgen(n):
    '''Inital Positions'''
    lin = np.linspace(0.5,-0.5,n)
    #lin = 0*np.ones((20,1))
    return lin

def RK4step(ti,zi,vi,h,dv,dz):             
    k11=dz(ti,zi,vi)
    k21=dv(ti,zi,vi,B_0)    
    k12=dz(ti+h/2,zi +(h/2)*k11,vi +(h/2)*k21)
    k22=dv(ti+h/2,zi +(h/2)*k11,vi +(h/2)*k21,B_0)    
    k13=dz(ti+h/2,zi +(h/2)*k12,vi +(h/2)*k22)
    k23=dv(ti+h/2,zi +(h/2)*k12,vi +(h/2)*k22,B_0)    
    k14=dz(ti+h,zi +(h)*k13,vi +(h)*k23)
    k24=dv(ti+h,zi +(h)*k13,vi +(h)*k23,B_0)    
    z1=zi+(h/6.0)*(k11+2.0*k12+2.0*k13+k14)  
    v1=vi+(h/6.0)*(k21+2.0*k22+2.0*k23+k24)               
    zi = z1
    vi = v1    
    return zi,vi

""" Physical & Atomic Constants """
kb=sc.Boltzmann
mu0 = sc.mu_0
muB = 9.2740099*10**-24
u=sc.proton_mass
hbar=sc.hbar
c=sc.c
pi=np.pi 
e=sc.e
M=86.9*u
wab=2*pi*384.23e12
G=38.11e6
Z =337
dip= 3.485e-29

''' Variable Orgy '''
Rabi = dip*E0/hbar
IoIs = 2*Rabi**2/G**2
IrE = c*8.85e-12/2*E0**2
w = wab - d
Lambda=2*pi*c/w
k = 2*pi/Lambda

i=0
zs=[]
vs=[]
ts=[]
a=0.000000000000001

print(mu0)
#print()
#Ir=power/(a**2*pi)

def dv(t,z,v,B_0):
    'Incremental Acceleration'
    fz = abs(z)
    O = w/(2*pi*c)-muB*B_0*fz/hbar
    c1 = IoIs**2+4*d**2/G**2
    c2 = O*d*8/G**2
    c3 = 4*O**2/G**2   
    rhoaa =  -IoIs**2/(c1+c2*v+c3*v**2) + IoIs**2/(c1-c2*v+c3*v**2)   
    return rhoaa*hbar*k*G/M

def dz(t,z,v):
    return v  

plt.close('all')
fig = plt.figure()
ax1 = plt.subplot2grid((2,1), (0,0))
ax2 = plt.subplot2grid((2,1), (1,0),sharex=ax1)
fig.subplots_adjust(hspace=0)


"""step size"""		
h=0.002
"""number of atoms"""
nj=20
"""number of steps"""   
ni=900
"""creation of our array of velocities"""
#vis=Jgen(T,nj,M)
vlin=vgen(nj)
zlin=zgen(nj)
"""this loop goes through all the atoms we've got and
applies the force dv to them for a number of steps ni"""
for j in range(nj):
 #   vi=vis[j]
    vi = vlin[j]
    zi = zlin[j]
    for i in range(ni):
        ti=a+h*i 
        zs.append(zi)
        vs.append(vi)   
        ts.append(ti)     
        z1=RK4step(ti,zi,vi,h,dv,dz)[0]
        v1=RK4step(ti,zi,vi,h,dv,dz)[1]
        zi = z1
        vi = v1
        
V = np.reshape(vs, (nj,ni))
Z = np.reshape(zs, (nj,ni))
tt = np.array(ts)
thet = np.split(tt, nj)[1]

ax1.plot(thet,V[0])
ax1.plot(thet,V[1])
ax1.plot(thet,V[2])
ax1.plot(thet,V[3])
ax1.plot(thet,V[4])
ax1.plot(thet,V[5])
ax1.plot(thet,V[6])
ax1.plot(thet,V[7])
ax1.plot(thet,V[8])
ax1.plot(thet,V[9])
ax1.plot(thet,V[10])
ax1.plot(thet,V[11])
ax1.plot(thet,V[12])
ax1.plot(thet,V[13])
ax1.plot(thet,V[14])
ax1.plot(thet,V[15])
ax1.plot(thet,V[16])
ax1.plot(thet,V[17])
ax1.plot(thet,V[18])
ax1.plot(thet,V[19])

ax2.plot(thet,Z[0])
ax2.plot(thet,Z[1])
ax2.plot(thet,Z[2])
ax2.plot(thet,Z[3])
ax2.plot(thet,Z[4])
ax2.plot(thet,Z[5])
ax2.plot(thet,Z[6])
ax2.plot(thet,Z[7])
ax2.plot(thet,Z[8])
ax2.plot(thet,Z[9])
ax2.plot(thet,Z[10])
ax2.plot(thet,Z[11])
ax2.plot(thet,Z[12])
ax2.plot(thet,Z[13])
ax2.plot(thet,Z[14])
ax2.plot(thet,Z[15])
ax2.plot(thet,Z[16])
ax2.plot(thet,Z[17])
ax2.plot(thet,Z[18])
ax2.plot(thet,Z[19])

ax1.set_title('Optical Molasses Simulation w/ Real Parameters', size=18)
ax1.set_ylabel('Velocity', size = 16)
ax2.set_ylabel('Coordinate', size = 16)
ax2.set_xlabel('Time', size = 16)
plt.legend(title='B-field Grad = {}G/cm\nIntensiy = {}W/cm2\nE0 = {}\nDetuning = {}'.format(B_0*10,IrE,d,E0), loc=8)
#ax1.set_xticks([])
plt.show()
print('Intensity = {}W/cm2'.format(IrE/1000))