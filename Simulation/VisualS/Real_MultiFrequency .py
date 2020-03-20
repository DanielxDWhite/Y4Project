import numpy as np 
import matplotlib.pyplot as plt
import scipy.constants as sc


""" Variables """
E0 = 90000
d = 2e7
B_0=0.000000001

def vgen(n):
    '''Initial Velocities'''
    lin = np.linspace(-0.015,0.017,n)
    #lin = 0*np.ones((20,1))
    return lin

def zgen(n):
    '''Inital Positions'''
    lin = np.linspace(0.00005,-0.00005,n)
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
IrR = 2*Rabi**2/G**2
IrE = c*8.85e-12/2*E0**2
w = wab - d
Lambda=2*pi*c/w
k = 2*pi/Lambda

i=0
zs=[]
vs=[]
ts=[]
a=1
b=5
#print()
#print()
#Ir=power/(a**2*pi)

def dv(t,z,v,B_0):
    fz = abs(z)
    O = w/(2*pi*c)-mu0*B_0*fz/hbar    # Approximate - needs g-factors & which F?
    c1 = 1+IrR+4*d**2/G**2
    c2 = O*d*8/G**2
    c3 = 4*O**2/G**2   
    rhoaa =  -IrR/(c1+c2*v+c3*v**2) + IrR/(c1-c2*v+c3*v**2)   
    return rhoaa*hbar*k*G/M

def dz(t,z,v):
    return v  

plt.close('all')
fig = plt.figure()
ax1 = plt.subplot2grid((2,1), (0,0))
ax2 = plt.subplot2grid((2,1), (1,0),sharex=ax1)
fig.subplots_adjust(hspace=0)


"""step size"""		
h=0.00005
"""number of atoms"""
nj=20
"""number of steps"""   
ni=300
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
        
V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20 = np.split(np.array(vs), nj)
Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,Z10,Z11,Z12,Z13,Z14,Z15,Z16,Z17,Z18,Z19,Z20 = np.split(np.array(zs), nj)
tt = np.array(ts)
thet = np.split(tt, nj)[1]
#print(thet)
ax1.plot(thet,V1)
ax1.plot(thet,V2)
ax1.plot(thet,V3)
ax1.plot(thet,V4)
ax1.plot(thet,V5)
ax1.plot(thet,V6)
ax1.plot(thet,V7)
ax1.plot(thet,V8)
ax1.plot(thet,V9)
ax1.plot(thet,V10)
ax1.plot(thet,V11)
ax1.plot(thet,V12)
ax1.plot(thet,V13)
ax1.plot(thet,V14)
ax1.plot(thet,V15)
ax1.plot(thet,V16)
ax1.plot(thet,V17)
ax1.plot(thet,V18)
ax1.plot(thet,V19)
ax1.plot(thet,V20)

ax2.plot(thet,Z1)
ax2.plot(thet,Z2)
ax2.plot(thet,Z3)
ax2.plot(thet,Z4)
ax2.plot(thet,Z5)
ax2.plot(thet,Z6)
ax2.plot(thet,Z7)
ax2.plot(thet,Z8)
ax2.plot(thet,Z9)
ax2.plot(thet,Z10)
ax2.plot(thet,Z11)
ax2.plot(thet,Z12)
ax2.plot(thet,Z13)
ax2.plot(thet,Z14)
ax2.plot(thet,Z15)
ax2.plot(thet,Z16)
ax2.plot(thet,Z17)
ax2.plot(thet,Z18)
ax2.plot(thet,Z19)

ax2.plot(thet,Z20)

ax1.set_title('Optical Molasses Simulation w/ Real Parameters', size=18)
ax1.set_ylabel('Velocity', size = 16)
ax2.set_ylabel('Coordinate', size = 16)
ax2.set_xlabel('Time', size = 16)
ax1.legend(title='B-field Grad = {}\nE0 = {}'.format(B_0,E0), loc=8)
#ax1.set_xticks([])
plt.show()
print('Intensity = {}W/cm2'.format(IrE/1000))