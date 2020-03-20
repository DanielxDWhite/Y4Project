import numpy as np 
import scipy.stats as sts
import matplotlib.pyplot as plt
import scipy.constants as sc

def Jgen(T,n,M):
    
    MB=sts.maxwell
    kb=1.38e-23
    a=abs(np.sqrt((kb*T)/(M)))
    vt=MB.rvs(loc=0, scale=a, size=n, random_state=None)    
    return vt
    '''
    for i in range(19):
        I1,I2,I3,I4,I5,I6,I7,I8,I9,I0,I11,I12,I13,I14,I15,I16,I17,I18,I19,I20
    '''
def vgen(n):
    lin = np.linspace(-1000,2000,n)
    return lin
    
def zgen(n):
    lin = np.linspace(100,-100,n)
    #lin = -10*np.ones((20,1))
    return lin
    
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
    zi =z1
    vi = v1    
    return zi,vi
    
"""Ordinary physical constants. unchanged"""
kb=1.38e-23
u=1.66e-27
h=6.626e-34
hbar=h/(2.0*np.pi)
c=3.00e8
pi=np.pi
"""unique constants of experiment. to be varied"""
n=20
T=200
Ir = 100
G=38.11e6
M=86.9*u
atoms = []
atoms1=[]
w=1.0e14
d = 9e8
wab= w - d
Lambda=1.0
i=0
"""non-BS constants"""
zs=[]
vs=[]
ts=[]

a=1
b=5
zi=1000 # This doesn't matter

#Ir=power/(a**2*pi)
def dv(t,z,v):        
    O = w/c
    c1 = 1+Ir+4*d**2/G**2
    c2 = O*2*d*4/G**2
    c3 = 4*O**2/G**2
    return -Ir/(c1+c2*v+c3*v**2)         
def dz(t,z,v):
    return v  
    
plt.close('all')
fig = plt.figure()
ax1 = plt.subplot2grid((2,1), (0,0))
ax2 = plt.subplot2grid((2,1), (1,0), sharex=ax1)

"""step size"""		
h=100    
"""number of atoms"""
nj=20
"""number of steps"""   
ni=500    
"""creation of our array of velocities"""
#vis=vgen(T,nj,M)
vlin=vgen(nj)
zlin=zgen(nj)
"""this loop goes through all the atoms we've got and applies the force dv to them for a number of steps ni"""
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
thet = np.split(tt, nj)[0]



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

#ax2 = plt.subplot2grid((2,1), (1,0), sharex=ax1)
 
fig.subplots_adjust(hspace=0)
ax1.set_ylabel('Velocity', size = 14)
ax2.set_ylabel('Coordinate', size = 14)
ax2.set_xlabel('Time', size = 14)
'''
def Force(v,T):
    Deff=(abs(w-wab))-(v*w/c)
    F=0.5*hbar*(w/c)*Lambda*Ir/(1.0+Ir+4*Deff**2/(Lambda**2))
    return F
'''
plt.show()
    
   