import numpy as np 
import scipy.stats as sts
import matplotlib.pyplot as plt

def vgen(T,n,M): 
    MB=sts.maxwell
    kb=1.38e-23
    a=abs(np.sqrt((kb*T)/(M)))
    vt=MB.rvs(loc=0, scale=a, size=n, random_state=None)    
    return vt

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
    
#print(RK4step(1,2,3,4,5,6)[0])
"""Ordinary physical constants. unchanged"""
kb=1.38e-23
u=1.66e-27
h=6.626e-34
hbar=h/(2.0*np.pi)
c=3.00e8
pi=np.pi
"""unique constants of experiment. to be varied"""
n=20
T=500
M=86.9*u
atoms = []
atoms1=[]
w=1.0e14
wab=1.1e14
Lambda=1.0
i=0
"""BS constants"""
zs=[]
vs=[]
ts=[]
c1=7
c2=.005
c3=.005
c4=.0005
a=1
b=5
zi=1

#Ir=power/(a**2*pi)
def dv(t,z,v):        
    return -c1/(c2+c3*v+c4*v**2)         
def dz(t,z,v):
    return v  


"""step size"""		
h=5    
"""number of atoms"""
nj=1 
"""number of steps"""   
ni=200    
"""creation of our array of velocities"""
vis=vgen(T,nj,M)
"""this loop goes through all the atoms we've got and applies the force dv to them for a number of steps ni"""

for j in range(nj):
    vi=vis[j]    
    for i in range(ni):
        ti=a+h*i 
        zs.append(zi)
        vs.append(vi)   
        ts.append(ti)     
        z1=RK4step(ti,zi,vi,h,dv,dz)[0]
        v1=RK4step(ti,zi,vi,h,dv,dz)[1]
        zi =z1
        vi = v1
    plt.plot(ts,vs) 



"""   
def Force(v,T):
    Deff=(abs(w-wab))-(v*w/c)
    F=0.5*hbar*(w/c)*Lambda*Ir/(1.0+Ir+4*Deff**2/(Lambda**2))
    return F
    
"""    