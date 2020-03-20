print('Simulation of Optical Bloch Equations')
import numpy as np
import matplotlib.pyplot as plt

def drho11(t,r12,r21,r22,gam,O_0,d):
    rh11 = 1j*0.5*O_0*np.conj(np.exp(-1j*d*t))*r12-1j*0.5*O_0*np.exp(-1j*d*t)*r21+2*gam*r22
    return rh11
    
def drho21(t,r11,r12,r22,gam,O_0,d):
    rh21 = -1j*0.5*O_0*np.conj(np.exp(-1j*d*t))*(r11-r22)-gam*r12
    return rh21
    
def drho12(t,r11,r12,r22,gam,O_0,d):
    rh21 = 1j*0.5*O_0*np.conj(np.exp(1j*d*t))*(r11-r22)-gam*r12
    return rh21    

def drho22(t,r12,r21,r22,gam,O_0,d):
    rh11 = -1j*0.5*O_0*np.conj(np.exp(-1j*d*t))*r12+1j*0.5*O_0*np.exp(1j*d*t)*r21-2*gam*r22
    return rh11
    
#------    A R C E T E C T U R E 

N=1000
t=np.linspace(0,0.25e-6,N)
dt=t[2]-t[1]
     
print(t[10])


print(dt)
#------    I N I T I A L  C O N D I T I O N S
rho11=[]
rho22=[]
rho12=[]
rho21=[]
rho11=np.append(rho11,float(1))
rho22=np.append(rho22,float(0))
rho12=np.append(rho12,float(0))
rho21=np.append(rho21,float(0))
O_0=2*np.pi*10*1e6 
d=float(0)
gam=0.1*O_0 # linewidth

print(rho11)



def RK4_OBE(t,rho11,rho12,rho21,rho22):            
    '''Runge-Kutta 4 for the interconnected density matrix elements'''

    k1 = drho22(t[i-1],rho12[i-1],rho21[i-1],rho22[i-1],gam,O_0,d)*dt
    l1 = drho11(t[i-1],rho12[i-1],rho21[i-1],rho22[i-1],gam,O_0,d)*dt 
    m1 = drho12(t[i-1],rho11[i-1],rho12[i-1],rho22[i-1],gam,O_0,d)*dt
    n1 = drho21(t[i-1],rho11[i-1],rho12[i-1],rho22[i-1],gam,O_0,d)*dt

    k2 = drho22(t[i-1]+0.5*dt,rho12[i-1]+0.5*m1,rho21[i-1]+0.5*n1,rho22[i-1]+0.5*k1,gam,O_0,d)*dt 
    l2 = drho11(t[i-1]+0.5*dt,rho12[i-1]+0.5*m1,rho21[i-1]+0.5*n1,rho22[i-1]+0.5*k1,gam,O_0,d)*dt 
    m2 = drho12(t[i-1]+0.5*dt,rho11[i-1]+0.5*l1,rho12[i-1]+0.5*m1,rho22[i-1]+0.5*k1,gam,O_0,d)*dt 
    n2 = drho21(t[i-1]+0.5*dt,rho11[i-1]+0.5*l1,rho12[i-1]+0.5*m1,rho22[i-1]+0.5*k1,gam,O_0,d)*dt 
   
    k3 = drho22(t[i-1]+0.5*dt,rho12[i-1]+0.5*m2,rho21[i-1]+0.5*n2,rho22[i-1]+0.5*k2,gam,O_0,d)*dt 
    l3 = drho11(t[i-1]+0.5*dt,rho12[i-1]+0.5*m2,rho21[i-1]+0.5*n2,rho22[i-1]+0.5*k2,gam,O_0,d)*dt 
    m3 = drho12(t[i-1]+0.5*dt,rho11[i-1]+0.5*l2,rho12[i-1]+0.5*m2,rho22[i-1]+0.5*k2,gam,O_0,d)*dt 
    n3 = drho21(t[i-1]+0.5*dt,rho11[i-1]+0.5*l2,rho12[i-1]+0.5*m2,rho22[i-1]+0.5*k2,gam,O_0,d)*dt 
   
    k4 = drho22(t[i-1]+dt,rho12[i-1]+m3,rho21[i-1]+n3,rho22[i-1]+k3,gam,O_0,d)*dt 
    l4 = drho11(t[i-1]+dt,rho12[i-1]+m3,rho21[i-1]+n3,rho22[i-1]+k3,gam,O_0,d)*dt 
    m4 = drho12(t[i-1]+dt,rho11[i-1]+l3,rho12[i-1]+m3,rho22[i-1]+k3,gam,O_0,d)*dt 
    n4 = drho21(t[i-1]+dt,rho11[i-1]+l3,rho12[i-1]+m3,rho22[i-1]+k3,gam,O_0,d)*dt 
         
    rho22[i] = rho22[i-1]+((k1+2*(k2+k3)+k4)/6)
    rho11[i] = rho11[i-1]+((l1+2*(l2+l3)+l4)/6)
    rho12[i] = rho12[i-1]+((m1+2*(m2+m3)+m4)/6)
    rho21[i] = rho21[i-1]+((n1+2*(n2+n3)+n4)/6)    
    
    return rho11[i], rho12[i], rho21[i], rho22[i]#,k1


#print(RK4_OBE(12,34,5,6,7,8)[4])

    


ts = []
a = 1
for i in range(2,N+2):
    a=1
    h=5
    ti=a+h*i 
    Rho11i = 1
    Rho12i = 1
    Rho21i = 1
    Rho22i = 1
    rho11 = np.append(rho11,Rho11i)
    rho12 = np.append(rho12,Rho12i)
    rho21 = np.append(rho21,Rho21i)
    rho22 = np.append(rho22,Rho22i)
#    vs.append(vi)   
    ts = np.append(ts,ti)     
    Rho11_1 = RK4_OBE(ti,rho11,rho12,rho21,rho22)[0]
    Rho12_1 = RK4_OBE(ti,rho11,rho12,rho21,rho22)[1]
    Rho21_1 = RK4_OBE(ti,rho11,rho12,rho21,rho22)[2]
    Rho22_1 = RK4_OBE(ti,rho11,rho12,rho21,rho22)[3]
#    z1=RK4step(ti,zi,vi,h,dv,dz)[0]
#    v1=RK4step(ti,zi,vi,h,dv,dz)[1]
    Rho11i = Rho11_1
    Rho12i = Rho12_1
    Rho21i = Rho21_1
    Rho22i = Rho22_1
#    zi =z1
#    vi = v1
    plt.plot(ts,Rho11_1)
    
    
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
    zi =z1
    vi = v1    
    return zi,vi'''