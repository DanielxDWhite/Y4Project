
print('Velocity Distribution')
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

T = 300

u = sc.proton_mass         # Proton Mass
m = 87*u                   # Mass of 87Rb

vy = 500
vz = 2000

def MaxBoltz_3D(vx,vy,vz, t):
    N = (m/(2*np.pi*sc.k*t))**(3/2)
    exp = -m*(vx**2+vy**2+vz**2)/(2*sc.k*T)
    return N*np.exp(exp)
'''
def Trap_3D(MaxBoltz_3D):
    np.trapz(np.trapz(np.trapz(np.trapz(MaxBoltz_3D,)),vy),vx)
    trap0 = MaxBoltz_3d(...)
    for i in range(3):
        trap1 = np.trapz(trap0)
        trap0 = trap1
    
    return trap1
    '''
'''
for i in range(len(vINC)):
    plt.plot(np.trapz(MaxBoltz_3D(np.linspace(1,vINC[i],100),vy,vz)) , np.linspace(1,vINC[i],len(vINC)))
'''

Data = []
BData = []

TT=30
t = np.linspace(3,300, TT)

vx = 5000
Nsteps = 300
vINC = np.linspace(0,vx,Nsteps)
plt.close('all')
fig = plt.figure()
ax1 = plt.subplot2grid((1,2), (0,1))
ax2 = plt.subplot2grid((1,2), (0,0), sharey=ax1)
n=-1
for i in t:
    n += 1
    n=float(n)
    for v in vINC:
        a = np.trapz(MaxBoltz_3D(np.split(vINC,TT)[n][v],vy,vz,i), vx)
        Data = np.append(Data, a)
    BData = np.append(BData, Data)
    
print(BData)

D = np.reshape(Data, (TT, Nsteps))
print(D)
plt.close()
#print(vINC)
for i in range(len(t)):
    plt.plot(vINC,D[i])
    

    
ax2.set_xlabel('Velocity m/s')
ax2.set_ylabel('Probability')
plt.show()