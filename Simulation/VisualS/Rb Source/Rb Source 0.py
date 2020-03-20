print('Rubidium Source')
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sts
import scipy.constants as sc
import timeit
import math
start = timeit.default_timer()

r_max = 0.02
tcc = 0.003

Col_zPos = 0.06
Col_Gap = 0.002
z_0 = 0.18
ymot = 0.004   # RADIUS

Natoms = 100000
Nsteps = 300
h = 0.0000005

# Realsitic Source Effusion
T = 100     #Temp in Celcius    
T_ = 273 + T
P = 133.322*(3e-7/(273+25))*T_
kb = sc.k
M = 87 * sc.proton_mass
A = np.pi*(Col_Gap/2)**2
Q = P*A/((2*np.pi*M*kb*T)**0.5)

def z_initial(Natoms):
    '''Initial conditions for z position'''
    z0 = np.random.rand(Natoms)
    return Col_zPos*z0

def y_initial(Natoms): 
    '''Initial conditions for y position'''
    r0 = np.random.rand(Natoms)-0.5
    return r_max*r0*2

pm=[]
kb = sc.Boltzmann
M = 87*sc.proton_mass

a=abs(np.sqrt((kb*(273+T))/(M)))
print('scale, a', a)
#pm_=-1*np.sign(y_initial(Natoms))


def zy_vel(Natoms,PM):
    z0 = np.random.rand(Natoms)
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

def dv(t,z,v):
    return 0

def dz(t,z,v):
    return v  

Y  = y_initial(Natoms)
Z  = z_initial(Natoms)
pm = np.random.rand(Natoms).round(0)
PM = np.power(np.linspace(-1,-1,Natoms),pm)  # This gives us a 1d array of +-1 randomly
#print(PM)
Vz, Vy = zy_vel(Natoms,PM)


"""    L O O P    """
zs,vs,ts=[],[],[]
ys,yvs=[],[]
for j in range(Natoms):
    vi = Vz[j] 
    zi = Z[j]
    yvi= Vy[j]
    yi = Y[j]
    for i in range(Nsteps):
        ti=h*i 
        zs.append(zi)
        vs.append(vi)
        ts.append(ti)
        ys.append(yi)
        yvs.append(yvi)
        z1,v1=RK4step(ti,zi,vi,h,dv,dz)
        y1,yv1=RK4step(ti,yi,yvi,h,dv,dz)
        yvi,yi,zi,vi = yv1,y1,z1,v1
        
Y_data = np.reshape(ys, (Natoms,Nsteps))
Vy_data = np.reshape(yvs, (Natoms,Nsteps))
Z_data = np.reshape(zs, (Natoms,Nsteps))

col = []
th = []
n = []
m = []

v_ = np.mean(zy_vel(Natoms, PM)[0])
vy = ((ymot-Col_Gap/2)*v_)/(z_0**2+(ymot-Col_Gap/2)**2)**0.5
print('vzmean, vymax, vyavg= {}  {}  {}'.format(v_,vy,np.mean(zy_vel(Natoms, PM)[1])))
#print(Vy)
y_ = 100
for j in range(Natoms):
    col.append('red')
    th.append(0.1)
    n.append(0)
    m.append(0)
    for i in range(Nsteps):
        nnn=False
        if ( Col_zPos+tcc > Z_data[j][i] > Col_zPos and abs(Y_data[j][i]) < Col_Gap/2 ):
            col[j] = 'green'
            th[j] = 2.0
            n[j] = 1
        if ( abs(Vy_data[j][i]) < vy and Col_zPos+tcc > Z_data[j][i] > Col_zPos and abs(Y_data[j][i]) < Col_Gap/2 ):
            col[j] = 'cyan'
            th[j] = 4.0
            m[j] = 1
        else:
            pass


Nn = np.sum(n)
Mm = np.sum(m)
#print(th)
#print(col)
print(Natoms,Nn,Mm)
stop = timeit.default_timer() 
print(round(stop - start, 3))

plt.close('all')
for i in range(Natoms):
    'A plot for each of the Natoms particles'
#    col = (0.1, float(i/(Natoms+1)+0.0001), 1-float(i/(Natoms+5)+0.0001))
    plt.plot(Z_data[i],Y_data[i],linewidth = th[i], color = col[i])
    
'''
def OoM(number):
    return math.floor(math.log(number, 10))
ToFx = Q*Mm/Nn
print('Total Reasonable Flux = ',ToFx, OoM(ToFx))
'''
#k = (np.random.rand, 0, np.random.rand)

plt.bar(-tcc,         2*r_max+tcc,       width=tcc,      bottom= -r_max-tcc, color='k') #back
plt.bar(-tcc,         tcc,               width=Col_zPos+tcc, bottom= r_max,      color='k') #Top
plt.bar(-tcc,         tcc,               width=Col_zPos+tcc, bottom= -r_max-tcc, color='k') #Bot
plt.bar(Col_zPos, r_max-Col_Gap+tcc,     width=tcc,      bottom= Col_Gap,    color='k') #Top Col
plt.bar(Col_zPos, r_max-Col_Gap+tcc,     width=tcc,      bottom= -r_max-tcc, color='k') #Bot Col
#plt.xlim(left=-2*tcc, right=Col_zPos+8*tcc)
#plt.ylim(bottom=-r_max-2*tcc, top=r_max+2*r_max)
stop = timeit.default_timer()
print(round(stop - start, 3))
leString = " nAtoms = {}\nnSteps/Size={}/{}\n Green = {}% [{}]\n Cyan = {}% [{}]\n Runtime = {}s".format(Natoms, Nsteps,h,round(Nn*100/Natoms,2),Nn, round(Mm*100/Natoms,3),Mm, round(stop - start, 3))
BIGString = " Pin Hole Diameter = {}mm \n Ratio of Simulated / Escaped = {}% ".format(Col_Gap*1000, round(Mm*100/Nn,3))
plt.text(-Col_zPos/2+0.02,-r_max*2.3, leString, fontsize=19, bbox = dict(boxstyle='round', fc=(0.99,0.87,0.12), alpha=0.3))
plt.text(-Col_zPos/2+0.0,r_max*1.6,BIGString,fontweight='bold',fontsize=19,  bbox = dict(boxstyle='round', fc=(0.79,0.98,0.6), alpha=0.2))

#plt.title('Rubidium Vapour - Trajectory & Collimation', fontsize=30)
plt.xlabel('$z$ -----> ', fontsize=30)
plt.ylabel('<----- $r$ ----->', fontsize=30)

#plt.scatter(z_initial(Natoms),y_initial(Natoms), s=20, marker=">", color='c')
plt.show()
