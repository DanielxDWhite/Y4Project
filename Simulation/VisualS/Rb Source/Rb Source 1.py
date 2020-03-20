print('Rubidium Source THE 1ND')
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sts
import scipy.constants as sc
import timeit
start = timeit.default_timer()
T = 20
r_max = 0.04
tcc = 0.001

Col_zPos = 0.04
Col_Gap = 0.03
z_0 = 0.09
ymot = 0.008

Natoms = 8
Nsteps = 30
h = 0.05

def z_initial(Natoms):
    z0 = np.random.rand(Natoms)
    return Col_zPos*z0

def y_initial(Natoms): 
    r0 = np.random.rand(Natoms)-0.5
    return r_max*r0*0.22

kb = sc.Boltzmann
M = 87*sc.proton_mass
a=abs(np.sqrt((kb*(273+T))/(M)))
def z_vel(Natoms):
    x = np.random.rand(Natoms)
    #y = np.random.rand(Natoms)   
    #for i in range(Natoms):
    #    pm.append((-1)**round(y[i]))
    #a=abs(np.sqrt((kb*(273+T))/(M)))
    
    X = sts.maxwell.isf(x, scale=.09)
   # return np.multiply(pm,X)
    return X

a=abs(np.sqrt((kb*(273+T))/(M)))
print('scale, a', a)


def y_vel(Natoms, PM):
    x = np.random.rand(Natoms)
    #y = np.random.rand(Natoms)
    #for i in range(Natoms):
    #    pm_.append((-1)**round(y[i]))
    #a=abs(np.sqrt((kb*(T+273)))/(M))
    X = sts.maxwell.isf(x, scale=0.07)
    return np.multiply(PM,X)

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

def dy(t,z,v):
    return v

def dz(t,z,v):
    return v  

Y  = y_initial(Natoms)
pm = np.sign(Y)
PM=np.multiply(pm,-1)
Vy = y_vel(Natoms, PM)
Vz = z_vel(Natoms)
Z  = z_initial(Natoms)

"""    L O O P    """
#col=np.zeros(Natoms)
col = []
#for i in range(Natoms):
col.append('red')
th=[]#np.zeros(Natoms)
n=[]#zeros(Natoms)
m=[]#np.zeros(Natoms)

#print(RK4step(3,5,2,5,4,6))
zs,vs,ts=[],[],[]
ys,yvs=[],[]
xD = []
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
        
        z1,v1 = RK4step(ti,zi,vi,h,dv,dz)
        y1,yv1 = RK4step(ti,yi,yvi,h,dv,dz)

        yvi = yv1
        yi = y1
        zi = z1
        vi = v1

        if ( Col_zPos+tcc > zi > Col_zPos and abs(yi) < Col_Gap/2 ):
            #  Green
            xD_ = 1
            xD.append(xD_)
            
        elif ( abs(yvi) < 10 and Col_zPos+tcc > zi > Col_zPos and abs(yi) < Col_Gap/2 ):
            # Cyan
            xD_ = -1
            xD.append(xD_)
        else:
            #  Red
            xD_ = 0
            xD.append(xD_)
print(xD)


#                 NEED TO MAKE IT SO THAT IT DOESNT OVERWRITE THE COLOR
x_D = np.reshape(xD, (Natoms,Nsteps))
print(x_D )
Y_data = np.reshape(ys, (Natoms,Nsteps))
Vy_data = np.reshape(yvs, (Natoms,Nsteps))
Z_data = np.reshape(zs, (Natoms,Nsteps))

v_ = np.mean(z_vel(Natoms))
vy = ((ymot-Col_Gap/2)*v_)/(z_0**2+(ymot-Col_Gap/2)**2)**0.5
print("vy' is ",vy)
#print(Vy)
#y_ = 1
Nn = np.sum(n)
Mm = np.sum(m)
#print(th)
#print(col)
print(Natoms,Nn,Mm)
for j in range(Natoms):
    if x_D[j][Natoms-1] == 0 :
        col.append('red')
        th.append(0.3)
    if max(x_D[j]) > 1:
        col.append('green')
        th.append(2.0)
        n.append(1)
    else:
        col.append('cyan')
        th.append(4.0)
        m.append(1)
            
for i in range(Natoms):
    'A plot for each of the Natoms particles'
#    col = (0.1, float(i/(Natoms+1)+0.0001), 1-float(i/(Natoms+5)+0.0001))
    plt.plot(Z_data[i],Y_data[i],linewidth = th[i], color = col[i])
    

#def v_r_inital(Natoms):

#k = (np.random.rand, 0, np.random.rand)
'''Illustrative Borders PC '''
plt.bar(-tcc,         2*r_max+tcc,       width=tcc,      bottom= -r_max-tcc, color='k') #back
plt.bar(-tcc,         tcc,               width=Col_zPos+tcc, bottom= r_max,      color='k') #Top
plt.bar(-tcc,         tcc,               width=Col_zPos+tcc, bottom= -r_max-tcc, color='k') #Bot
plt.bar(Col_zPos, r_max-Col_Gap+tcc,     width=tcc,      bottom= Col_Gap,    color='k') #Top Col
plt.bar(Col_zPos, r_max-Col_Gap+tcc,     width=tcc,      bottom= -r_max-tcc, color='k') #Bot Col
plt.xlim(left=2*Col_zPos/3, right=Col_zPos+8*tcc)
plt.ylim(bottom=-r_max-0.5*r_max, top=r_max+0.5*r_max)
stop = timeit.default_timer()
leString = " nAtoms = {}\nnSteps/Size={}/{}\n Escaped = {} %\n <vy' = {}%\n Runtime = {}s".format(Natoms, Nsteps,h,round(Nn*100/Natoms,2), round(Mm*100/Natoms,3), round(stop - start, 3))
BIGString = " Pin Hole Diameter = {}mm \n Ratio of Simulated / Escaped = {}%".format(Col_Gap*1000, round(Mm*100/Nn),2)
plt.text(-Col_zPos/2+0.02,-r_max*2, leString, fontsize=19, bbox = dict(boxstyle='round', fc=(0.57,0.44,0.86), alpha=0.2))
plt.text(-Col_zPos/2+0.02,r_max*2,BIGString,fontweight='bold',fontsize=19,  bbox = dict(boxstyle='round', fc=(0.57,0.44,0.86), alpha=0.2))
plt.title('Ruubidium Vapour - Trajectory & Collimation', fontsize=30)
plt.xlabel('z -----> ', fontsize=25)
plt.ylabel('<----- r ----->', fontsize=25)
#plt.scatter(z_initial(Natoms),y_initial(Natoms), s=20, marker=">", color='c')
plt.show()

