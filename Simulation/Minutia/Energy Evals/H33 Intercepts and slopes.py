print('H33 -  x Intercepts & Slopes')

import matplotlib.pyplot as plt
import numpy as np
import timeit
start = timeit.default_timer()

n = 3998                 # No. Steps (Even please)
t = np.linspace(-100,100,n)  #time / x axis

d = 300                   #Detuning

Xi = 1.816742               # Ξ -> Rabi Freqency (Metcalf uses 1,3,5 for Edbif Criss-Cross)

PhaR = 0
PhaB = 2.5                  #Initial Phase Difference

OrF = 1.3                 #Freq of Red
ObF = 1.3               #Freq of Blue

O_0 = d * Xi            # Raw Rabi Freqency

if PhaR == 0:
   PhaseR = 0
else: 
   PhaseR = np.pi/PhaR

if PhaB == 0:
   PhaseB = 0
else: 
   PhaseB = np.pi/PhaB

def Or(O_0, OrF, PhaseR, i):
    '''Rabi 'Influence' of Red Electric Field. By defalt Red has the phase difference'''
    o_r = O_0*np.sin(i*OrF*np.pi/(n)+PhaseR)
    return float(o_r)
    
def Ob(O_0, ObF, PhaseB, i):
    '''Rabi 'Influence' of Blue Electric Field'''
    ob = O_0*np.sin((i-0.0001)*ObF*np.pi/(n)+PhaseB)   
    return float(ob)

E11,E12,E13,E14,E15,E16,E17,E18,E19,E20,E21,E22,E23 = np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0))
Er, Eb = np.zeros((1,0)),np.zeros((1,0))

def H33(Or, Ob, d, t):# 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3
    '''33x33 Hamiltonian taken from Source "EbifExpansion"'''
    Haux = np.array([[16*d,Ob/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#1
                     [Ob/2,15*d,Or/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#2-1
                     [0,Or/2,14*d,Ob/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#3
                     [0,0,Ob/2,13*d,Or/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#4-2
                     [0,0,0,Or/2,12*d,Ob/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#5
                     [0,0,0,0,Ob/2,11*d,Or/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#6-3
                     [0,0,0,0,0,Or/2,10*d,Ob/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#7
                     [0,0,0,0,0,0,Ob/2,9*d,Or/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#8-4
                     [0,0,0,0,0,0,0,Or/2,8*d,Ob/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#9
                     [0,0,0,0,0,0,0,0,Ob/2,7*d,Or/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#10-5
                     [0,0,0,0,0,0,0,0,0,Or/2,6*d,Ob/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#11
                     [0,0,0,0,0,0,0,0,0,0,Ob/2,5*d,Or/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#12-6
                     [0,0,0,0,0,0,0,0,0,0,0,Or/2,4*d,Ob/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#13
                     [0,0,0,0,0,0,0,0,0,0,0,0,Ob/2,3*d,Or/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#14-7
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,Or/2,2*d,Ob/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#15
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,Ob/2,d,Or/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#16-8
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Or/2,0,Ob/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#17
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Ob/2,-d,Or/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#18-8
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Or/2,-2*d,Ob/2,0,0,0,0,0,0,0,0,0,0,0,0,0],#19
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Ob/2,-3*d,Or/2,0,0,0,0,0,0,0,0,0,0,0,0],#20-7
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Or/2,-4*d,Ob/2,0,0,0,0,0,0,0,0,0,0,0],#21
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Ob/2,-5*d,Or/2,0,0,0,0,0,0,0,0,0,0],#22-6
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Or/2,-6*d,Ob/2,0,0,0,0,0,0,0,0,0],#23
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Ob/2,-7*d,Or/2,0,0,0,0,0,0,0,0],#24-5
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Or/2,-8*d,Ob/2,0,0,0,0,0,0,0],#25
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Ob/2,-9*d,Or/2,0,0,0,0,0,0],#26-4
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Or/2,-10*d,Ob/2,0,0,0,0,0],#27
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Ob/2,-11*d,Or/2,0,0,0,0],#28-3
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Or/2,-12*d,Ob/2,0,0,0],#29
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Ob/2,-13*d,Or/2,0,0],#30-2
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Or/2,-14*d,Ob/2,0],#31
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Ob/2,-15*d,Or/2],#32-1
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Or/2,-16*d]])#33           #Ox/2 is from EdbifExp
    return np.linalg.eigvals(Haux)
    
#print( H13(1,2,3,4) )       Or/2,*d,Ob/2,
                                                                                                                                                 

for i in range(int(-n/2),int(n/2)):   
    i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,i32,i33 = np.split(H33(Or(O_0, OrF, PhaseR, i), Ob(O_0, ObF, PhaseB, i), d, t), 33)
  
    E16 = np.append(E16, i16)
    E17 = np.append(E17, i17)
    E18 = np.append(E18, i18)
    
    Er = np.append(Er, Or(O_0, OrF, PhaseR, i))
    Eb = np.append(Eb, Ob(O_0, ObF, PhaseB, i))

plt.close('all')
 
fig = plt.figure()
ax1 = plt.subplot2grid((4,1), (0,0), rowspan=3)
ax2 = plt.subplot2grid((4,1), (3,0), sharex=ax1)
fig.subplots_adjust(hspace=0)
th = 1.6

ax1.set_title('Energies of States in 33x33 Bichromatic Dressed Hamiltonian')
ax1.plot(t, E16, color='yellowgreen',linewidth=th)
ax1.plot(t, E17, color='seagreen',linewidth=th)
#ax1.plot(t, E18, color='turquoise',linewidth=th)

ax1.legend(title_fontsize=15,
           title='Ξ = {}\nno. steps = {}\nInit. Phase of:\n  Blue=π/{} \n  Red=π/{}\nΔFreq = {}'.format(Xi,n,PhaB,PhaR,np.abs(OrF-ObF)),
           shadow=True,
           fancybox=True,
           loc=3)
#ax1.axvline(x=22.058956)

''''
x_int = np.zeros((1,0))
for i in range(1, n+1):
    if Or(O_0, i) == 0 or Ob(O_0, i) == 0:    
        x_int = np.append(x_int, i)
        aa = i
        print(aa)                      FAILURE TO GET X INTERCEPT
print(x_int)               
y0 = np.zeros((1,0))       
y0= np.append(y0,np.split(x_int, len(x_int))) 
print(y0)
#for j in range(0,len(x_int)):
'''

ax2.plot(t, Er, color='r')
ax2.plot(t, Eb, color='b')
ax2.fill_between(t, 0, Er, color='red', alpha=0.2)
ax2.fill_between(t, 0, Eb, color='blue', alpha=0.2)

plt.show
ax2.set_xlabel('Distance / a.u.', size = 14)
ax1.set_ylabel('Energy / a.u.', size = 14)
ax2.set_ylabel('Electric Field of Standing Wave', size = 14)

#ax2.set_yticks(minor=True)
plt.show

plt.close

print('H33 Close Up Plot:')
print('                   Time Steps = {}'.format(n))
#print('                     detuning = {}'.format(d))
stop = timeit.default_timer()
print('                     Run Time =',round(stop - start, 3),'sec')
#print('              Or = {}'.format(o_r))
#print('              Ob = {}'.format(ob))