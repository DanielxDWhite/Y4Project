print('H33 - Focused on a Spike  w/ Intercepts & Slope Analysis')
print('      -  -  -  -  -  -  -  -  -  -  -  -  - ')
import matplotlib.pyplot as plt
import numpy as np
import math
import timeit
import scipy.constants as sc
start = timeit.default_timer()

'''Archetecture'''
n = 4000                                           # No. Steps (Even please)
l_of_t = 200                                        # Time Range (Even please)
t = np.linspace(-l_of_t/2,l_of_t/2,n)                # Time / x axis

'''Atom/Light Paramters'''
d = 10**7                                            # Detuning
Xi = 0.0001         #***                            # Ξ -> Rabi Freqency (Metcalf uses 1,3,5 for Edbif Criss-Cross)
PhaR, PhaB = 0,3                                   # Initial Phase Differences
OrF, ObF = 384.2*10**12-d, 384.2*10**12+d         # Freq of Red & Blue
p = 3.584*10**-29                                # Electric Dipole = <r>e

h = sc.hbar                                    # Reduced Planck's Constant
O_0 = d * Xi                                  # d-less Rabi Freq
O_ = d*(Xi**2+1)**0.5                        # Omega'
E0 = h*O_0/p                                # Elecric Field

'''Intercepts'''
Int_Thresh = 0.0000000001

if PhaR == 0:
   PhaseR = 0
else: 
   PhaseR = np.pi/PhaR

if PhaB == 0:
   PhaseB = 0
else: 
   PhaseB = np.pi/PhaB

def Or(O_, OrF, PhaseR, i):
    '''Rabi 'Influence' of Red Electric Field. By defalt Red has the phase difference'''
    o_r = O_*np.sin(i*OrF*np.pi*(n/10000000000)/(n)+PhaseR)
    return float(o_r)
    
def Ob(O_, ObF, PhaseB, i):
    '''Rabi 'Influence' of Blue Electric Field'''
    ob = O_*np.sin((i-0.000001)*ObF*np.pi*(n/10000000000)/(n)+PhaseB)   
    return float(ob)

E11,E12,E13,E14,E15,E16,E17,E18,E19,E20,E21,E22,E23 = np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0)), np.zeros((1,0))
Er, Eb = np.zeros((1,0)),np.zeros((1,0))

def H33(Or, Ob, d, t):# 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3
    '''33x33 Hamiltonian taken from Source "EbifExpansion"'''
    Haux =np.array([[16*d,Ob/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],#1
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
    return h*np.linalg.eigvals(Haux)
    
for i in range(int(-n/2),int(n/2)):   
    i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,i32,i33 = np.split(H33(Or(O_, OrF, PhaseR, i), Ob(O_, ObF, PhaseB, i), d, t), 33)
  
    E16 = np.append(E16, i16)
    E17 = np.append(E17, i17)
    E18 = np.append(E18, i18)
    
    Er = np.append(Er, Or(O_, OrF, PhaseR, i))
    Eb = np.append(Eb, Ob(O_, ObF, PhaseB, i))

plt.close('all')
 
fig = plt.figure()
ax1 = plt.subplot2grid((4,1), (0,0), rowspan=3)
ax2 = plt.subplot2grid((4,1), (3,0), sharex=ax1)
fig.subplots_adjust(hspace=0)
th = 1.6

ax1.set_title('Energies of States in 33x33 Bichromatic Dressed Hamiltonian')   # TOP PlOT
ax1.plot(t, E16/(h*d), color='yellowgreen',linewidth=th)
ax1.plot(t, E17/(h*d), color='seagreen',linewidth=th)
#ax1.plot(t, E18, color='turquoise',linewidth=th)

EpR = np.square(Er)                                                             # Intercepts
EpB = np.square(Eb)        
xWb = np.where(EpR < Int_Thresh)
xWr = np.where(EpB < Int_Thresh)
NoIntR = len(xWr[0])
NoIntB = len(xWb[0])
print('No. intercepts = {}          when looking for E^2 < {} with the biggest = {} (Blue)'.format(NoIntB, Int_Thresh,  max(EpB[xWb[0]])))  
for i in range(0, NoIntR):
    ax1.axvline(x = ((float(xWr[0][i])-n/2)*l_of_t/n), linestyle=':',color='red')
    ax2.axvline(x = ((float(xWr[0][i])-n/2)*l_of_t/n), linestyle=':',color='red')
for i in range(0, NoIntB):
    ax1.axvline(x = ((float(xWb[0][i])-n/2)*l_of_t/n), linestyle=':',color='b', label='Node of Red Field')
    ax2.axvline(x = ((float(xWb[0][i])-n/2)*l_of_t/n), linestyle=':',color='b')
    

print(' -     -      -      -      -      -      -')
ax2.plot(t, Er, color='r')                                                     # BOTTOM PlOT
ax2.plot(t, Eb, color='b')
ax2.fill_between(t, 0, Er, color='red', alpha=0.2)
ax2.fill_between(t, 0, Eb, color='blue', alpha=0.2)

plt.show
ax2.set_xlabel('Distance / a.u.', size = 14)
ax1.set_ylabel('Energy by 'r'$\hbar\delta$', size = 14)
ax2.set_ylabel('Electric Field of Standing Wave', size = 14)

#ax2.set_yticks(minor=True)
plt.show

plt.close

'''Slope Finding'''          # Firstlly, need to find correct Int__Tresh (line 24) so that len xWb is 1                        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

fur = 40                     # index left of grey region / nodal point
NodeInx = int(xWb[0])

'''                         #ENERGY LEVEL SLOPE
s1,s2,s3,s4,s5,s6,s7,s8,s9,s10 = E16[NodeInx-fur],E16[NodeInx-fur+4],E16[NodeInx-fur+8],E16[NodeInx-fur+12],E16[NodeInx-fur+16],E16[NodeInx-fur+20],E16[NodeInx-fur+24],E16[NodeInx-fur+28],E16[NodeInx-fur+32],E16[NodeInx-fur+36]
dE = np.sum(np.diff([s1,s2,s3,s4,s5,s6,s7,s8,s9,s10]))
dt = 36*l_of_t/n
dEdt = np.abs(dE/dt)
                             ELECTRIC FIELD SLOPE

s1,s2,s3,s4,s5,s6,s7,s8,s9,s10 = Eb[NodeInx-fur],Eb[NodeInx-fur+4],Eb[NodeInx-fur+8],Eb[NodeInx-fur+12],Eb[NodeInx-fur+16],Eb[NodeInx-fur+20],Eb[NodeInx-fur+24],Eb[NodeInx-fur+28],Eb[NodeInx-fur+32],Eb[NodeInx-fur+36]
dE = np.sum(np.diff([s1,s2,s3,s4,s5,s6,s7,s8,s9,s10]))
dt = 36*l_of_t/n

dEdt = np.abs(dE/dt)
'''

dOdt = O_*2*np.pi*ObF*np.sin(PhaseB)
Pna = np.exp((-np.pi*O_**2)/(2*dOdt))

#print("From {} index's away, the slope for 10 index's is dE / dt  =  {}".format(fur,dOdt))

ax1.axvspan((NodeInx-fur-n/2)*l_of_t/n, (NodeInx-fur-n/2+36)*l_of_t/n, alpha=0.2, color='gray', label='Region of Slope Calulation')
ax2.axvspan((NodeInx-fur-n/2)*l_of_t/n, (NodeInx-fur-n/2+36)*l_of_t/n, alpha=0.2, color='gray')
ax1.legend(title_fontsize=15,
           title='Ξ = {}\nno. steps = {}\nInit. Phase of:\n  Blue=π/{} \n  Red=π/{}\nΔFreq = {}\n Prob(Non-Adiabatic Transition) = {}'.format(Xi,n,PhaB,PhaR,np.abs(OrF-ObF),Pna),
           shadow=True,
           fancybox=True,
           loc=3)
def OoMag(x):
    xOoM= int(math.log10(x))
    return xOoM
print(' _       _       _      _      _')
print('H33 Close Up Plot:')
print('                   Time Steps = {}'.format(n))
#print('                     detuning = {}'.format(d))
stop = timeit.default_timer()
print('                     Run Time =',round(stop - start, 3),'sec')
#print('No. intercepts = {} when looking for E^2 < {} with the biggest = {} (Blue)'.format(NoIntB, Int_Thresh,  max(EpB[xWb[0]])))  
#print('              Or = {}'.format(o_r))
#print('              Ob = {}'.format(ob))
print('')
rR = O_**2/dOdt
print("O'^2/dOdt = {} >> 1 if in Adiabatic Limit".format(rR))
print('Probability of a Non-Adiabatic Transition = {}'.format(Pna))
print('')
print('Intensity = {}  [{}]'.format(E0**2,OoMag(E0**2)))
print("Omega'^2 = {}  [{}]".format(O_**2,OoMag(O_**2)))
print("dO/dt = {}  [{}]".format(dOdt, OoMag(dOdt)))