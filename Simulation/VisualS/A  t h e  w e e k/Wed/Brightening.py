'''  G O O D   D A T A '''
import numpy as np 
import matplotlib.pyplot as plt
fig = plt.figure()
ax1 = plt.subplot2grid((2,1), (0,0), colspan=1)
ax2 = plt.subplot2grid((2,1), (1,0),sharex=ax1 )
# a

Beamless = [2.8,2.7,2.66,2.69,2.74,2.79,2.8,2.77]
for i in range(1,len(Beamless)):
    ax1.axhline(Beamless[i],c = 'grey', linewidth=5, alpha=0.3)
ax1.axhline(Beamless[0],c = 'grey',label = 'Contol (No Slowing Beam)', linewidth=4, alpha=0.5)


ax1.axvline(115, c='pink',linewidth=14,label='Transition Resonance')
ax2.axvline(115, c='pink',linewidth=14,label='Transition Resonance')

a = [2.85,2.59,2.63,2.71,2.7,2.58,2.49, 3.43,3.36,3.32,3.52,3.28,3.5]
af= [112,112,112,110,109,110,115, 91,91,92,90,90,90]

b = [3.62,3.47,3.55,3.69,3.56,3.58,3.48,3.47,3.59,3.55,3.79,3.79,3.34,3.23,3.3, 2.85,3.1,3.24,3.27]
bf= [85,85,85,82,80,81,77,77,77,73,72,73,56,52,53, 40,43,42,40]

c = [2.73,2.76,2.83,2.87,2.96, 2.71,2.81,2.87,2.79,2.79,2.83,2.77,2.7,2.87]
cf= [23,23,22,10,8, -4,-3,-1,-24,-30,-22,-34,-42,-41]

d = [2.84,2.88,2.74,3.3,3.24,3.2,2.76,3.28,2.99,2.78,3.46,3.38,3.32,3.41,3.27,3.54,3,3.05,2.8,2.83,2.83,2.81,2.57,3.09]
df= [107,106,107,94,90,93,107,94,95,105,98,85, 86,86,86,87,100,100,100,103,102,103,112,96]


A = [2.71,2.74,2.74,2.84,2.87,2.76,2.84,2.79,2.82]
B = [2.96,2.91,2.77,2.91,2.94,2.95,2.78,2.92,2.94]
C = [2.94,3.21,3.17,3.14,3.14,2.92,2.94]
D = [3.54,3.42,3.41,3.18,3.48,3.6,3.06,3.26]
Af= [113,115,116,113,105,101,97,93,91]
Bf= [90,90,91,87,86,81,81,78,77]
Cf= [72,58,57,65,66,73,72]
Df= [42,20,20,2,42,37,-20,26]

#print(len(a),len(af),len(b),len(bf),len(c),len(cf),len(d),len(df))

ax1.scatter(af,a,c='r',s=240)
ax1.scatter(bf,b,c='r',s=240)
ax1.scatter(cf,c,c='r',s=240)
ax1.scatter(df,d,c='r',s=240)

ax1.scatter(Af,A,c='b',s=200)
ax1.scatter(Bf,B,c='b',s=200)
ax1.scatter(Cf,C,c='b',s=200)
ax1.scatter(Df,D,c='b',s=200)

ax1.set_title('The first sign of MOT brightness increasing as a result of our slower',size=21)
ax1.set_xlabel('Frequency / MHz + const.',size=20)
ax1.set_ylabel('# Trapped Atoms x10^6',size=20)


ax1.grid()
ax2.grid()
fig.subplots_adjust(hspace=0)
# # # # # # # # # # #  S I M U L A T I O N  # # # # # # # # # # 
s1_= [-40,-20,-10,0,5,10,15,20,25,28,30,33,35,38,40,43,45,50,60,70,80,100]
S1 = [5,7,11,18,22,27,35,44,58,71,81,89,93,95,99,88,71,35,15,10,7,5]
s2_= [-10,0,10,20,30,40,50,60,70,80,90,100, 35,45,48,55,57,65]
S2 =[21,27,33,40,46,58,72,50,40,16,9,7, 51,64,70,71,64,52]

s1 = np.negative(s1_)+115
ax2.scatter(s1,S1,c='r',s=240,label='3Vp-p  15MHz  53.3mW')
s2 = np.negative(s2_)+115
ax2.scatter(s2,S2,c='b',s=200,label='4.2Vp-p  15MHz  53.3mW')
ax2.axhline(0,c = 'grey',label = 'Contol (No Slowing Beam)', linewidth=4, alpha=0.5)
ax2.legend(loc='upper left',fontsize=15)#,title='3Vp-p  15MHz  53.3mW')
plt.show()