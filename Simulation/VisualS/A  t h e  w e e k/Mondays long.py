import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
fig = plt.figure()
ax1 = plt.subplot2grid((2,1), (0,0), colspan=1)
ax2 = plt.subplot2grid((2,1), (1,0),sharex=ax1 )
ax1.grid()
ax2.grid()
fig.subplots_adjust(hspace=0)

Raw = pd.read_excel('allDansdata1.xlsx') 


F1,N1,C1 = [],[],[]
for index,row in Raw[[48]].iterrows():
    N1.append(row[0])
for index,row in Raw[[47]].iterrows():
    F1.append(row[0])
for index,row in Raw[[49]].iterrows():
    C1.append(row[0])
N1 = np.asarray(N1)
F1 = np.asarray(F1)
C1 = np.asarray(C1)
N1 = N1[~ np.isnan(N1)]
F1 = F1[~ np.isnan(F1)]
C1 = C1[~ np.isnan(C1)]
print('1',len(N1),len(F1),len(C1))
for i in range(1,len(C1)):
    ax1.axhline(C1[i],c='grey', linewidth=7, alpha=0.3)
ax1.axhline(C1[0],c='grey', linewidth=7, alpha=0.3)
ax1.axvline(1115, linewidth=10, c='pink',alpha=0.7)
ax1.scatter(F1,N1,c='b',s=300,label='4.9Vp-p  5MHz  7.1mW')
ax1.legend(loc = 'lower left',prop={'size': 14})
ax1.set_title('The Best Atom Number Curve',size=19)

F2,N2,C2 = [],[],[]
for index,row in Raw[[52]].iterrows():
    N2.append(row[0])
for index,row in Raw[[51]].iterrows():
    F2.append(row[0])
for index,row in Raw[[53]].iterrows():
    C2.append(row[0])
N2 = np.asarray(N2)
F2 = np.asarray(F2)
C2 = np.asarray(C2)
N2 = N2[~ np.isnan(N2)]
F2 = F2[~ np.isnan(F2)]
C2 = C2[~ np.isnan(C2)]
print('2',len(N2),len(F2),len(C2))
for i in range(1,len(C2)):
    ax2.axhline(C2[i],c='grey', linewidth=7, alpha=0.3)
ax2.axhline(C2[0],c='grey', linewidth=7, alpha=0.3, label='Beamless Controls')
ax2.axvline(1115, linewidth=10, c='pink',alpha=0.7, label='Transition Resonance')
ax2.scatter(F2,N2,c='purple',s=300,label='2Vp-p  5MHz  30mW')

ax2.legend(loc = 'upper left',prop={'size': 14})

ax2.set_xlabel('( Frequency + 384227 ) / MHz',size=17)
plt.ylabel('Number of Trapped Atoms',size=17)





plt.show()