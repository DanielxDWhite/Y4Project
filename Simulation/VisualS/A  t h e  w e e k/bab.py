import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


Raw = pd.read_excel('allDansdata.xlsx')
F,N,C = [],[],[]

for index,row in Raw[[44]].iterrows():
    N.append(row[0])
for index,row in Raw[[43]].iterrows():
    F.append(row[0])
for index,row in Raw[[45]].iterrows():
    C.append(row[0])
N = np.asarray(N)
F = np.asarray(F)
C = np.asarray(C)
N = N[~ np.isnan(N)]
F = F[~ np.isnan(F)]
C = C[~ np.isnan(C)]
print(len(N),len(F),C)
for i in range(1,len(C)):
    plt.axhline(C[i],c='grey', linewidth=7, alpha=0.3)
plt.axhline(C[0],c='grey', linewidth=7, alpha=0.3, label='Beamless Controls')
plt.axvline(115, linewidth=10, c='pink',alpha=0.7, label='Transition Resonance')
plt.scatter(F,N,c='b',s=300,label='3Vp-p  5MHz  220mW')
plt.legend(loc = 'upper left')
plt.title('The Best Atom Number Curve',size=19)
plt.xlabel('( Frequency + 384228 ) / MHz',size=17)
plt.ylabel('Number of Trapped Atoms',size=17)
plt.show()