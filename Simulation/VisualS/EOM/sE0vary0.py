import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# Natoms = 2000
E0 = [15,20,25,30,35,40,45,50]
#)path = r"C:\\Users\\Daniel White\\E0_data500000.csv"  # 5 0's

L = len(E0)

def Data(a):
    '''[0] is data, [1] is E0 value'''
    df = pd.read_excel('sE0_data{}0.xlsx'.format(a))
    return df,df.columns[1]

HallD,BallD = [],[]
for i in range(L):
    B = []
    H = []
    for index,row in Data(E0[i])[0].iterrows():
        B.append(row[0])
        H.append(row[1])
    BallD.append(B)
    HallD.append(H)

plt.title('Density Distribution of Increasing E0\nw/ 2000 Atoms, 150MHz AOM, 15MHz E0, Beta 1.5',size=19)
plt.xlabel('Distance from source / m',size=17)
plt.ylabel('Particle Density',size=17)

# # # # # # # # # # ##
'''  F W H M  ''' # ##
FWHMs = []
for i in range(L):
    HallD_ = HallD[i]
    BallD_ = BallD[i]
    Hpeak = [max(HallD_), int(round( np.median(np.where(HallD_==max(HallD_))) ))]

    Lmin = max(np.where(HallD_[:Hpeak[1]] == min(HallD_[:Hpeak[1]] ) )[0])

    Lmax = max(np.where(HallD_[Hpeak[1]:] == min(HallD_[Hpeak[1]:] ) )[0])
    #print('Index Distance from center to edge R L= ', BallD[Hpeak[1]]-BallD[Lmin], BallD[Lmax]-BallD[Hpeak[1]])
    #vLmin,vLmax = BinFactor*  IS IT POSSIBLE TO CONVERT INTO ORGINAL VELOCITY = not needed right now
    FWi = Lmax-Lmin
    Bot = max(HallD_[Lmax],HallD_[Lmin])

    HM = Bot + (Hpeak[0]+Bot)/2

    lHM = np.abs(HallD_[:Hpeak[1]]-HM).argmin()
    rHM = np.abs(HallD_[Hpeak[1]:]-HM).argmin()+Hpeak[1]

    #print(lHM,rHM)
    Skew = -1*(BallD_[Hpeak[1]]-BallD_[lHM]-BallD_[rHM]+BallD_[Hpeak[1]])
    #print('Skew =',Skew,'  +=MB')
    FWHM = BallD_[rHM]-BallD_[lHM]

    FWHMs.append(FWHM)
#print(FWHMs)
c = 3e8
def IrE(b):
    xD = c*8.85e-12/2*b**2/10000
    return xD
for i in range(1,L):
    col = ( float((i/(L+1)+0.0001)), float((i/(L+1)+0.0001)**2), 0.9-0.5*float((i/(L+5)+0.0001)**0.7))
    plt.plot(BallD[i], HallD[i],c=col,linewidth=5,label='E0 = '+str(E0[i]*10)+', I = '+str(round(IrE(E0[i]*10),3)*1000)+'mW/cm2  {}cm'.format(round(FWHMs[i]*100,4)))
    plt.fill_between(BallD[i], HallD[i],color=col,alpha=0.2)
plt.legend(fontsize=25)
plt.show()

#            FWHM is a bit broken, need way more bins

''' 
 3 D   A t t e m p t 
 #  #  #  #  #  #  #
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_suBallD(1,1,1,projection='3d')
'''
#ax.plot_wireframe(BallD,HallD,E0)
