import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# Natoms = 2000
Beta = [0,3,6,9,12,15,18,21]#26,30,35,40]
#)path = r"C:\\Users\\Daniel White\\Beta_data500000.csv"  # 5 0's

L = len(Beta)

def Data(a):
    '''[0] is data, [1] is Beta value'''
    df = pd.read_excel('EOM_Beta_data{}.xlsx'.format(a))
    return df,df.columns[1]

HallD,BallD = [],[]
for i in range(L):
    B = []
    H = []
    for index,row in Data(Beta[i])[0].iterrows():
        B.append(row[0])
        H.append(row[1])
    BallD.append(B)
    HallD.append(H)

plt.title('Density Distribution of Varying Beta\nw/ 1000 Atoms, 335MHz AOM, RF 5MHz ',size=19)
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

for i in range(0,L):
    col = ( float((i/(L+1)+0.0001)**1.5), float((i/(L+1)+0.0001)**0.5), 1-0.2*float((i/(L+5)+0.0001)**0.7))
    plt.plot(BallD[i], HallD[i],c=col,linewidth=5,label='Beta = '+str(Beta[i]/10)+', Vp-p = '+str(Beta[i]/20)+'V  {}cm'.format(round(FWHMs[i]*100,4)))
    plt.fill_between(BallD[i], HallD[i],color=col,alpha=0.4)
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
#ax.plot_wireframe(BallD,HallD,Beta)
