import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Natoms = 1500
RF = [5,10,25,50,100,150,200,250,300,400,500,600]
#path = r"C:\\Users\\Daniel White\\RF_data500000.csv"  # 5 0's

def Data(a):
    '''[0] is data, [1] is RF value'''
    df = pd.read_excel('RF_data{}00000.xlsx'.format(a))
    return df,df.columns[1]


HallD,BallD = [],[]
for i in range(len(RF)):
    B = []
    H = []
    for index,row in Data(RF[i])[0].iterrows():
        B.append(row[0])
        H.append(row[1])
    BallD.append(B)
    HallD.append(H)

plt.title('Density Distribution of Increasing RF\nw/ 150 Atoms, 150MHz AOM, Beta 1, Intens = 170 (4mW/cm2)',size=19)
plt.xlabel('Distance from source / m',size=17)
plt.ylabel('Particle Density',size=17)

# # # # # # # # # # ##
'''  F W H M  ''' # ##
FWHMs = []
for i in range(len(RF)):
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

col_ = ( float(0.3), float((0/(len(RF)+1)+0.0001)**2), 1-float((0/(len(RF)+5)+0.0001)**0.7))
plt.fill_between(BallD[0], HallD[0],color=col_,alpha=0.2)
plt.plot(BallD[0], HallD[0],c=col_,linewidth=3,label='0.5MHz {}cm = FWHM'.format(round(FWHMs[0]*100,4)))
for i in range(1,len(RF)):
    col = ( float(0.3), float((i/(len(RF)+1)+0.0001)**2), 1-float((i/(len(RF)+5)+0.0001)**0.7))
    plt.plot(BallD[i], HallD[i],c=col,linewidth=5,label=str(RF[i]/10)+'MHz {}cm'.format(round(FWHMs[i]*100,4)))
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
#ax.plot_wireframe(BallD,HallD,RF)





