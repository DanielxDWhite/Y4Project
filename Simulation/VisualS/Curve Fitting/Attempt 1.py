import numpy as np 
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.special as scp
'''
Bessel for EOM
'''
RF = 1
WHAT = 1
Beta = WHAT * RF
x = np.linspace(0,5,100)


j0,j1,j2 = scp.jv(0,Beta),scp.jv(1,Beta),scp.jv(2,Beta)
J0,J1,J2 = scp.jv(0,x), scp.jv(1,x), scp.jv(2,x)
plt.plot(x,abs(J0))#, label=j0)
plt.plot(x,abs(J1))#, label=j1)
plt.plot(x,abs(J2))#, label=j2)

EfromY = 1
E0,E1,E2 = EfromY*j0,EfromY*j1,EfromY*j2

plt.axvline(x = RF, color = 'orange')
plt.axhline(y = j0,color= 'blue')
plt.axhline(y = j1,color= 'green')
plt.axhline(y = j2,color= 'red')


STR = 'Beta = {}\n J0 = {}\n J1 = {}\n J2 = {}'.format(Beta,round(j0,2),round(j1,2),round(j2,2))
plt.text(3,0.8,STR)
strI = '\n I0 = {}\n I1 = {}\n I2 = {}'.format(round(E0,2)**2,round(E1,2)**2,round(E2,2)**2)
plt.text(4,0.8,strI)
plt.show()