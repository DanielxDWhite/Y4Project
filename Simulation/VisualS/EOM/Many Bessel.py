# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 22:24:38 2020

@author: white
"""

import numpy as np 
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.special as scp
'''
Bessel for EOM
'''
RF = np.linspace(0,10,100)#0.2
WHAT = 1
Beta = WHAT * RF
x = np.linspace(0,5,100)


j0,j1,j2,j3,j4,j5,j6,j7,j0a,j1a,j2a,j3a,j4a,j5a,j6a,j7a,j0b,j1b,j2b = scp.jv(0,Beta),scp.jv(1,Beta),scp.jv(2,Beta),scp.jv(3,Beta),scp.jv(4,Beta),scp.jv(5,Beta),scp.jv(6,Beta),scp.jv(7,Beta),scp.jv(8,Beta),scp.jv(9,Beta),scp.jv(10,Beta),scp.jv(11,Beta),scp.jv(12,Beta),scp.jv(13,Beta),scp.jv(14,Beta),scp.jv(15,Beta),scp.jv(16,Beta),scp.jv(17,Beta),scp.jv(18,Beta)
J0,J1,J2 = scp.jv(0,x), scp.jv(1,x), scp.jv(2,x)

#plt.plot(x,abs(J0))#, label=j0)
#plt.plot(x,abs(J1))#, label=j1)
#plt.plot(x,abs(J2))#, label=j2)

plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2)+2*(j6**2)+2*(j7**2)+2*(j0a**2)+2*(j1a**2)+2*(j2a**2)+2*(j3a**2)+2*(j4a**2)+2*(j5a**2)+2*(j6a**2)+2*(j7a**2)+2*(j0b**2)+2*(j1b**2)+2*(j2b**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2)+2*(j6**2)+2*(j7**2)+2*(j0a**2)+2*(j1a**2)+2*(j2a**2)+2*(j3a**2)+2*(j4a**2)+2*(j5a**2)+2*(j6a**2)+2*(j7a**2)+2*(j0b**2)+2*(j1b**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2)+2*(j6**2)+2*(j7**2)+2*(j0a**2)+2*(j1a**2)+2*(j2a**2)+2*(j3a**2)+2*(j4a**2)+2*(j5a**2)+2*(j6a**2)+2*(j7a**2)+2*(j0b**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2)+2*(j6**2)+2*(j7**2)+2*(j0a**2)+2*(j1a**2)+2*(j2a**2)+2*(j3a**2)+2*(j4a**2)+2*(j5a**2)+2*(j6a**2)+2*(j7a**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2)+2*(j6**2)+2*(j7**2)+2*(j0a**2)+2*(j1a**2)+2*(j2a**2)+2*(j3a**2)+2*(j4a**2)+2*(j5a**2)+2*(j6a**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2)+2*(j6**2)+2*(j7**2)+2*(j0a**2)+2*(j1a**2)+2*(j2a**2)+2*(j3a**2)+2*(j4a**2)+2*(j5a**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2)+2*(j6**2)+2*(j7**2)+2*(j0a**2)+2*(j1a**2)+2*(j2a**2)+2*(j3a**2)+2*(j4a**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2)+2*(j6**2)+2*(j7**2)+2*(j0a**2)+2*(j1a**2)+2*(j2a**2)+2*(j3a**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2)+2*(j6**2)+2*(j7**2)+2*(j0a**2)+2*(j1a**2)+2*(j2a**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2)+2*(j6**2)+2*(j7**2)+2*(j0a**2)+2*(j1a**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2)+2*(j6**2)+2*(j7**2)+2*(j0a**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2)+2*(j6**2)+2*(j7**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2)+2*(j6**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2)+2*(j5**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2)+2*(j4**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2)+2*(j3**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2)+2*(j2**2),linewidth=2 )
plt.plot(RF, j0**2+2*(j1**2),linewidth=2 )
plt.plot(RF, j0**2,linewidth=2 )

plt.xlabel(r'$ \beta $',size=26)
plt.ylabel('$Total$ $Output$ $Intensity$',size=23)

EfromY = 1
E0,E1,E2 = EfromY*j0,EfromY*j1,EfromY*j2
plt.show()
'''
plt.axvline(x = RF, color = 'orange')
plt.axhline(y = j0,color= 'blue')
plt.axhline(y = j1,color= 'green')
plt.axhline(y = j2,color= 'red')
'''
#print(Beta, j0**2+j1**2+j2**2)
'''
STR = 'Beta = {}\n J0 = {}\n J1 = {}\n J2 = {}'.format(Beta,round(j0,2),round(j1,2),round(j2,2))
plt.text(3,0.8,STR)
strI = '\n I0 = {}\n I1 = {}\n I2 = {}\n SUM = {}'.format(round(E0,2)**2,round(E1,2)**2,round(E2,2)**2, j0**2+j1**2+j2**2)
plt.text(4,0.8,strI)
plt.show()
'''