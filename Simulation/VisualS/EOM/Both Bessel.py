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


plt.plot(RF, abs(j0)+2*abs(j4),linewidth=2 )
plt.plot(RF, abs(j0)+2*abs(j2),linewidth=2 )
plt.plot(RF, abs(j0),linewidth=2 )

# plt.plot(RF, abs(j2),linewidth=2 )
# plt.plot(RF, abs(j1),linewidth=2 )
# plt.plot(RF, abs(j0),linewidth=2 )

plt.xlabel('$Beta$',size=19)
plt.ylabel('$INCOMING$ $LaTeX$ x 2',size=19)

EfromY = 1
E0,E1,E2 = EfromY*j0,EfromY*j1,EfromY*j2
plt.show()