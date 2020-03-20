import matplotlib.pyplot as plt
import numpy as np

Pa_1 = [3.1,3.2,0.0]#
Pa_3 = [3.8,5.5,6.8]
Pa_5 = [5.5,4.2,5.5]
Pb_1 = [4.2,1.6,4.8]#
Pb_3 = [3.5,5.0,4.8]
Pb_5 = [6.4,7.8,6.3]
Pc_1 = [2.0,4.25,2.0]#
Pc_3 = [1.8,3.1,3.3]
Pc_5 = [4.8,7.3,4.4]

P15 = [3.5,6.77,4.563]
z = [0.01,0.03,0.05]

z1 = [0.01,0.01,0.01]
z3 = [0.03,0.03,0.03]
z5 = [0.05,0.05,0.05]

plt.scatter(z1,Pa_1,s=60)
plt.scatter(z1,Pb_1,s=100,c='c')
plt.scatter(z1,Pc_1,s=180,c='y')


plt.scatter(z3,Pa_3,s=60)
plt.scatter(z3,Pb_3,s=100,c='c')
plt.scatter(z3,Pc_3,s=180,c='y')


plt.scatter(z5,Pa_5,s=60, label = '1250 = Satoms ')
plt.scatter(z5,Pb_5,s=100,c='c', label = '2000')
plt.scatter(z5,Pc_5,s=180,c='y', label = '5000')


plt.scatter(z,P15, s = 400, c='r', label = '15000')
plt.xlabel('PinHole Diameter', fontsize=22)
plt.ylabel('Total % Captured',fontsize=22)
plt.title('Unification for PinHole - Natoms const. @ 500')
plt.legend(loc=2)
plt.show()