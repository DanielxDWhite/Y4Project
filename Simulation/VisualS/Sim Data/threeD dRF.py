'''
3D test
'''
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np


fig = plt.figure()
ax = fig.add_subplot(1,1,1,projection='3d')

x,y,z = axes3d.get_test_data(0.05)

print(np.shape(x),np.shape(y),np.shape(z))
ax.plot_surface(x,y,z,rstride=4,cstride=5)





plt.show()


 