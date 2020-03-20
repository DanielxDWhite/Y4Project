print('Basic Simulataion Archaatectuure')
import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0,15, 15)

x0 = -5

x= np.array([x0])
neXt = 1
print(x)

def k(i):
    return i**4




for i in range(0,len(t)-1):
    
    y = np.append(x,x0+np.sum(k(i)))
    x=y          #Assayag's Equality
    print(y)
    
plt.plot(t,y)

