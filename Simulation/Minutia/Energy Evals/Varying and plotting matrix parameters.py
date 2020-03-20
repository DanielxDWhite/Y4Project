import matplotlib.pyplot as plt
import numpy as np


#x = [0,1,2,3,4,5,6,7]



def H7(x):
    H7a = np.array([[x,0],
                   [0,x**2]])
    return np.linalg.eigvals(H7a)
    
print(H7(10))


#print(H7(x))

E1 = np.zeros((1,0))
E2 = np.zeros((1,0))
g=10
t = (0,1,2,3,4,5,6,7,8,9)

for i in range(0,g):
    i1, i2 =  np.split(H7(i),2) 
    E1 = np.append(E1, i1)
    E2 = np.append(E2, i2)
print(E1)
print(E2)

plt.plot(t, E1)
plt.plot(t, E2)
plt.show()



'''
print(H7(1))
E1, E2 = np.split(H7(i), 2)
print(E1)


#,print(H7())
'''