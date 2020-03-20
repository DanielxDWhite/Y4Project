import matplotlib.pyplot as plt
import numpy as np

t = np.linspace(0,10,10)

''' Atom & Light Parameters'''
c = 3#*10**8
f = 100#*10**15             # Frequency
y = 11                        # Gamma
m =1
Ir = 2
d = 0.1
                  

''' Initial Conditions'''
v0 = -1
x0  = 0

def x(t, F):
    return F*t**2/(2*m)+v0*t + x0

def v(t, F):
    return F*t/m+v0

def F(Ir, v(t,F(Ir, v(t,F(Ir, v(t,F(Ir, v(t,...   ), d, f, y) ), d, f, y) ), d, f, y) ), d, f, y): ...  # ):
    DMaa = Ir/(1+Ir+4*(d-v*f/c)/y**2)
    return DMaa*y*2 #hbaar k to come

x = [x0]
v = [v0]
F0 = F(Ir, v0, d, f, y)
F = [F0]   
    
for i in range(0,len(t)-1):
    
    Fi = np.append(F, F(Ir, v(i, F(Ir, v(i, F), d, f, y))))
    F = Fi
                #FUNCEPTION
    
    xi = np.append(x, x(F))
    x=xi
    
    vi = np.append(v, v(i,F))
    v=vi                           
    print(v)#Assayag's Equality
    


print(x(10,F(1,2,3,f,y)))