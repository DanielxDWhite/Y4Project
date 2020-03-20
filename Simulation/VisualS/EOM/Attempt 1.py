'''
        Curve Fit Attempt 1
'''

import numpy as numpy
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

a = 2
c = 0
def fun(x, a ,c):
    return a*x + c

x = [-2,-1,0,  1,2]
y = [-4,-2,0.1,2,4]
f = [-4,-2,0,  2,4]

b = curve_fit(fun(x, a, c), x,f)

print(b[0])