import numpy as np
import scipy.integrate as si 
import scipy.optimize as so

def integrand(t):
    return 3*np.exp((-2)/(8.314*(996.74 + (1037.1*np.exp(-0.2696*t)))))

def func(x):
    return si.quad(integrand, 0, x)[0] - 7    # func is integral minus 7

sol = so.fsolve(func, 1.0)                    # equated func to 0
print(sol)