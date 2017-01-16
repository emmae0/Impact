from __future__ import division
#### using collocation method, find coefficients of the Chebyshev functions for a known function
import numpy as np
import math
import scipy.integrate as integrate

def get_Cheby_coeff(k,u):
    
    def pick_Chebyshev(n):
        T_0 = lambda x: 1
        T_1 = lambda x: x

        T_funcs = []

        T_funcs.append(T_0)
        T_funcs.append(T_1)

        for i in range(2,int(n+1)):
            T_funcs.append(lambda x, i=i: (2*x*T_funcs[i-1](x)-T_funcs[i-2](x)))
        return T_funcs

    T_functions = pick_Chebyshev(k)

    T_k = T_functions[k]

    if k == 0: c_k = 2
    elif k >= 1: c_k = 1

    first_fun = integrate.quad(lambda x: u(x)*T_k(x),0,1)

    a_k = 2/(3.14*c_k)*(integrate.quad(lambda x: (u(1/2*(x-1))*T_k(x))/(1-x**2)**(1/2),-1,1))[0]

    return a_k
