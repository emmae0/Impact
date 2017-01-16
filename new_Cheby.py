from __future__ import division
import numpy as np
import math
import scipy.integrate
def define_Chebyshev(r_coeffs,z_coeffs,n,no_s):
    s_list = np.linspace(0,1,no_s).tolist()
    def pick_Chebyshev(n):
        
        import functools

        T_0 = lambda x: 1
        T_1 = lambda x: x

        U_0 = lambda x: 1
        U_1 = lambda x: 2*x

        T_funcs = []
        U_funcs = []

        T_funcs.append(T_0)
        U_funcs.append(U_0)
        T_funcs.append(T_1)
        U_funcs.append(U_1)

        for i in range(2,int(n+1)):
            T_funcs.append(lambda x, i=i: (2*x*T_funcs[i-1](x)-T_funcs[i-2](x)))

        for i in range(2,int(n+1)):
            U_funcs.append(lambda x, i=i: x*U_funcs[i-1](x)+T_funcs[i](x))

        return T_funcs,U_funcs

    T_functions,U_functions = pick_Chebyshev(n)

    T_func_prime = [lambda x: 0]
    for i in range(1,int(n+1)):
        T_func_prime.append(lambda x, i=i: i*U_functions[i-1](x))

    def func_output(n,coeffs,T_functions,s):
        x_is = 2*s-1
        rest_of_sum = 0
        for i in range(0,int(n+1)):
            rest_of_sum = rest_of_sum + coeffs[i]*T_functions[i](float(x_is))
        return rest_of_sum

    z_out = [0 for i in range(no_s)]
    r_out = [0 for i in range(no_s)]
    for i in range(no_s):
        z_out[i] = func_output(n,z_coeffs,T_functions,s_list[i])
        r_out[i] = func_output(n,r_coeffs,T_functions,s_list[i])

    return r_out,z_out

## test

#n = 1
#r_coeffs = [5,-5]
#z_coeffs = [-5,-5]
#no_s = 5

#r_out_mat, z_out_mat = define_Chebyshev(r_coeffs,z_coeffs,n,no_s)

#print r_out_mat
#print z_out_mat
