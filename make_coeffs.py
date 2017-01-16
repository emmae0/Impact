from __future__ import division
import numpy as np

def make_coeffs(n,r_1,r_2,z_1,z_2,r_extra_coeffs,z_extra_coeffs):

    a1 = r_2/2-r_1/2
    b1 = z_2/2-z_1/2
    
    for i in range(n-1):
        if i % 2 != 0:
            a1 += -r_extra_coeffs[i]
            b1 += -z_extra_coeffs[i]

    a0 = r_1 + a1
    b0 = z_1 + b1
    
    for i in range(n-1):
        if i % 2 == 0:
            a0 += -r_extra_coeffs[i]
            b0 += -z_extra_coeffs[i]
        elif i % 2 != 0:
            a0 += r_extra_coeffs[i]
            b0 += z_extra_coeffs[i]


    r_coeffs = []
    z_coeffs = []
    r_coeffs.append(a0)
    r_coeffs.append(a1)
    z_coeffs.append(b0)
    z_coeffs.append(b1)
    for i in range(n-1):
        r_coeffs.append(r_extra_coeffs[i])
        z_coeffs.append(z_extra_coeffs[i])

    return r_coeffs, z_coeffs

## test

#n = 2
#R = 10
#H = 10
#r_1_extra_coeffs = [0]
#z_1_extra_coeffs = [0]


#r_coeff_mat,z_coeff_mat = make_coeffs(n,R,0,0,-H,r_1_extra_coeffs,z_1_extra_coeffs)

#print r_coeff_mat
#print z_coeff_mat
    
            
