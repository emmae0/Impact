from __future__ import division
### How to calculate volume, and how to alter geom to get a particular volume
import numpy as np
import scipy.integrate

from new_Cheby import define_Chebyshev

R = 10
H = 10
n = 3
kinks = 0
no_s = 10

s_list = np.linspace(0,1,no_s).tolist()

[r1,r2,r3] = [0,0,0]
[z1,z2,z3] = [0,0,0]
[r_1_extra_coeffs,z_1_extra_coeffs,r_2_extra_coeffs,z_2_extra_coeffs,r_3_extra_coeffs,z_3_extra_coeffs,r_4_extra_coeffs,z_4_extra_coeffs] = [[0],[0],[0],[0],[0],[0],[0],[0]]
### if kinks > 0, need to specify r1 and z1 (etc). if n > 2, need to specify r_1_extra_coeffs (etc if kinks > 0)

r_1_extra_coeffs = [1,-3]
z_1_extra_coeffs = [1,1]

from make_all_coeffs import make_all_coeffs
r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs = make_all_coeffs(kinks,n,R,H,r1,r2,r3,z1,z2,z3,r_1_extra_coeffs,z_1_extra_coeffs,r_2_extra_coeffs,z_2_extra_coeffs,r_3_extra_coeffs,z_3_extra_coeffs,r_4_extra_coeffs,z_4_extra_coeffs)

### limit to no kinks right now
r_of_s,z_of_s,z_prime_of_s = define_Chebyshev(r_1_coeffs,z_1_coeffs,n,no_s)

inside_integral = [0 for i in range(no_s)]
for i in range(no_s):
    inside_integral[i] = r_of_s[i]**2*z_prime_of_s[i]

volume = abs(3.14*scipy.integrate.trapz(inside_integral,s_list))

print volume
