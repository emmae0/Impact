from __future__ import division
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.integrate
import sys
import decimal
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import operator
import os
from gen_control_function import gen_control_function
make_graph = 'no'
calc_res = 'no'
inc_L = 'yes'

### constants
rho = 1000
g = 9.81
pi = 3.14

### wave/sea conditions
H_s = 1
T_p = 10
k_p = 4*pi**2/(g*T_p**2)
lambda_p = 2*pi/k_p
omega_p = 2*pi/T_p
N_beta = 1
frequencies = np.linspace(0.08,0.8,99).tolist()

Spectrum = lambda w: (1.25/4)*(omega_p**4/w**5)*(H_s**2)*math.exp((-1.25*((omega_p/w)**4)))
Spectrum_star = lambda w: (1.25/4)*(1/(w**5))*math.exp(-1.25*(1/w)**4)

L_panel_factor = 5
L_lower_bound = 20
L_upper_bound = 30
N_L = 15
L_list = np.linspace(L_lower_bound,L_upper_bound, N_L)
k_pL = [k_p*l for l in L_list]

W_S_over_L_list = []
W_S_over_SA_list = []
W_S_over_R_list = []

### size/shape of body
for L in L_list:

    R = L/((2/3*pi)**(1/3))
    H = R

    SA = 2*pi*(R**2)

    n = 5
    kinks = 0

    [r1,r2,r3] = [R,0,0]
    [z1,z2,z3] = [-H,0,0]
    [r_1_extra_coeffs,z_1_extra_coeffs,r_2_extra_coeffs,z_2_extra_coeffs,r_3_extra_coeffs,z_3_extra_coeffs,r_4_extra_coeffs,z_4_extra_coeffs] = [[0],[0],[0],[0],[0],[0],[0],[0]]

    r_1_extra_coeffs = [0,0,0,0]
    z_1_extra_coeffs = [0,0,0,0]

    from make_all_coeffs import make_all_coeffs
    r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs = make_all_coeffs(kinks,n,R,H,r1,r2,r3,z1,z2,z3,r_1_extra_coeffs,z_1_extra_coeffs,r_2_extra_coeffs,z_2_extra_coeffs,r_3_extra_coeffs,z_3_extra_coeffs,r_4_extra_coeffs,z_4_extra_coeffs)


### to define hemisphere

    from get_Cheby_coeff import get_Cheby_coeff
    r_function = lambda x: -R*math.cos(3.14/2*(x-1))
    z_function = lambda x: R*math.sin(3.14/2*(x-1))
    r_1_coeffs = [get_Cheby_coeff(i,r_function) for i in range(n+1)]
    z_1_coeffs = [get_Cheby_coeff(i,z_function) for i in range(n+1)]

############# GET CAPTURE WIDTH ##################

    W_S = gen_control_function(make_graph,H_s,T_p,N_beta,Spectrum,Spectrum_star,L_panel_factor,L,n,kinks,R,H,r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs,frequencies,calc_res,inc_L,r1,r2,r3,z1,z2,z3)
            
    W_S_over_L_list.append(W_S/L)
    W_S_over_SA_list.append(W_S/(SA**(1/2)))
    W_S_over_R_list.append(W_S/R)


######### DETERMINING OPTIMAL  ########

import operator
ind_L, val_L = max(enumerate(W_S_over_L_list),key=operator.itemgetter(1))
ind_SA,val_SA = max(enumerate(W_S_over_SA_list),key=operator.itemgetter(1))
ind_R, val_R = max(enumerate(W_S_over_R_list),key=operator.itemgetter(1))

L_opt_L = L_list[ind_L]
W_S_over_L_opt = W_S_over_L_list[ind_L]

L_opt_SA = L_list[ind_SA]
W_S_over_SA_opt = W_S_over_SA_list[ind_SA]

L_opt_R = L_list[ind_R]
W_S_over_R_opt = W_S_over_R_list[ind_R]


print ('L',L,'SA',SA,'W_S',W_S)

print ('W_S/L')
print (L_opt_L)
print (W_S_over_L_opt)
print ('W_S/(SA)^(1/2)')
print (L_opt_SA)
print (W_S_over_SA_opt)
print ('W_S/R')
print (L_opt_R)
print (W_S_over_R_opt)


######### PLOTTING ################

plt.plot(k_pL,W_S_over_L_list,label='W_S/L')
plt.plot(k_pL,W_S_over_SA_list,label='W_S/(SA)^(1/2)')
plt.plot(k_pL,W_S_over_R_list,label='W_S/R')
plt.legend()
plt.show()
plt.close()




