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
L_lower_bound = 5
L_upper_bound = 20
N_L = 15
L_list = np.linspace(L_lower_bound,L_upper_bound, N_L)
k_pL = [k_p*l for l in L_list]

g_RH_list = np.linspace(1,5,5)

all_W_S_over_L_list = []
all_g_list = []
all_L_list = []

all_W_S_over_SA_list = []
all_W_S_over_R_list = []

for g_RH in g_RH_list:
	g_W_S_over_L_list = []
	g_W_S_over_SA_list = []
	g_W_S_over_R_list = []

	for L in L_list:
		H = (3/(pi*g_RH**2))**(1/3)*L
		R = H*g_RH

		SA=pi*R*(R+(H**2+R**2)**(1/2))

		n = 2
		kinks = 0

		[r1,r2,r3] = [0,0,0]
		[z1,z2,z3] = [0,0,0]
		[r_1_extra_coeffs,z_1_extra_coeffs,r_2_extra_coeffs,z_2_extra_coeffs,r_3_extra_coeffs,z_3_extra_coeffs,r_4_extra_coeffs,z_4_extra_coeffs] = [[0],[0],[0],[0],[0],[0],[0],[0]]

		from make_all_coeffs import make_all_coeffs
		r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs = make_all_coeffs(kinks,n,R,H,r1,r2,r3,z1,z2,z3,r_1_extra_coeffs,z_1_extra_coeffs,r_2_extra_coeffs,z_2_extra_coeffs,r_3_extra_coeffs,z_3_extra_coeffs,r_4_extra_coeffs,z_4_extra_coeffs)

############# GET CAPTURE WIDTH ##################
		W_S = gen_control_function(make_graph,H_s,T_p,N_beta,Spectrum,Spectrum_star,L_panel_factor,L,n,kinks,R,H,r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs,frequencies,calc_res,inc_L,r1,r2,r3,z1,z2,z3)
	
		all_W_S_over_L_list.append(W_S/L)
		g_W_S_over_L_list.append(W_S/L)
		all_L_list.append(L)
		all_g_list.append(g_RH)
		all_W_S_over_SA_list.append(W_S/(SA**(1/2)))
		g_W_S_over_SA_list.append(W_S/(SA**(1/2)))
		all_W_S_over_R_list.append(W_S/R)
		g_W_S_over_R_list.append(W_S/R)  
		
	plt.plot(k_pL,g_W_S_over_R_list,label='R/H=%.1f'%(g_RH))
	
plt.legend()
plt.show()
plt.close()
sys.exit()
	

######### DETERMINING OPTIMAL  ########

import operator
ind_L, val_L = max(enumerate(all_W_S_over_L_list),key=operator.itemgetter(1))
ind_SA,val_SA = max(enumerate(all_W_S_over_SA_list),key=operator.itemgetter(1))
ind_R, val_R = max(enumerate(all_W_S_over_R_list),key=operator.itemgetter(1))

L_opt_L = all_L_list[ind_L]
g_opt_L = all_g_list[ind_L]
W_S_over_L_opt = all_W_S_over_L_list[ind_L]

L_opt_SA = all_L_list[ind_SA]
g_opt_SA = all_g_list[ind_SA]
W_S_over_SA_opt = all_W_S_over_SA_list[ind_SA]

L_opt_R = all_L_list[ind_R]
g_opt_R = all_g_list[ind_R]
W_S_over_R_opt = all_W_S_over_R_list[ind_R]

print ('W_S/L')
print (L_opt_L)
print (g_opt_L)
print (W_S_over_L_opt)
print ('W_S/(SA)^(1/2)')
print (L_opt_SA)
print (g_opt_SA)
print (W_S_over_SA_opt)
print ('W_S/R')
print (L_opt_R)
print (g_opt_R)
print (W_S_over_R_opt)

##### Change to 2D arrays for 3D ploting #######

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

k_pL_list = [k_p*l for l in L_list]

X = np.array(k_pL_list)
Y = np.array(g_RH_list)
X_graph,Y_graph = np.meshgrid(X,Y)
Z_L = np.array(all_W_S_over_L_list).reshape(len(g_RH_list),len(L_list))
Z_SA = np.array(all_W_S_over_SA_list).reshape(len(g_RH_list),len(L_list))
Z_R = np.array(all_W_S_over_R_list).reshape(len(g_RH_list),len(L_list))

from matplotlib import cm
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(X_graph,Y_graph,Z_L,rstride=1,cstride=1,cmap=cm.coolwarm)

ax.set_xlabel('k_p*L')
ax.set_ylabel('R/H')
ax.set_zlabel('W_S/L')

plt.savefig('figures/cone/e_and_e_W_S_over_L.png')
plt.show()
plt.close()

###
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(X_graph,Y_graph,Z_SA,rstride=1,cstride=1,cmap=cm.coolwarm)

ax.set_xlabel('k_p*L')
ax.set_ylabel('R/H')
ax.set_zlabel('W_S/SA^(1/2)')

plt.savefig('figures/cone/e_and_e_W_S_over_SA.png')
plt.show()
plt.close()

###
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(X_graph,Y_graph,Z_R,rstride=1,cstride=1,cmap=cm.coolwarm)

ax.set_xlabel('k_p*L')
ax.set_ylabel('R/H')
ax.set_zlabel('W_S/R')

plt.savefig('figures/cone/e_and_e_W_S_over_R.png')
plt.show()
plt.close()

#### WRITE TO FILE #####
np.savetxt('matrices_etc/cone/kpL.out',X,delimiter=',')
np.savetxt('matrices_etc/cone/g.out',Y,delimiter=',')
np.savetxt('matrices_etc/cone/W_S_over_L.out',Z_L,delimiter=',')
np.savetxt('matrices_etc/cone/W_S_over_SA.out',Z_SA,delimiter=',')
np.savetxt('matrices_etc/cone/W_S_over_R.out',Z_R,delimiter=',')







