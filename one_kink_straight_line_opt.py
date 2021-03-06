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
sys.path.append('/usr/local/lib/python2.7/site-packages')
import nlopt
from gen_control_function import gen_control_function
import scipy.optimize
make_graph = 'yes'
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

L_panel_factor = 10

### size/shape of body
def power_opt(variables,grad):

    L = variables[0]
    g_RH = variables[1]
    g_1 = variables[2]
    kappa = variables[3]
    
    def to_get_r1(r1):
        R = kappa*r1
        H = R/g_RH
        z1 = r1/g_1

        n = 2
        kinks = 1

        [r1,r2,r3] = [r1,0,0]
        [z1,z2,z3] = [-z1,0,0]
        [r_1_extra_coeffs,z_1_extra_coeffs,r_2_extra_coeffs,z_2_extra_coeffs,r_3_extra_coeffs,z_3_extra_coeffs,r_4_extra_coeffs,z_4_extra_coeffs] = [[0],[0],[0],[0],[0],[0],[0],[0]]

        from make_all_coeffs import make_all_coeffs
        r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs = make_all_coeffs(kinks,n,R,H,r1,r2,r3,z1,z2,z3,r_1_extra_coeffs,z_1_extra_coeffs,r_2_extra_coeffs,z_2_extra_coeffs,r_3_extra_coeffs,z_3_extra_coeffs,r_4_extra_coeffs,z_4_extra_coeffs)

        from Volume import Volume_fun
        Vol_from_fun = Volume_fun(n,kinks,r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs)

        return Vol_from_fun-L**3

    initial_value = L/2
    r1 = scipy.optimize.newton(to_get_r1,initial_value, tol=1e-2, maxiter=100)

    #### Got z1--now run.

    R = kappa*r1
    H = R/g_RH
    z1 = r1/g_1

    n = 2
    kinks = 1

    [r1,r2,r3] = [r1,0,0]
    [z1,z2,z3] = [-z1,0,0]
    [r_1_extra_coeffs,z_1_extra_coeffs,r_2_extra_coeffs,z_2_extra_coeffs,r_3_extra_coeffs,z_3_extra_coeffs,r_4_extra_coeffs,z_4_extra_coeffs] = [[0],[0],[0],[0],[0],[0],[0],[0]]

    from make_all_coeffs import make_all_coeffs
    r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs = make_all_coeffs(kinks,n,R,H,r1,r2,r3,z1,z2,z3,r_1_extra_coeffs,z_1_extra_coeffs,r_2_extra_coeffs,z_2_extra_coeffs,r_3_extra_coeffs,z_3_extra_coeffs,r_4_extra_coeffs,z_4_extra_coeffs)


############# GET CAPTURE WIDTH ##################

    W_S_over_L = gen_control_function(make_graph,H_s,T_p,N_beta,Spectrum,Spectrum_star,L_panel_factor,L,n,kinks,R,H,r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs,frequencies,calc_res,inc_L,r1,r2,r3,z1,z2,z3)

    sys.exit()

    return W_S_over_L

########### OPTIMIZATION ##################

lower_bounds = [1, 0.5, 0.5, 0.5] #  L R/H r1/z1 R/r1
upper_bounds = [40, 5, 5, 5]
initial_vals = [37.5, 2.1, 1.2, 1.4] #[20, 1, 1, 1]
#initial_vals = [(upper_bounds[i]-lower_bounds[i])/2 for i in range(len(lower_bounds))]

def bobyqa_opt():
    opt_bobyqa = nlopt.opt(nlopt.LN_BOBYQA,len(lower_bounds))
    opt_bobyqa.set_lower_bounds(lower_bounds)
    opt_bobyqa.set_upper_bounds(upper_bounds)
    opt_bobyqa.set_max_objective(power_opt)
    opt_bobyqa.set_xtol_rel(0.5)
    x_bobyqa = opt_bobyqa.optimize(initial_vals)
    maxf_bobyqa = opt_bobyqa.last_optimum_value()
    return x_bobyqa, maxf_bobyqa

[x, ma] = bobyqa_opt()

print 'L', x[0]
print 'R/H', x[1]
print 'r1/z1', x[2]
print 'R/r1', x[3]
print 'W_S/L, L', ma





