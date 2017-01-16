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

#def gen_control_function(make_graph,H_s,T_p,N_beta,Spectrum,Spectrum_star,L_panel,L,R,H,n,kinks,r1,r2,r3,z1,z2,z3,r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs,frequencies):
def gen_control_function(make_graph,H_s,T_p,N_beta,Spectrum,Spectrum_star,L_panel_factor,L_maybe,n,kinks,R,H,r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs,frequencies,calc_res,inc_L,r1,r2,r3,z1,z2,z3):

    if inc_L == 'yes':
        L = L_maybe

### constants
    rho = 1000
    g = 9.81
    pi = 3.14

### wave/sea conditions
    k_p = 4*pi**2/(g*T_p**2)
#T_p = 2*pi/(g*k_p)**(1/2)
    lambda_p = 2*pi/k_p
    omega_p = 2*pi/T_p


### Determine arclengths and volume
    from Cheby_with_lengths import define_Cheby_lengths

    L_1,V_1 = define_Cheby_lengths(r_1_coeffs,z_1_coeffs,n)
    Arclengths = [L_1]
    Volume = V_1
    if kinks > 0:
        L_2,V_2 = define_Cheby_lengths(r_2_coeffs,z_2_coeffs,n)
        Arclengths.append(L_2)
        Volume = Volume + V_2
        if kinks > 1:
            L_3,V_3 = define_Cheby_lengths(r_3_coeffs,z_3_coeffs,n)
            Arclengths.append(L_3)
            Volume = Volume + V_3
            if kinks > 2:
                L_4,V_4 = define_Cheby_lengths(r_4_coeffs,z_4_coeffs,n)
                Arclengths.append(L_4)
                Volume = Volume + V_3

    from Volume import Volume_fun
    Volume_from_fun = Volume_fun(n,kinks,r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs)
    print ('from fun', Volume_from_fun)

    print ('from L', L**3)

#    sys.exit()

    if inc_L == 'no':
        L = Volume**(1/3)

    vol_round = len(str(int(round(Volume))))
    print (vol_round)

    if round(Volume,-(vol_round-1)) != round(L**3,-(vol_round-1)): 
        print (round(Volume,-(vol_round-1)), round(L**3,-(vol_round-1)))
        sys.exit('Volume error')

    if L < lambda_p:
        L_panel = L/L_panel_factor
    else:
        L_panel = lambda_p/L_panel_factor

    for_L_panel = [R,H,r1,r2,r3,abs(z1),abs(z2),abs(z3)]
    L_panel_min = min(x for x in for_L_panel if x > 0)

    if L_panel > L_panel_min/2:
        L_panel = L_panel_min/2

    print ('L_panel', L_panel)
        

### determining number of panels

    N_t = int(math.ceil(2*pi*R/L_panel))
    N_s = [int(math.ceil(LL/L_panel)) for LL in Arclengths]


### frequencies/periods/etc

    specific_k = [o**2/g for o in frequencies]
    N_T = len(frequencies)
    T_specific = [2*pi/(g*k)**(1/2) for k in specific_k]
    periods = T_specific[::-1]
    frequencies = [2*3.14/T for T in periods] ### This gives biggest to smallest--how it is ouput from WAMIT
    k = [(o**2)/9.81 for o in frequencies]
    freq_nd = [(o/(g/L)**(1/2)) for o in frequencies]
    kL = [L*kk for kk in k]
    
########### DELETE OLD FILES ####################
        
    if os.path.isfile('general.p2f'):
        os.remove('general.p2f')
    if os.path.isfile('general.out'):
        os.remove('general.out')
    if os.path.isfile('general.1'):
        os.remove('general.1')
    if os.path.isfile('general.2'):
        os.remove('general.2')
    if os.path.isfile('general.3'):
        os.remove('general.3')
    if os.path.isfile('general.4'):
        os.remove('general.4')
    if os.path.isfile('general.8'):
        os.remove('general.8')
    if os.path.isfile('general.gdf'):
        os.remove('general.gdf')
    if os.path.isfile('general.pot'):
        os.remove('general.pot')
    if os.path.isfile('general.frc'):
        os.remove('general.frc')


######### MAKE GDF AND POT FILES ############
    from make_gdf import make_gdf
    from make_pot import make_pot
    from make_frc import make_frc

    with open('general.gdf','w') as f:
        [gdf_string,z_of_s,r_of_s,warnings] = make_gdf(R,H,r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs,n,make_graph,N_t,N_s)
        f.write(gdf_string)

    with open('general.pot','w') as h:
        pot_string = make_pot(specific_k,N_beta)
        h.write(pot_string)

    with open('general.frc','w') as gg:
        frc_string = make_frc(H,R)
        gg.write(frc_string)

    if len(warnings) > 0:
        sys.exit('Error with function')

###################### RUN WAMIT ########################
    os.system('/impact_workspace1/emmae/geom_WEC/spring_2017/poten')
    os.system('/impact_workspace1/emmae/geom_WEC/spring_2017/force')

    #import subprocess as sub
    #sub.call('/impact_workspace1/emmae/geom_WEC/spring_2017/poten',shell=True)
    #sub.call('/impact_workspace1/emmae/geom_WEC/spring_2017/force',shell=True)


####### GET ADDED MASS AND DAMPING COEFFICIENTS AND EXCITING FORCE OUT ########
 
    f = open('general.1','r+')
    hh = list(f)
    A22_W = [0 for i in range(N_T)]
    A22 = [0 for i in range(N_T)]
    B22_W = [0 for i in range(N_T)] ### Note we will need to change this if we consider more than one mode of motion!
    B22 = [0 for i in range(N_T)]
    for i in range(N_T):
        A22_W[i] = float(hh[int((i%N_T)+1)][27:41])
        A22[i] = A22_W[i]/(L**3)
        B22_W[i] = float(hh[int((i%N_T)+1)][42:55])
        B22[i] = B22_W[i]*frequencies[i]/((L**3)*((g/L)**(1/2)))

#h = open('general.2','r+')
    h = open('general.3','r+')
    j = list(h)
    X2_W = [0 for i in range((N_T)*N_beta)] ### We will have to change the X2 matrix if there is more than one direction of propagation
    A = 1 #incoming wave amplitude (always at 1)
    X2 = [0 for i in range((N_T)*N_beta)]
    for i in range((N_T)*N_beta):
        X2_W[i] = float(j[int((i%N_T)+1)][34:49])
        #X2_W[i] = float(j[int((i%N_T)+1)][63:77])
        X2[i] = X2_W[i]/(A*(L**2))
        

    l = open('general.4','r+')
    m = list(l)
    chsi_W = [0 for i in range((N_T)*N_beta)]### must change if more than one dir
    chsi = [0 for i in range((N_T)*N_beta)]
    for i in range((N_T)*N_beta):
        chsi_W[i] = float(m[int((i%N_T)+1)][34:49])
        chsi[i] = abs(chsi_W[i])

### Let's also have a dimensional form in case

    X2_dim = [x*rho*g*L**2 for x in X2]
    A22_dim = [a*rho*L**3 for a in A22]
    B22_dim = [b*rho*L**(5/2)*g**(1/2) for b in B22]

############## GEOMETRIC CALCULATIONS ###############

### S
    S_dim = 3.14*(R**2)
    S = S_dim/(L**2)

    mass_dim = rho*Volume
    mass = mass_dim/(rho*(L**3))

############ FIND RESONANT FREQUENCY #################

    from added_mass_function import added_mass
    import scipy.optimize

    def resonance_function(x):
        AM = added_mass(x,gdf_string,L)
        return rho*g*S_dim-x**2*(mass_dim+AM*rho*L**3)

    if calc_res == 'yes':
        initial_value = ((rho*9.81*2*R)/(mass_dim+0.25*rho*3.14*R**3))**(1/2)
        resonance_value = scipy.optimize.newton(resonance_function,initial_value, tol=1e-08,maxiter=100)

        from get_res_values import get_res_values
        [A22_res_W,A22_res,B22_res_W,B22_res,X2_res_W,X2_res,chsi_res] = get_res_values(resonance_value,gdf_string,L,A)

    if calc_res == 'no':
        B22_res = 0
        resonance_value = 0

########## DETERMINE POWER ###############

    from power_calculations_short import power_calculations
    W_S = power_calculations(N_T, A22, B22, X2, frequencies, mass, S,Spectrum,Spectrum_star,rho,g,L,k_p,H_s)
    
    return W_S


