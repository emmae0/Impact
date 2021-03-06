from __future__ import division
import scipy.integrate
import numpy as np
import matplotlib.pyplot as plt
import sys
from heapq import nsmallest

###### this is JUST finding W_S^, which is the max considering beta to be any value

def power_calculations(N_T, A22, B22, X2, frequencies, mass, S,Spectrum,Spectrum_star,rho,g,L,k_p,H_s):

    omega_p = (k_p*g)**(1/2)
    T_p = 2*3.14/omega_p

    power_warnings = []

    freq_star = [o/omega_p for o in frequencies]
    k = [o**2/g for o in frequencies]

    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return idx

    Spec_dim = [Spectrum(o) for o in frequencies]
    Spec_star = [Spectrum_star(o) for o in freq_star]

####### COMPARING TO INFLUX #########
    
    ### P_F = amount of energy available per unit crest length = 1/2*rho*g*cg
    P_F = [1/4*rho*g**2/o for o in frequencies]
    P_F_star = [p/(rho*g**(1/2)*L**(3/2)*omega_p**2) for p in P_F]

    P_F_star_2 = [1/4*g**(3/2)*L**(-3/2)*(1/omega_p**3)*(1/o) for o in freq_star]
    ### quick non-dim check
    for i in range(N_T):
        dif = P_F_star[i]-P_F_star_2[i]
        if round(dif,5) != 0:
            print ('error non-dim 1', dif)
            sys.exit()


############# to calculate P_OPT ###############

    B_function = [0 for i in range(N_T)]
    for i in range(N_T):
        B_function[i] = (B22[i]**2+(g**(1/2)*L**(-1/2)*S/omega_p-L**(1/2)*g**(-1/2)*(mass+A22[i])*freq_star[i]**2*omega_p)**2/(freq_star[i]**2*omega_p**2))**(1/2)

    beta_min = 0.1*B_function[50]
    beta_max = 3*B_function[50]

    betas = np.linspace(beta_min,beta_max,500) 
    
    Script_P_star_betas = [0 for p in range(len(betas))]
    
    for j in range(len(betas)):

        P_star_betas = [0 for q in range(N_T)]

        for i in range(N_T):

            ### fixed at beta[j]
            P_star_betas[i] = 1/2*betas[j]*X2[i]**2*freq_star[i]**2/((S-(mass+A22[i])*freq_star[i]**2*omega_p**2*L/g)**2+freq_star[i]**2*omega_p**2*L/g*(betas[j]+B22[i])**2)

        P_star_betas_times_spec_star = [P_star_betas[i]*Spec_star[i] for i in range(N_T)]
        Script_P_star_betas[j] = 2*scipy.integrate.trapz(P_star_betas_times_spec_star[::-1],freq_star[::-1])

    import operator
    max_ind,Script_P_star_OPT = max(enumerate(Script_P_star_betas), key=operator.itemgetter(1))

    ### compare to P_F

    Spec_star_times_inv_freq_star = [Spec_star[i]*(1/freq_star[i]) for i in range(N_T)]

    Script_P_F_star = 1/2*scipy.integrate.trapz(Spec_star_times_inv_freq_star[::-1],freq_star[::-1])*g**(3/2)*(1/omega_p)**3*L**(-3/2)

    ### checking this is right
    Spec_dim_times_P_F_dim = [Spec_dim[i]*P_F[i] for i in range(N_T)]
    P_F_dim_1 = 2*scipy.integrate.trapz(Spec_dim_times_P_F_dim[::-1],frequencies[::-1])
    P_F_dim_2 = Script_P_F_star*H_s**2*rho*g**(1/2)*omega_p**2*L**(3/2) 

    ### quick non-dim check
    for i in range(N_T):
        dif = P_F_dim_1 - P_F_dim_2
        if round(dif,5) != 0:
            print ('error non-dim 3', dif)
            sys.exit()

    W_S_OPT = L*Script_P_star_OPT/Script_P_F_star

    ### quick non-dim check

    Script_P_dim_OPT = Script_P_star_OPT*H_s**2*rho*g**(1/2)*L**(5/2)*omega_p**2
    Script_P_dim_F = Script_P_F_star*H_s**2*rho*g**(1/2)*omega_p**2*L**(3/2)

    W_S_dim = Script_P_dim_OPT/(Script_P_dim_F)

    dif = W_S_dim - W_S_OPT
    if round(dif,5) != 0:
        print ('error non-dim 4', dif)
        sys.exit()

    return W_S_OPT



