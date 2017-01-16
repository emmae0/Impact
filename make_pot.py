from __future__ import division
### Make POT file for WAMIT with specific frequencies and wave directions
import numpy as np
#def make_pot(N_T,T_min,N_beta,T_max,R,spacing,specific_kL,L):
def make_pot(specific_k,N_beta):

    import datetime

    N_T = len(specific_k)

    ### Header - today's date
    H = datetime.date.today().strftime("%B %d, %Y")
    ### Infinite water dept
    L_1 = '-1 0. 0. 0. 0.'
    ### solve radiation and diffraction only for specified mode
    L_2 = '0 0 IRR=0'
    ### IMODE
    L_3 = '0 0 1 0 0 0'
    ### Number of wave periods
    L_4 = str(N_T)
    ### Number of wave headings
    L_5 = str(N_beta)

#    delta_T = (T_max-T_min)/N_T

#    ka_min = R*4*3.14**2/(9.81*T_max**2)
#    ka_max = R*4*3.14**2/(9.81*T_min**2)

    T_list_specific = [2*3.14/(k*9.81)**(1/2) for k in specific_k]
    flip_T_specific = T_list_specific[::-1]

#    ka_list = np.linspace(ka_min,ka_max,N_T)

#    T_list = [R**(1/2)*2*3.14/(9.81*i)**(1/2) for i in ka_list]
#    flip_T = T_list[::-1]

    Per_str = '%.10f' % min(T_list_specific) + '\n'
    for i in range(1,N_T):
        Per_str = Per_str + '%.10f' % + (flip_T_specific[i]) + '\n' 
#        elif spacing == 'T':
#            Per_str = Per_str + '%.10f' % + (T_min+(i+1)*delta_T) + '\n'
#        elif spacing == 'kL':
#            Per_str = Per_str + '%.10f' % + (flip_T[i]) + '\n'
    Per_str = Per_str[:-1]


    Beta_str = ''
    for i in range(N_beta):
        Beta_str = Beta_str + str((i)*(360/N_beta)) + '\n'

    pot_string =  "\n".join([H,L_1,L_2,L_3,L_4,Per_str,L_5,Beta_str])
    
    return pot_string

