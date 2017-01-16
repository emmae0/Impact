from __future__ import division
from Cheby_with_lengths import define_Cheby_lengths

def Volume_fun(n,kinks,r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs):

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

    return Volume
