from __future__ import division
from make_coeffs import make_coeffs

def make_all_coeffs(kinks,n,R,H,r1,r2,r3,z1,z2,z3,r_1_extra_coeffs,z_1_extra_coeffs,r_2_extra_coeffs,z_2_extra_coeffs,r_3_extra_coeffs,z_3_extra_coeffs,r_4_extra_coeffs,z_4_extra_coeffs):

    if kinks < 3:
        r_4_coeffs = [0 for i in range(n)]
        z_4_coeffs = [0 for i in range(n)]
        if kinks < 2:
            r_3_coeffs = [0 for i in range(n)]
            z_3_coeffs = [0 for i in range(n)]
            if kinks < 1:
                r_2_coeffs = [0 for i in range(n)]
                z_2_coeffs = [0 for i in range(n)]

        
    if kinks == 0:
        r_1_coeffs, z_1_coeffs = make_coeffs(n,R,0,0,-H,r_1_extra_coeffs,z_1_extra_coeffs)
    elif kinks == 1:
        r_1_coeffs, z_1_coeffs = make_coeffs(n,R,r1,0,z1,r_1_extra_coeffs,z_1_extra_coeffs)
        r_2_coeffs, z_2_coeffs = make_coeffs(n,r1,0,z1,-H,r_2_extra_coeffs,z_2_extra_coeffs)
    elif kinks == 2:
        r_1_coeffs, z_1_coeffs = make_coeffs(n,R,r1,0,z1,r_1_extra_coeffs,z_1_extra_coeffs)
        r_2_coeffs, z_2_coeffs = make_coeffs(n,r1,r2,z1,z2,r_2_extra_coeffs,z_2_extra_coeffs)
        r_3_coeffs, z_3_coeffs = make_coeffs(n,r2,0,z2,-H,r_3_extra_coeffs,z_3_extra_coeffs)
    elif kinks == 3:
        r_1_coeffs, z_1_coeffs = make_coeffs(n,R,r1,0,z1,r_1_extra_coeffs,z_1_extra_coeffs)
        r_2_coeffs, z_2_coeffs = make_coeffs(n,r1,r2,z1,z2,r_2_extra_coeffs,z_2_extra_coeffs)
        r_3_coeffs, z_3_coeffs = make_coeffs(n,r2,r3,z2,z3,r_3_extra_coeffs,z_3_extra_coeffs)
        r_4_coeffs, z_4_coeffs = make_coeffs(n,r3,0,z3,-H,r_4_extra_coeffs,z_4_extra_coeffs)


    return r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs
