from __future__ import division
def added_mass(omega,gdf_string,L):
    import os
    period_am = 2*3.14/omega
    k_am = (omega**2)/9.81
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
    if os.path.isfile('general.pot'):
        os.remove('general.pot')
        
    from make_pot import make_pot

    with open('general.pot','w') as g:
        pot_string = make_pot([k_am],1)
        g.write(pot_string)
    with open('general.gdf','w') as h:
        h.write(gdf_string)
    os.system('/usr/home/active/emmae/geom_WEC/WAMIT/general/poten')
    os.system('/usr/home/active/emmae/geom_WEC/WAMIT/general/force')
    f_am = open('general.1','r+')
    g_am = list(f_am)
    A22_W = float(g_am[int(1)][27:40])
    A22 = A22_W/(L**3)
    return A22
