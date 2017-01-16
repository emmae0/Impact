from __future__ import division
def get_res_values(omega,gdf_string,L,A):
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
    f_dam = open('general.1','r+')
    g_dam = list(f_dam)
    A22_W = float(g_dam[1][27:40])
    A22 = A22_W/(L**3)
    B22_W = float(g_dam[1][42:54])
    B22 = B22_W*omega/((L**3)*((9.81/L)**(1/2)))
    h_dam = open('general.2','r+')
    j_dam = list(h_dam)
    #X2_W = float(j_dam[1][63:77])
    X2_W = float(j_dam[1][34:49])
    X2 = X2_W/(A*(L**2))
    k_dam = open('general.4','r+')
    l_dam = list(k_dam)
    chsi = float(l_dam[1][34:49])
    return A22_W,A22,B22_W,B22,X2_W,X2,chsi
