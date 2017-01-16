from __future__ import division
def make_gdf(R,H,r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs,n,make_graph,N_t,N_s_list):

    ### Make a GDF for a particular geometry 
    ### General geometry, described by Chevyshev polynomials

    import numpy as np
    import math
    import datetime

    from new_Cheby import define_Chebyshev

    ### Header
    Header = datetime.date.today().strftime("%B %d, %Y")
    ### Characteristic length scale, gravitational constant
    L_1 = '1.0 9.81'
    ### Are x=0 and y=0 planes of symmetry?
    L_2 = '0 0'

    kinks = len(N_s_list)-1
    N_s = sum(N_s_list)-kinks

    NNODES = N_t*N_s
    N_Panel = N_t*(N_s-1)

    ### Array of node indices: cylindrical coordinates and Cartesian coordinates
    node_ind_cyl = [[0 for j in range(3)] for i in range(NNODES)] # theta r z
    node_ind_car = [[0 for j in range(3)] for i in range(NNODES)] # x y z

    ### here we must describe the geometry of the body using Chebyshev poly's
    delta_t = np.linspace(0,(2*3.14*(1-1/N_t)),N_t).tolist()# 2*3.14/N_t
    delta_s = np.linspace(0,1,N_s).tolist()
    r_of_s_1, z_of_s_1 = define_Chebyshev(r_1_coeffs,z_1_coeffs,n,N_s_list[0])
    r_of_s = r_of_s_1
    z_of_s = z_of_s_1
    if kinks > 0:
        r_of_s_2, z_of_s_2 = define_Chebyshev(r_2_coeffs,z_2_coeffs,n,N_s_list[1])
        r_of_s = r_of_s + r_of_s_2[1:]
        z_of_s = z_of_s + z_of_s_2[1:]
        if kinks > 1:
            r_of_s_3, z_of_s_3 = define_Chebyshev(r_3_coeffs,z_3_coeffs,n,N_s_list[2])
            r_of_s = r_of_s + r_of_s_3[1:]
            z_of_s = z_of_s + z_of_s_3[1:]
            if kinks > 2:
               r_of_s_4, z_of_s_4 = define_Chebyshev(r_4_coeffs,z_4_coeffs,n,N_s_list[3])
               r_of_s = r_of_s + r_of_s_4[1:]
               z_of_s = z_of_s + z_of_s_4[1:] 

    for i in range(NNODES):
        theta = delta_t[i%N_t]
        z = z_of_s[int(math.ceil((i+1)/N_t)-1)]
        r = r_of_s[int(math.ceil((i+1)/N_t)-1)]
        node_ind_cyl[i][0] = theta
        node_ind_cyl[i][1] = r
        node_ind_cyl[i][2] = z

    for i in range(NNODES):
        node_ind_car[i][0] = node_ind_cyl[i][1]*math.cos(node_ind_cyl[i][0])
        node_ind_car[i][1] = node_ind_cyl[i][1]*math.sin(node_ind_cyl[i][0])
        node_ind_car[i][2] = node_ind_cyl[i][2] 

    ### Warnings about parametric troubles

    how_many = 1000
    how_close = 0.01
    r_of_s_more = []
    z_of_s_more = []

    r_of_s_1_m, z_of_s_1_m = define_Chebyshev(r_1_coeffs,z_1_coeffs,n,how_many)
    r_of_s_m = r_of_s_1_m
    z_of_s_m = z_of_s_1_m
    if kinks > 0:
        r_of_s_2_m, z_of_s_2_m = define_Chebyshev(r_2_coeffs,z_2_coeffs,n,how_many)
        r_of_s_m = r_of_s_m + r_of_s_2_m
        z_of_s_m = z_of_s_m + z_of_s_2_m
        if kinks > 1:
            r_of_s_3_m, z_of_s_3_m = define_Chebyshev(r_3_coeffs,z_3_coeffs,n,how_many)
            r_of_s_m = r_of_s_m + r_of_s_3_m
            z_of_s_m = z_of_s_m + z_of_s_3_m
            if kinks > 2:
               r_of_s_4_m, z_of_s_4_m = define_Chebyshev(r_4_coeffs,z_4_coeffs,n,how_many)
               r_of_s_m = r_of_s_m + r_of_s_4_m
               z_of_s_m = z_of_s_m + z_of_s_4_m 
    warnings = []

    for i in range(N_s):
        if r_of_s_m[i] < 0 and i%N_s != 0:
            #print 'r < 0 at s = ',delta_s[i]
            warnings.append('r < 0')
        if z_of_s_m[i] > 0 and i%N_s != 0:
            #print 'z > 0 at s = ',delta_s[i]
            warnings.append('z > 0 at s = %s, z = %s'%(delta_s[i],z_of_s_m[i]))

    for i in range(len(r_of_s_more)):
        r_1 = r_of_s_more[i]
        r_1_high = float("{0:.5f}".format(r_1+how_close))
        r_1_low = float("{0:.5f}".format(r_1-how_close))
        z_1 = float("{0:.5f}".format(z_of_s_more[i]))
        z_1_low = float("{0:.5f}".format(z_1-how_close))
        z_1_high = float("{0:.5f}".format(z_1+how_close))
        for j in range(how_many):
            #s_2 = float("{0:.5f}".format(delta_s_more[j]))
            r_2 = float("{0:.5f}".format(r_of_s_more[j]))
            z_2 = float("{0:.5f}".format(z_of_s_more[j]))
            if j < (i-10) or j > (i+10):
                if r_1_low < r_2 < r_1_high:
                #print 'R crosses at s = ', delta_s_more[i],' ',delta_s_more[j]
                    if z_1_low < z_2 < z_1_high:
                        warnings.append('crosses')
                        #print 'Z crosses at s = ', delta_s_more[i],' ',delta_s_more[j]

    ### Now must turn these points into panels
    ### Define new list
    panel_list = [[0 for j in range(4)]for i in range(N_Panel)] # 1 2 3 4
    for i in range(N_Panel):
        if i < (N_t*N_s-N_t):
            panel_list[i][0] = (N_t)*(math.ceil((i+1)/N_t)-1)+(i+1)%N_t+1
            panel_list[i][1] = (N_t)*(math.ceil((i+1)/N_t)-1)+i%N_t+1
            panel_list[i][2] = (N_t)*(math.ceil((i+1)/N_t))+i%N_t+1
            panel_list[i][3] = (N_t)*(math.ceil((i+1)/N_t))+(i+1)%N_t+1
        elif i < N_Panel-N_t:
            panel_list[i][0] = (math.floor((i-(N_s*N_t-N_t))/N_t+1))*N_t + (i+1)%N_t + (N_s*N_t-N_t)+1
            panel_list[i][1] = (math.floor((i-(N_s*N_t-N_t))/N_t+1))*N_t + i%N_t+ (N_t*N_s-N_t)+1
            if i < (N_s-1)*N_t+N_t:
                panel_list[i][2] = i%N_t + 1
                panel_list[i][3] = (i+1)%N_t + 1
            else:
                panel_list[i][2] = (N_t*N_s-N_t) + N_t*(math.floor((i-(N_s*N_t-N_t))/N_t)) + i%N_t+1
                panel_list[i][3] = N_s*N_t + 2 + (math.floor((i-N_s*N_t)/N_t))*N_t + (i+1)%N_t -1
        else:
            panel_list[i][0] = N_Panel+1
            panel_list[i][1] = N_Panel+1
            panel_list[i][2] = N_Panel - N_t + i%N_t+1
            panel_list[i][3] = N_Panel - N_t + (i+1)%N_t+1

    ### Now we have the indices for panels. Must now import the correct string of coordinates for each node for each panel
    Panel_coord_str = ''
    X_coord_list = [[0 for j in range(3)]for i in range(N_Panel*4)]
    Y_coord_list = [[0 for j in range(3)]for i in range(N_Panel*4)]
    Z_coord_list = [[0 for j in range(3)]for i in range(N_Panel*4)]
    for i in range(N_Panel):
        str_1 = ''
        for j in range(4):
            str_2 = ''
            for k in range(3):
                str_2 = str_2 + str(node_ind_car[int(panel_list[i][j]-1)][k]) + ' '
            str_1 = str_1 + str_2 + '\n'
            X_coord_list[4*i+j] = float(node_ind_car[int(panel_list[i][j]-1)][0])
            Y_coord_list[4*i+j] = float(node_ind_car[int(panel_list[i][j]-1)][1])
            Z_coord_list[4*i+j] = float(node_ind_car[int(panel_list[i][j]-1)][2])
        Panel_coord_str = Panel_coord_str + str_1

    gdf_string = "\n".join([Header,L_1,L_2,str(N_Panel),Panel_coord_str])

    ### Now we want to graph the body with panels
    if make_graph == 'yes':
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import axes3d, Axes3D
        fig = plt.figure()
        ax = Axes3D(fig)
    # plot nodes
        ax.scatter(X_coord_list,Y_coord_list,Z_coord_list)

        max_r = max(max(np.absolute(X_coord_list)),max(np.absolute(Y_coord_list)))
        max_z = max(np.absolute(Z_coord_list))

    # plot function at each theta and lines from origin for panels on top
        for i in range(N_t):
            theta = delta_t[i]
            z = z_of_s
            r = r_of_s
            x = [r_c*math.cos(theta) for r_c in r]
            y = [r_c*math.sin(theta) for r_c in r]
            ax.plot(x,y,z,'b')

    # plot at each z for function
        delta_t_for_plot = delta_t + [2*3.14]
        for i in range(N_s):
            z = z_of_s[i]
            r = r_of_s[i]
            x = [r*math.cos(theta_c) for theta_c in delta_t_for_plot]
            y = [r*math.sin(theta_c) for theta_c in delta_t_for_plot]
            z_r = [z]*(N_t+1)
            ax.plot(x,y,z_r,'b')
    
    # set limits
        max_r_z = max(max_r,max_z)
        ax.set_xlim3d([-1.25*max_r_z,1.25*max_r_z])
        ax.set_xlabel('r')
        ax.set_ylim3d([-1.25*max_r_z,1.25*max_r_z])
        ax.set_zlabel('z')
        ax.set_zlim3d([-1.25*max_r_z,1.25*max_r_z])

    # draw axis
        ax.plot([-1.5*max_r_z,1.5*max_r_z],[0,0],[0,0],'k')
        ax.plot([0,0],[-1.5*max_r_z,1.5*max_r_z],[0,0],'k')
        ax.plot([0,0],[0,0],[-1.5*max_r_z,0.5*max_r_z],'k')

        ax.view_init(elev=60., azim=160)

        plt.show()
        plt.close()

    # Draw 2D version to see shape of body
        #fig2 = plt.figure()
        theta = 0
        z = z_of_s
        r = r_of_s
        #x = [r_c*math.cos(theta) for r_c in r]
        label_string = 'n = %s'%(str(n))
        plt.plot(r,z,label=label_string)
        plt.plot([-1.5*max_r_z,1.5*max_r_z],[0,0],'k')
        plt.plot([0,0],[-1.5*max_r_z,0.5*max_r_z],'k')
        plt.xlim([-1.25*max_r_z,1.25*max_r_z])
        plt.ylim([-1.25*max_r_z,0.25*max_r_z])
        plt.xlabel('r')
        plt.ylabel('z')
        plt.show()
        plt.close()
    

    return gdf_string,z_of_s,r_of_s,warnings




#R = 1
#H = 2
#n = 2
#r_1_coeffs = [0,1,-0.5]
#r_2_coeffs = [0,2,0]
#r_3_coeffs = [0,0.5,0.25]
#r_4_coeffs = [0,0,0]
#z_1_coeffs = [0,-3,0.5]
#z_2_coeffs = [0,-4,0.5]
#z_3_coeffs = [0,-1,0.01]
#z_4_coeffs = [0,0,-0.01]
#no_s = 5
#N_r = 0
#N_t = 5
#H0 = 0
#kinks = 0
#make_graph = 'yes'

#test = make_gdf(R,H,r_1_coeffs,r_2_coeffs,r_3_coeffs,r_4_coeffs,z_1_coeffs,z_2_coeffs,z_3_coeffs,z_4_coeffs,s_1,s_2,s_3,N_t,no_s,N_r,H0,n,make_graph,kinks)
