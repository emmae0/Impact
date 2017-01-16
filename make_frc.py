from __future__ import division
import numpy as np
def make_frc(H,R):

    import datetime

    ### Header - today's date
    Header = datetime.date.today().strftime("%B %d, %Y")

    ### option indices
    L_1 = '1  1  1  1  0  0  0  1  0 0 0 0 0 0 0 0'

    ### VCG
    #L_2 = '0 0 %s'%str(-H/2)
    L_2 = '0'

    ### Matrix of gyration
    #L_3 = '%s 0.0 0.0'%str(((9*R**2+12*H**2)**(1/2))/6)
    #L_4 = '0.0  %s  0.0'%str(((9*R**2+12*H**2)**(1/2))/6)
    #L_5 = '0.0 0.0 %s'%str(2**(1/2)*R/2)

    L_3 = '1 0 0'
    L_4 = '0 1 0'
    L_5 = '0 0 1'

    ###NBETAH
    L_6 = '0'

    ### NFIELD for external flow (?)
    L_7 = '0'

    frc_string = "\n".join([Header,L_1,L_2,L_3,L_4,L_5,L_6,L_7])

    return frc_string
