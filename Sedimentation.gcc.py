'''
Computation of terminal fall speed
'''
#GCCif (KERNEL == 1 && LongK_Options == 2)  || (COLUMN == 1 && PROCESS != 2) 
# this block is active, if Long kernel values are explicitly computed (and not read from a file) or if sedimentation is switched on in the column model

import sys
import Misc as FK
import numpy as np
import math


#converted from Botts F77 code using Beards parametrisation
def Fallg(r):
    #r      IN: radius in m
    #winf  OUT: fall speed in m/s

    #internal computations use units cm for r and cm/s for winf
    # final return with winf converted to m/s
    
    b=np.array([-0.318657e1,0.992696,-0.153193e-2,-0.987059e-3,-0.578878e-3,0.855176e-4,-0.327815e-5])
    c=np.array([-0.500015e1,0.523778e1,-0.204914e1,0.475294,-0.542819e-1,0.238449e-2])
    eta=1.818*1e-4
    xlamb=6.62*1e-6
    rhow=1
    rhoa=1.223*1e-3
    grav=980.665
    cunh=1.257*xlamb
    t0=273.15
    sigma=76.1-0.155*(293.15-t0)
    stok=2*grav*(rhow-rhoa)/(9*eta)
    stb=32*rhoa*(rhow-rhoa)*grav/(3*eta*eta)
    phy=sigma*sigma*sigma*rhoa*rhoa/((eta**4)*grav*(rhow-rhoa))
    py=phy**(1./6)
    n=r.size
    winf=np.zeros(n)

    rr=r*1e2 # rr in cm

    for j in range(0,n):
        if (rr[j]<1e-3):
            winf[j]=stok*(rr[j]*rr[j]+cunh*rr[j])
        elif( rr[j]>(1e-3) and rr[j]<(5.35*10**(-2))):
            x=math.log(stb*rr[j]*rr[j]*rr[j])
            y=0.
            #print('x',x)
            for i in range(1,8):
                y=y+b[i-1]*(x**(i-1))
                #print('y',i,y)
            xrey=(1+cunh/rr[j])*math.exp(y)
            winf[j]=xrey*eta/(2*rhoa*rr[j])
            #print('c',xrey,winf[j])
        elif(rr[j]>5.35e-2):
            bond=grav*(rhow-rhoa)*rr[j]*rr[j]/sigma;
            if (rr[j]>0.35):
                bond=grav*(rhow-rhoa)*0.35*0.35/sigma
            x=math.log(16*bond*py/3)
            y=0.
            for i in range(1,7):
                y=y+c[i-1]*(x**(i-1))
            xrey=py*math.exp(y)
            winf[j]=xrey*eta/(2*rhoa*rr[j])
            if (rr[j]>0.35):
                winf[j]=xrey*eta/(2*rhoa*0.35)
    return winf*1e-2 # convert winf from cm/s to m/s

##############################################################################
#GCCendif /* (KERNEL == 1 && LongK_Options == 2)  || (COLUMN == 1 && PROCESS != 2) */
