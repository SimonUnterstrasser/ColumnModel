import numpy as np
import math
import sys
#GCCif (KERNEL == 1 && LongK_Options == 2)
import Sedimentation as SD
#GCCendif /* (KERNEL == 1 && LongK_Options == 2) */

#GCCif (KERNEL == 1)

def LongKernel(rho_w):
    #rho_w      IN: density of water in kg/m^3
    #cck       OUT: kernel values in m^3/s (WELLMIXED == 0) bzw m^2 (WELLMIXED > 0)
    #mass_grid OUT: mass values in kg

    #### CHANGE: include absolute path of file location
    fp="data/Long_Kernel_fein/"

    #Einlesen des Radius in Einheit um
    file=np.loadtxt(fp+"Long_radius.txt")
    #print(file.shape)
    n=400
    radius_grid=np.zeros(400)
    radius_grid[0:400:3]=file[0:134,0]
    radius_grid[1:400:3]=file[0:133,1]
    radius_grid[2:400:3]=file[0:133,2]
    print('kernel values given for radius range: min ', radius_grid[0], ' max ', radius_grid[399])
    #Berechnung der Masse
    mass_grid=radius_grid**3*1.3333*math.pi*rho_w*1e-18

#GCCif (LongK_Options == 2)
    sys.exit('check this option carefully')
    #Parameter
    k1=4.5*1.e-4 #in 1/um**2
    k2=3.            #in um
    #Berechnung der Kollisionseffizienzen nach Long (1973)
    ec=np.zeros([n,n])
    for i in range(0,n):
        for j in range(0,i+1):
            if(radius_grid[j]<=50):
                ec[i,j]=k1*radius_grid[j]*radius_grid[j]*(1-k2/(max([radius_grid[i],3.0])+0.01))
            else:
                ec[i,j]=1.
            ec[j,i]=ec[i,j]

    #Berechnen der Fallgeschwindigkeiten
    rr=radius_grid*1e-6 # Input Radius in m
    winf=SD.Fallg(rr)*1e2  # Umrechnen der Einheiten von m/s in cm/s

    rr=radius_grid*1e-4 # im weiteren Radius in cm
    cck=np.zeros(n,n)
    #Berechnen der Look-up Tabelle des Kernels
    for i in range(0,n):
        for j in range(i+1,n):
            cck[i,j]=math.pi*(rr[j]+rr[i])*(rr[j]+rr[i])*ec[i,j]*abs(winf[j]-winf[i])
            cck[j,i]=cck[i,j]
    cck=cck*1e-6
#GCCendif

#GCCif (LongK_Options == 1)
        #GCCif (WELLMIXED == 0)
    #Einlesen der Kernel-Werte
    file_cck=np.loadtxt(fp+"Values_Longkernel.txt") # given in cm^3/s
    nxn=file_cck.size-2
    file_cck=file_cck.flatten()[0:nxn]
    cck=np.reshape(file_cck,[n,n])
    cck=cck*1e-6
        #GCCendif /* (WELLMIXED == 0)*/
        #GCCif (WELLMIXED > 0)
    #Einlesen der Efficiency-Werte
    file_eck=np.loadtxt(fp+"Efficiency_Longkernel.txt") # given in 1
    nxn=file_eck.size-2
    file_eck=file_eck.flatten()[0:nxn]
    eck=np.reshape(file_eck,[n,n])
    #Berechnen der 2D-Kernel-Werte in Einheit m^2
    rr=radius_grid*1e-6
    cck=np.zeros([n,n])
    for i in range(0,n):
        for j in range(i+1,n):
            cck[i,j]=math.pi*(rr[j]+rr[i])*(rr[j]+rr[i])*eck[i,j]
            cck[j,i]=cck[i,j]
        #GCCendif /* (WELLMIXED > 0)*/
#GCCendif

    return cck, mass_grid

#GCCendif /* (KERNEL == 1) */


#GCCif (KERNEL == 2)

def HallKernel(rho_w):
    #rho_w      IN: density of water in kg/m^3
    #cck       OUT: kernel values in m^3/s
    #mass_grid OUT: mass values in kg


#GCCif (KERNEL_INTPOL <= 1) /* logarithmic mass bin*/
    #### CHANGE: include absolute path of file location
    fp="data/Hall_Kernel_fein/"
    #Einlesen des Radius in Einheit um
    file=np.loadtxt(fp+"Hall_radius.txt")
    n=400
    radius_grid=np.zeros(n)
    radius_grid[0:n:3]=file[0:134,0]
    radius_grid[1:n:3]=file[0:133,1]
    radius_grid[2:n:3]=file[0:133,2]
    print('kernel values given for radius range: min ', radius_grid[0], ' max ', radius_grid[399])
    clip = 2
#GCCelse  /* linear mass bin*/
    #### CHANGE: include absolute path of file location
    fp="data/Hall_Kernel_Alfonsogrid/"
    #Einlesen des Radius in Einheit um
    file=np.loadtxt(fp+"Hall_radius.txt")
    n=42
    radius_grid=np.zeros(n)
    radius_grid[0:n:3]=file[0:14,0]
    radius_grid[1:n:3]=file[0:14,1]
    radius_grid[2:n:3]=file[0:14,2]
    print('kernel values given for radius range: min ', radius_grid[0], ' max ', radius_grid[39])
    clip = 0
#GCCendif /* (KERNEL_INTPOL > 1) */

    #Berechnung der Masse
    mass_grid=radius_grid**3*1.3333*math.pi*rho_w*1e-18

    #Einlesen der Werte
    file_cck=np.loadtxt(fp+"Values_Hallkernel.txt") # given in cm^3/s
    nxn=file_cck.size - clip
    file_cck=file_cck.flatten()[0:nxn]
    cck=np.reshape(file_cck,[n,n])

    return cck*1e-6, mass_grid

#GCCendif /* (KERNEL == 1) */
