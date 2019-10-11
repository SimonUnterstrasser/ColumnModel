import numpy as np
import Misc as FK

#GCCifdef COMP
#GCCdefine ICOMPPLOT 1
#the following statement includes a parameter file via a preprocessor directive
#GCCinclude "params.txt"
#GCCendif

#GCCifdef PLOT
#the following statement includes a parameter file via a preprocessor directive
#GCCinclude "params_plot.txt"
#GCCdefine ICOMPPLOT 2
#GCCendif


#GCCif (IREF == 1)
    #GCCif (KERNEL == 0)

# analytische Referenzloesung fuer Golovin-Kernel
def PlotGV_0D(ZPit_vec, ifirstGV=None, ilastGVonly=None):
    import matplotlib.pyplot as plt

    #Einlesen der Referenzloesung
    print('Einlesen der Golovin Loesung')
    # Momente festlegen, Berechnet aus analytischer Loesung
    #       ;analytische Loesung: Massenverteilung einlesen alle 200s
    #         Analytical_Numerical_Solutions/Golovin_N240/Golovin_analyticsolution_3600s_i200s.txt ;  10um Mean radius
    #        Analytical_Numerical_Solutions/Golovin/Golovin_analyticsolution_3600s_i200s.txt       ; 9.3um Mean radius

    #### CHANGE: include absolute path of file location
    fp_data='data/Golovin/'

    file=open(fp_data + 'Golovin_analyticsolution_3600s_i200s_p.txt','r')
    p=np.loadtxt(file)
    nt=int(p[0])
    nbin=int(p[1])
    print('nt,nbin',nt,nbin)
    file.close()

    file=open(fp_data + 'Golovin_analyticsolution_3600s_i200s_t.txt','r')
    time_vec=np.loadtxt(file) # contains times in seconds
    file.close()
    print('time_vec',time_vec)

    file=open(fp_data + 'Golovin_analyticsolution_3600s_i200s_r.txt','r')
    r=np.loadtxt(file)
    radius_vec=r.ravel()
    #print('radius_vec',radius_vec)

    file=open(fp_data + 'Golovin_analyticsolution_3600s_i200s_g.txt','r')
    listoflines = file.readlines()
    items = ' '.join(listoflines).split()
    file.close()
    print(len(items))
    g_vec=np.zeros([nt,nbin])
    for i in range(0,nt):
        g_vec[i,:]=np.array(items[nbin*i:nbin*(i+1)])

    #Plotten der Referenzloesung
    if(ifirstGV == 1):
        plt.plot(radius_vec[0:nbin],g_vec[0,0:nbin]*1e3,"k--")
    iPlot=1
    ntPlotREF=nt
    ntPlot=len(ZPit_vec)
    print(ntPlotREF,ntPlot)
    print(time_vec)
    print(ZPit_vec)
    for i in range(0,nt):
        #Zum Plotten wird der Radius in um umgerechnet
        print(i,iPlot,time_vec[i], ZPit_vec[iPlot])
        if abs(time_vec[i] - ZPit_vec[iPlot]) < 0.001:
            print('plotte Zeitschritt',time_vec[i], ZPit_vec[iPlot],i,iPlot)
            plt.plot(radius_vec[0:nbin],g_vec[i,0:nbin]*1e3,"--", color='k', label='Golovin')
            iPlot += 1
            iPlot =min([iPlot,ntPlot-1])
            print('yes', iPlot)

def defineMoments_0D():
    #solution for 9.3 um init
    print("read Golovin moments")
    nt_ref=7
    GolovinMomente=np.zeros([nt_ref,4])
    timeGolovin=np.arange(nt_ref)*600 # in sekunden  
    GolovinMomente[:,0]=np.array([296.8e6, 120.7e6, 490.8e5, 199.6e5, 811.4e4, 329.9e4, 134.1e4])  # ; Anzahl in m^-3
    
    GolovinMomente[:,1]=np.array([1.00008e-3 ,1.00183e-3 , 1.00293e-3,1.00184e-3, 1.00101e-3,1.00063e-3, 1.00046e-3])/1.00008e-3 
    GolovinMomente[:,1]=np.zeros(7)+1.0
    #; Masse
    
    GolovinMomente[:,2]=np.array([6.739e-15 ,4.094e-14 , 2.472e-13 , 1.493e-12 , 9.029e-12 ,5.462e-11 , 3.304e-10])# ; Moment 2 in kg^2/m^3
    GolovinMomente[:,3]=np.array([6.812e-26 , 3.993e-24 ,1.676e-22, 6.457e-21, 2.411e-19, 8.898e-18, 3.267e-16] )#; ;Moment 3 in kg^3/m^3
    return timeGolovin, GolovinMomente
    #Ende Golovin-Kernel

    #GCCendif /* (KERNEL == 0) */


    #GCCif ( KERNEL == 1 || KERNEL == 2 )

def read_GVdata_Wang():

# kn_name entweder 'Hall' oder 'Long'

    #GCCif (KERNEL == 1 )
    kn_name = 'Long'
    #GCCendif /* (KERNEL == 1 ) */
    #GCCif (KERNEL == 2 )
    kn_name = 'Hall'
    #GCCendif /* (KERNEL == 2 ) */
    # 572 ist die maximale Binanzahl der eingelesenen Daten
    g_wang=np.zeros([7,572])
    r_wang=np.zeros([7,572])
    nr_wang=np.zeros(7,dtype='int')

    # Spalte 1: r ist gegeben in mm
    # Spalte 5: g(ln(r)) ist gegeben in g/m^3

    fn_part_vec = ("00","10","20","30","40","50","60")
    nr_time=7
    #### CHANGE: include absolute path of file location
    fp = "data/Wang_Reference_OutputData/"
    for i_time,fn_part in enumerate(fn_part_vec):
        fn="GQ_" +fn_part + "_"+ kn_name +"_s16.dat"
       # print('read file: ', fn)
        read_tmp=np.genfromtxt(fp+fn)
        nr_bins, nr_col = read_tmp.shape
        print('Wang ', kn_name, ', nr_bins, nr_col: ', nr_bins, nr_col, fn)
        g_wang[i_time,:nr_bins]=read_tmp[:,4]
        r_wang[i_time,:nr_bins]=read_tmp[:,0]
        nr_wang[i_time] = nr_bins

    print('nr_wang: ', nr_wang)
    return r_wang,g_wang,nr_wang


#def PlotGV_0D(ZPit_vec,ifirstGV=None,ilastGVonly=None):

    #import matplotlib.pyplot as plt
    ##Einlesen der Referenzloesung
    #r_wang,g_wang,nr_wang = read_GVdata_Wang()
    ##Plotten der Referenzloesung
    #nr_GVplot= len(ZPit_vec)
    #iGV_select = range(1,nr_GVplot)
    #if (ifirstGV is not None): iGV_select = range(0,nr_GVplot)
    #iGV_select = range(0,nr_GVplot)
    #if (ilastGVonly == 1): iGV_select=[nr_GVplot-1]
    #for i in iGV_select:
        ##Zum Plotten wird der Radius in um umgerechnet
        #if abs(int(ZPit_vec[i]/600) - ZPit_vec[i]/600) < 0.001:
            #iii=int(ZPit_vec[i]/600)
            #print('plotte Zeitschritt',ZPit_vec[i],iii)
            #plt.plot(r_wang[iii,0:nr_wang[iii]]*1e3,g_wang[iii,0:nr_wang[iii]], color='k', label='Wang')

def defineMoments_0D():

    #verwende tabellierte Werte aus Wang-Paper
    nt_ref=7
    WangMomente=np.zeros([nt_ref,4])
    timeWang=np.arange(nt_ref)*600 # in Sekunden

    #GCCif (KERNEL == 1 )
    WangMomente[:,0]=np.array([295.4e6     , 287.4e6     , 278.4e6     ,  264.4e6       , 151.7e6, 13.41e6, 1.212e6])  # ; Anzahl in m^-3
    WangMomente[:,1]=np.array([0.999989e-3 , 0.999989e-3 , 0.999989e-3 ,  0.999989e-3   ,0.999989e-3, 0.999989e-3, 0.999989e-3])/0.999989e-3 #; Masse
    WangMomente[:,2]=np.array([6.739e-15   , 7.402e-15   , 8.72e-15,      3.132e-13      , 3.498e-10, 1.068e-8,3.199e-8])# ; Moment 2 in kg^2/m^3
    WangMomente[:,3]=np.array([6.813e-26   , 9.305e-26   ,5.71e-25,       3.967e-20       , 1.048e-15, 2.542e-13, 1.731e-12] )#; ;Moment 3 in kg^3/m^3
    #GCCendif /* (KERNEL == 1 ) */
    #GCCif (KERNEL == 2 )
    WangMomente[:,0]=np.array([295.4e6 , 287.8e6, 279.9e6, 270.2e6, 231.7e6, 124.5e6, 73.66e6])  # ; Anzahl in m^-3
    WangMomente[:,1]=np.array([0.999989e-3 ,0.999989e-3 ,0.999989e-3,0.999989e-3,0.999989e-3, 0.999989e-3, 0.999989e-3])/0.999989e-3  #; Masse
    WangMomente[:,2]=np.array([6.739e-15 , 7.184e-15, 7.999e-15, 7.827e-14, 1.942e-11, 7.928e-10,6.997e-9])# ; Moment 2 in kg^2/m^3
    WangMomente[:,3]=np.array([6.813e-26 , 8.282e-26,3.801e-25,2.531e-21, 6.107e-18, 2.108e-15, 1.221e-13] )#; ;Moment 3 in kg^3/m^3
    #GCCendif /* (KERNEL == 2 ) */
    return timeWang, WangMomente

    #GCCendif /* ( KERNEL == 1 || KERNEL == 2 ) */
#GCCendif /* (IREF == 1) */

#GCCif (IREF == 9 || BIN == 1)

##########################################################
#
# Routines for reading BinColumnModel SimulationData
#
##########################################################

import gzip
import os, sys

def get_RefMetaData(isimREF=0, fp=None):

    if fp is None:
        if type(fp_ref) is list:
            fp = fp_ref[isimREF]
            print('read from REF folder: ', fp)
        else:
            fp = fp_ref
    else:
        print('take provided fp: ', fp)

    dat = FK.openfile(fp + 'Meta.txt',options='rb')
    strline = dat.readline().split()
    [nbin, nt, nz, dz, dt, Tsim]=np.array(strline,dtype='float')
    strline = dat.readline().split()
    [LWC, Ntot, rmean, mmean, scal, dlnr, ikernel, i_init, i_process, i_bc_periodic] = \
        np.array(strline,dtype='float')
    strline = dat.readline().split()
    #Zeitlich konstante inflow_SD am Oberrand in kg/m^3
    g_init_top = np.array(strline, dtype='float')
    strline = dat.readline().split()
    #Radius grid in m
    rgrid = np.array(strline, dtype='float')
    #Mass grid in kg
    strline = dat.readline().split()
    mgrid = np.array(strline, dtype='float')
    dat.close()

    return int(nt), int(nz), dz, dt, Tsim, int(nbin), int(scal), dlnr, rgrid, mgrid

def get_RefProfileData(nz, nt, isimREF=0, fp=None, iswap=0):
#>>>>>>>>>>>>>> read profile data of 1D Bott bin model >>>>>>>>>>>>>>>>>>>>>>

    Mom_bin = np.zeros([4, nt, nz])
    fn = 'Profiles.txt'
    if fp is None:
        if type(fp_ref) is list:
            fp = fp_ref[isimREF]
        else:
            fp = fp_ref
    else:
        print('take provided fp: ', fp)

    fn = FK.filename(fn, fp)
    for icol in range(4):
         Mom_bin[icol,:,:] = np.loadtxt(fp+fn, usecols=[icol]).reshape(nt,nz)

    if (iswap == 1):
        Mom_bin = np.expand_dims(np.transpose(Mom_bin, (1,2,0)), axis=0)
        # returns array of shape (1,nt,nz,4)

    return Mom_bin

def get_RefFluxData(nt, nbin):
#>>>>>>>>>>>>>> read flux data of 1D Bott bin model >>>>>>>>>>>>>>>>>>>>>>

    fn_vec = ('Fluxes_in.txt', 'Fluxes_in_acc.txt', 'Fluxes_out.txt', 'Fluxes_out_acc.txt')
    fluxesREF = np.zeros([4, 4, nt])  # 4 types of fluxes, 4 Moments, nt times
    for ii, fn in enumerate(fn_vec):
        fn_tmp = FK.filename(fn, fp_ref)
        flux_input = np.genfromtxt(fp_ref + fn_tmp, comments='----------')
        #print(flux_input.shape)
        for i in range(nt):
            for iMom in range(4):
                fluxesREF[ii, iMom, i] = flux_input[nbin*(i+1) + i, iMom]
                #if (ii == 0) and (iMom == 0): print(nt,i,nbin*(i+1)+i)
    return fluxesREF

def get_RefSDdata(nt, nz, nbin):
#>>>>>>>>>>>>>> read size distribution data of 1D Bott bin model >>>>>>>>>>>>>>>>>>>>>>

    dat = FK.openfile(fp_ref + 'SD.txt')

    SizeDistrREF=np.zeros([nt, nz+1, nbin])
    time=np.zeros(nt)

    for t in range(0,nt):
        string_time = dat.readline()
        print('time in SD: ', string_time)
        time[t] = float(string_time)
        for gb in range(0,nz+1):
            strline = dat.readline().split()
            SizeDistrREF[t,gb,:] = np.array(strline, dtype='float')
    dat.close()
    return SizeDistrREF
#GCCendif /* (IREF == 9) */