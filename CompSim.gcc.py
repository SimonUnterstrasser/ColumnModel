#!/usr/bin/env python3
# -*- coding: utf-8 -*-
########################################################################
#
# implementation of AON algorithm for collisional growth of hydrometeors.
# can be employed in a box model or column model
# various variants of the algorithm are covered here (LinSamp, explicit overtakes)
#
########################################################################

#Python modules
import os, sys
import numpy as np
import math
import warnings
import random
import time

#Aggregation modules
import SIPinit as SI
import AON_Alg as AON
import Misc as FK
#GCCif (KERNEL == 1 || KERNEL == 2)
import Kernel as K
#GCCendif /* (KERNEL == 1 || KERNEL == 2) */

#Sedimentation module
#GCCif (COLUMN == 1 && PROCESS != 2)
import Sedimentation as SD
#GCCendif /* (COLUMN == 1 && PROCESS != 2) */

#Plotting module
import PlotSim as PS

#GCCif (WELLMIXED > 0)
    #GCCif (PROCESS != 0)
print('if 2D-WELLMIXED approach is chosen, Sedimentation and Aggregation must be both activated (PROCESS = 0)')
sys.exit('STOPPED')
    #GCCendif  /* (PROCESS != 0) */
    #GCCif (KERNEL == !1 && KERNEL != 2)
print('if 2D-WELLMIXED approach is chosen, hydrodynmic kernel must be chosen')
sys.exit('STOPPED')
    #GCCendif  /* (KERNEL == !1 && KERNEL != 2) */
    #GCCif (LINEAR != 0)
print('if 2D-WELLMIXED approach is chosen, linear sampling option is not possible')
sys.exit('STOPPED')
    #GCCendif  /* (LINEAR != 0) */
    ####GCCif (AGG_MC != 0)
###print('if 2D-WELLMIXED approach is chosen, multiple collection option must be switched off')
###sys.exit('STOPPED')
    ####GCCendif  /* (AGG_MC != 0) */
    #GCCendif  /* (WELLMIXED > 0) */

warnings.filterwarnings("ignore",category =RuntimeWarning)

dV_skal = 1
nr_SIPs_skal = 1

#the following statement includes a parameter file via a preprocessor directive
#GCCinclude "params.txt"

#>>>>>>>>>> derived parameters >>>>>>>>>>>>>>>>>>>>

#GCCif (INITV == 1)
    #initial size distribution
        #Mean mass in kg
xf0=FK.r2m(r0,const_mass2rad)
        #initial concentration in m^-3
        #N0=239e6
N0=LWC/xf0
        #number of bins
n=n10*r10
#GCCendif /* (INITV == 1) */

#GCCif (COLUMN == 1)
    #GCCif (PROCESS == 0)
i_process = 0
    #GCCendif /* (PROCESS == 0) */
    #GCCif (PROCESS == 1)
i_process = 1
    #GCCendif /* (PROCESS == 1) */
    #GCCif (PROCESS == 2)
i_process = 2
    #GCCendif /* (PROCESS == 2) */
#GCCendif /* (COLUMN == 1) */

    #total number of time steps
iend=int(Tsim/dt) + 1
#<<<<<<<<<<< derived parameters <<<<<<<<<<<<<<<<<<<<<


#GCCif (KERNEL == 0 || KERNEL == 3)
cck=0
m_low=0
eta_indize=0
m_kernel = 0
ikernel=0
#GCCendif /* (KERNEL == 0 || KERNEL == 3) */

#GCCif (KERNEL == 1)
#Aufruf Long Kernel
ikernel=1
[cck,m_kernel]=K.LongKernel(rho_w)
eta_indize=m_kernel[1]/m_kernel[0]
#m_low=m_kernel[0]/eta_indize
m_low=m_kernel[0]

#GCCendif /* (KERNEL == 1) */

#GCCif (KERNEL == 2)
#Aufruf Hall-Kernel
ikernel=2
[cck,m_kernel]=K.HallKernel(rho_w)
eta_indize=m_kernel[1]/m_kernel[0]
m_low=m_kernel[0]
#GCCendif /* (KERNEL == 2) */

fLog_currColls= None
fLog_accColls = None
fLog_currColls_WM2D = None
fLog_Combs = None
fLog_p = None
#GCCif (COUNT_COLLS == 1)
fLog_currColls= open('log_currColls.dat','w')
fLog_accColls = open('log_accColls.dat','w')
fLog_Combs = open('log_Combs.dat','w')
fLog_p = open('log_p.dat','w')
    ##GCCif (WELLMIXED > 0)
fLog_currColls_WM2D= open('log_currColls_WM2D.dat','w')
    ##GCCendif /* (WELLMIXED > 0)*/
#GCCendif  /* (COUNT_COLLS == 1) */

#GCCif (COLUMN == 0)
#================================ BOX-MODEL ===========================================
#>>>>>>>>>>>> start of box model section>>>>>>>>>>>>>>>>>>>>
nz = 1

#GCCif (DISCRETE == 0)
# definition of bin grid for size distribution plots
mfix_plot=np.zeros(nplot)
mdelta_plot=np.zeros(nplot)
for i in range(0,nplot-1):
    mfix_plot[i+1]=10**(i/n10_plot+min10_plot)
    mdelta_plot[i]=mfix_plot[i+1]-mfix_plot[i]
#GCCendif /* (DISCRETE == 0) */

imod_GVplot = int(t_intervall_GVplot/dt)
nr_GVplot = int(1 + (Tsim-t_start_GVplot)/t_intervall_GVplot)
t_end_GVplot = t_start_GVplot + (nr_GVplot-1)*t_intervall_GVplot
t_vec_GVplot = np.arange(t_start_GVplot,t_end_GVplot+1,t_intervall_GVplot) #times at which SIP data is saved and size distributions can be plotted

nEK_sip_plot=np.zeros([nr_inst,nr_GVplot,nr_sip_max])
    #stores SIP data of all realisations at the times defined in t_vec_GVplot

#GCCif (DISCRETE == 0)
mEK_sip_plot=np.zeros([nr_inst,nr_GVplot,nr_sip_max])
#GCCendif /* (DISCRETE == 0) */
#GCCif (DISCRETE >= 1)
mEK_sip_plot=np.zeros([nr_inst,nr_GVplot,nr_sip_max],dtype='int')
#GCCendif /* (DISCRETE >= 1) */
nr_SIPs_plot=np.zeros([nr_inst,nr_GVplot],dtype='int')

imod_MOMsave=int(t_intervall_MOMsave/dt)
nr_MOMsave = int(1 + (Tsim-t_start_MOMsave)/t_intervall_MOMsave)
t_end_MOMsave = t_start_MOMsave + (nr_MOMsave-1)*t_intervall_MOMsave
t_vec_MOMsave=np.arange(t_start_MOMsave,t_end_MOMsave+1,t_intervall_MOMsave)
    #usually a finer time grid is used for saving the moment data

print(t_vec_MOMsave)
print(nr_MOMsave)

nr_MOMs=4  # evaluate moments of order 0 to nr_MOMs-1
MOMsave=np.zeros([nr_inst,nr_MOMsave,nr_MOMs])

#GCCif (DISCRETE == 0)
outp_format='%1.6e'
outp_format_long='%1.6e'
#GCCendif /* (DISCRETE == 0) */
#GCCif (DISCRETE >= 1)
outp_format='%4i'
outp_format_long='%6i'
#GCCendif /* (DISCRETE >= 1) */

fMom = open('Moments.dat', 'wb')
fGV  = open('SIP.dat', 'wb')
fLog= open('log.txt','a')

starttime = time.time()
localtime = time.asctime( time.localtime(starttime) )
print("Start time computation:", localtime)

fLog.write(os.getcwd()+ '\n')
fLog.write("Start time computation: "+ localtime+ '\n')
fLog.close()

# loop over all realisations
for k in range(0,nr_inst):
    if (iPM >= 1): print(os.getcwd())

    if (k%50 == 0):
        print('-------------------------------------- new instance ',k,'-------')
    fLog= open('log.txt','a')
    fLog.write('instance: ' + str(k)+ '\n')
    fLog.close()

    i_MOMsave=1
    i_GVplot=1
    random.seed(32+(k+nr_ins_start)*123433)
    count_colls=0
    #GCCif (COUNT_COLLS == 1)
    count_colls=np.zeros(20,dtype=np.uint64)
    #GCCendif  /* (COUNT_COLLS == 1) */
    ###########################################################################
    #
    # Initialise SIP ensemble
    #
    ###########################################################################

#GCCif (INITV == 1)
    [nr_SIPs,nEK_sip_ins,mEK_sip_ins,EK_krit,mdelta,mfix]=SI.InitSIP_ExpVert_singleSIP_WS(imlow,n10,r10,min10,eta_nu,xf0,N0,dV,nr_sip_max)
#GCCendif /* if (INITV == 1) */
#GCCif (INITV == 2)
    nr_SIPs, nEK_sip_ins, mEK_sip_ins = SI.InitSIP_Alfonso(dV_skal=dV_skal) #20 droplets with 17um and 10 droplets with 21.4um or similar setups
#GCCendif /* if (INITV == 2) */

    MOMsave[k,0,:] = FK.Moments_k0_3(nEK_sip_ins,mEK_sip_ins)
    if (iPM >= 1): print('initial moments: ', MOMsave[k,0,:])
    if (iPM >= 1): print('nr_SIPs: ', nr_SIPs)
    np.savetxt(fMom, MOMsave[k,0,:].reshape(1,4), fmt=outp_format_long) # reshape necessary to write all 4 values in one row
    nr_SIPs_plot[k,0 ]= nr_SIPs
    nEK_sip_plot[k,0,0:nr_SIPs]=nEK_sip_ins
    mEK_sip_plot[k,0,0:nr_SIPs]=mEK_sip_ins
    np.savetxt(fGV, [nr_SIPs], fmt='%5i')
    np.savetxt(fGV, nEK_sip_ins[0:nr_SIPs].reshape(1,nr_SIPs), fmt=outp_format)
    np.savetxt(fGV, mEK_sip_ins[0:nr_SIPs].reshape(1,nr_SIPs), fmt=outp_format)
    #np.savetxt(fGV, nEK_sip_ins[0:nr_SIPs].reshape(nr_SIPs,1), fmt=outp_format)
    #np.savetxt(fGV, mEK_sip_ins[0:nr_SIPs].reshape(nr_SIPs,1), fmt=outp_format)

    ibreak = 0
    iibreak = 0
    #>>>>>>>>>>>>> time iteration
    for it in range(0,iend):
        t=it*dt
        if (nr_SIPs > 1):
            #print('nr_SIPs',nr_SIPs,type(nr_SIPs))
            ibreak = AON.Aggregation(nr_SIPs,nEK_sip_ins,mEK_sip_ins,count_colls,cck,m_low,eta_indize,m_kernel)

        if (it%imod_MOMsave ==0) & (it != 0):
            MOMsave[k,i_MOMsave,:]=FK.Moments_k0_3(nEK_sip_ins,mEK_sip_ins)
            if (iPM >= 2): print(it,i_MOMsave,MOMsave[k,i_MOMsave,:])
            np.savetxt(fMom, MOMsave[k,i_MOMsave,:].reshape(1,4), fmt=outp_format_long)
            i_MOMsave = i_MOMsave+1

        #>>>>>>>>>>>>>>>>>>remove zero weight SIPs, if present at all>>>>>>>>>>>>>>>>>>>>>
        if (min(nEK_sip_ins) == 0):
            index_list=nEK_sip_ins.nonzero()
            #print('nr_SIPs old: ',nr_SIPs)
            nEK_sip_ins=nEK_sip_ins[index_list]
            mEK_sip_ins=mEK_sip_ins[index_list]
            nr_SIPs = nEK_sip_ins.size
            #print('nr_SIPs new: ',nr_SIPs)
            if (nr_SIPs == 1):
                print('only one SIP remains, stop computation of current instance, proceed with next instance' )
                print('nu: ', nEK_sip_ins, 'm: ', mEK_sip_ins, 'time:', t, 'it: ', it)
                ibreak = 1

        #>>>>>>>>save SIP data at specified points in time>>>>>>>>>>>>>
        if (it%imod_GVplot == 0) & (it != 0):
            nr_SIPs_plot[k,i_GVplot ]= nr_SIPs
            nEK_sip_plot[k,i_GVplot,0:nr_SIPs]=nEK_sip_ins
            mEK_sip_plot[k,i_GVplot,0:nr_SIPs]=mEK_sip_ins
            if (iPM >= 2):
                print('SIP output, #, it, time:', i_GVplot, it, it*dt)
            np.savetxt(fGV, [nr_SIPs], fmt='%5i')
            np.savetxt(fGV, nEK_sip_plot[k,i_GVplot,0:nr_SIPs].reshape(1,nr_SIPs), fmt=outp_format)
            np.savetxt(fGV, mEK_sip_plot[k,i_GVplot,0:nr_SIPs].reshape(1,nr_SIPs), fmt=outp_format)
            #np.savetxt(fGV, nEK_sip_ins[0:nr_SIPs].reshape(nr_SIPs,1), fmt=outp_format)
            #np.savetxt(fGV, mEK_sip_ins[0:nr_SIPs].reshape(nr_SIPs,1), fmt=outp_format)
            i_GVplot = i_GVplot+1

        #if (ibreak == 1):
            #print('break condition met at iteration it = ', it)

    #>>>>>>end of time iteration of a single instance>>>>>>>>>

    #>>>>>>>>>>analyse computing time>>>>>>>>>>>>>>>>>>>>
    currenttime = time.time()
    currenttime_str = time.asctime( time.localtime(currenttime))
    endtime_expected=starttime+ (nr_inst/(k+1)) * (currenttime-starttime)
    endtime_expected_str = time.asctime( time.localtime(endtime_expected))
    if (iPM >= 2):
        print("Instance {} of {} finished".format(k+1,nr_inst))
        print("Start time/Current time/Expected end time: ")
    #print(localtime, ' --- ', currenttime_str, ' --- ', endtime_expected_str)
    fLog= open('log.txt','a')
    fLog.write('total computing time in sec: '+ str(int(currenttime-starttime)) + '\n')
    fLog.write(currenttime_str+ ' --- '+ endtime_expected_str+ '\n')

#GCCif (COUNT_COLLS == 1)
    cc_sum=count_colls.sum()
    fLog.write('a '+"{}".format(cc_sum)+'\n')
    fLog.write('b '+" ".join("{}".format(x) for x in count_colls)+'\n')
    cc_frac=count_colls/cc_sum
    fLog.write('b '+" ".join("{:.4}".format(x) for x in cc_frac)+'\n')
#GCCendif  /* (COUNT_COLLS == 1) */
    fLog.close()
#<<<<<<<<<<<<<<<<<end of loop of all realisations<<<<<<<<<<<<<<<<<

fMom.close()
fGV.close()
localtime = time.asctime( time.localtime(time.time()) )
print("End time computation:", localtime)

#>>>>>>>>>>>>>>>>Generate meta data>>>>>>>>>>>>>
fMom = open('Moments_meta.dat', 'wb')
np.savetxt(fMom,np.array([dV,skal_m]),fmt = '%10e')
np.savetxt(fMom,np.array([nr_inst,nr_MOMsave]),fmt = '%4d')
np.savetxt(fMom, t_vec_MOMsave.reshape(1,nr_MOMsave),fmt = '%6d')
fMom.close()
fGV  = open('SIP_meta.dat', 'wb')
np.savetxt(fGV,np.array([nr_inst,nr_GVplot]),fmt = '%4d')
np.savetxt(fGV, t_vec_GVplot.reshape(1,nr_GVplot),fmt = '%6d')
fGV.close()
#>>>>>>>Generate data file with mean moments>>>>>>>>>>>>>>>
FK.CIO_MOMmean(data_in=MOMsave,fp_out='',skal_m=skal_m,dV=dV)

fLog= open('log.txt','a')
currenttime = time.time()
fLog.write('total computing time in sec: '+ str(int(currenttime-starttime)) + '\n')
fLog.write('finalised\n')
fLog.close()

#----------------------------------------------------------------------------------
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>generate plot for a first analysis>>>>>>>>>>>>>>>>>>>>>>>


for i in range (4): MOMsave[:,:,i] = MOMsave[:,:,i]*skal_m**i/dV
PS.PlotMoments(MOMsave,t_vec_MOMsave)
PS.PlotGV(nEK_sip_plot,mEK_sip_plot,nr_SIPs_plot,t_vec_GVplot,skal_m=skal_m)
#GCCif (IREF >= 2)
PS.PlotRelStdDev(mEK_sip_plot, t_vec_GVplot)
#GCCendif /* (IREF >= 2) */

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<end of box model section<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#=============================== BOX-MODEL =============================================
#GCCendif /* (COLUMN == 0) */



#GCCif (COLUMN == 1)
# start of column model section
#=============================== COLUMN-MODEL =============================================


imod_GVplot=int(t_intervall_GVplot/dt)
nr_GVplot = int(1 + (Tsim-t_start_GVplot)/t_intervall_GVplot)
t_end_GVplot = t_start_GVplot + (nr_GVplot-1)*t_intervall_GVplot
t_vec_GVplot=np.arange(t_start_GVplot,t_end_GVplot+1,t_intervall_GVplot) #times at which SIP data is saved and size distributions can be plotted

imod_MOMsave=int(t_intervall_MOMsave/dt)
nr_MOMsave = int(1 + (Tsim-t_start_MOMsave)/t_intervall_MOMsave)
t_end_MOMsave = t_start_MOMsave + (nr_MOMsave-1)*t_intervall_MOMsave
t_vec_MOMsave=np.arange(t_start_MOMsave,t_end_MOMsave+1,t_intervall_MOMsave)
    #usually a finer time grid is used for saving the moment data

print(t_vec_MOMsave)
print(nr_MOMsave)

nr_MOMs=4  # evaluate moments of order 0 to nr_MOMs-1
MOMsave=np.zeros([nr_inst,nr_MOMsave,nz,nr_MOMs])

# track fluxes of moments 0 to nr_MOMs-1 together with SIP number at the lower and upper boundary
# "instantaneous" flux: actually the flux is an average over a time period of length t_intervall_GVplot
# in the "SIP world" with discrete crossings of SIPs across the boundaries it is not reasonable/possible to determine instantaneous fluxes
# at time t=0 all fluxes are set to zero and not output, hence array of length nr_GVplot-1 is used
FluxIn    =np.zeros([nr_inst,nr_GVplot-1,nr_MOMs+1])
FluxOut   =np.zeros([nr_inst,nr_GVplot-1,nr_MOMs+1])
FluxInAcc =np.zeros([nr_inst,nr_GVplot-1,nr_MOMs+1])
FluxOutAcc=np.zeros([nr_inst,nr_GVplot-1,nr_MOMs+1])

#GCCif (TRACKCENTER == 1)
# last index has 3 components: before aggregation, after aggregation, after sedimentation
zCenterIC_full=np.zeros([nr_inst,nr_MOMsave,nz,nr_MOMs,3])
zCenterSIP_full=np.zeros([nr_inst,nr_MOMsave,nz,nr_MOMs,3])
zCenterIC_av=np.zeros([nr_MOMsave,nr_MOMs,3])
zCenterIC_av_denom=np.zeros([nr_MOMsave,nr_MOMs,3])
zCenterSIP_av=np.zeros([nr_MOMsave,nr_MOMs,3])
zCenterSIP_av_denom=np.zeros([nr_MOMsave,nr_MOMs,3])
nr_SIPs_GB_save=np.zeros([nr_inst,nr_MOMsave,nz])
#GCCendif  /* (TRACKCENTER == 1) */


outp_format='%1.6e'
outp_format_long='%1.6e'
outp_format_flux=('%1.4e','%1.4e','%1.4e','%1.4e','%9d')


#generate the general meta data file
fMeta = open('Meta.dat','w')
fMeta.write("{:5d} {:7.3f} {:7.1f} {:5d}\n".format(nz,dz,Tsim,nr_inst))
fMeta.write("{:7.3e} {:7.3e} {:7.3e} {:3d} {:3d} {:3d} \n".format(LWC,r0,xf0,ikernel,i_init_1D,i_process))
fMeta.close()

#generate the meta data files for particular meta data files
fMomMeta = open('Moments_meta.dat', 'wb')
np.savetxt(fMomMeta,np.array([dV,skal_m]),fmt = '%10e')
np.savetxt(fMomMeta,np.array([nr_inst,nr_MOMsave]),fmt = '%4d')
np.savetxt(fMomMeta, t_vec_MOMsave.reshape(1,nr_MOMsave),fmt = '%6d')
fMomMeta.close()
fGVMeta  = open('SIP_meta.dat', 'wb')
np.savetxt(fGVMeta,np.array([nr_inst,nr_GVplot]),fmt = '%4d')
np.savetxt(fGVMeta, t_vec_GVplot.reshape(1,nr_GVplot),fmt = '%6d')
fGVMeta.close()

#open output files
fMom = open('Moments.dat', 'wb')
fGV  = open('SIP.dat', 'wb')
fLog= open('log.txt','a')

fFluxIn     = open('Fluxes_in.dat'     , 'wb')
fFluxInAcc  = open('Fluxes_in_acc.dat' , 'wb')
fFluxOut    = open('Fluxes_out.dat'    , 'wb')
fFluxOutAcc = open('Fluxes_out_acc.dat', 'wb')

#GCCif (TRACKCENTER == 1)
fCenterIC_full = open('centerIC_full.dat' , 'wb')
fCenterSIP_full= open('centerSIP_full.dat', 'wb')
fCenterIC_av   = open('centerIC_av.dat' , 'wb')
fCenterSIP_av  = open('centerSIP_av.dat', 'wb')
#GCCendif  /* (TRACKCENTER == 1) */

#GCCif (TRACKOUT == 1)
fGVout  = open('SIPout.dat', 'wb')
fGVoutMeta  = open('SIPout_meta.dat', 'wb')
fGVoutMeta.close()
#GCCendif  /* (TRACKOUT == 1) */

#GCCif (INFLUX_TOP == 1)
fGVin  = open('SIPin.dat', 'wb')
fGVinMeta  = open('SIPin_meta.dat', 'wb')
fGVinMeta.close()
#GCCendif /* (INFLUX_TOP == 1)*/

starttime = time.time()
localtime = time.asctime( time.localtime(starttime) )
print("Start time computation:", localtime)

fLog.write(os.getcwd()+ '\n')
fLog.write("Start time computation: "+ localtime+ '\n')
fLog.close()

#Instanzschleife
for k in range(0,nr_inst):
    if (iPM >= 1): print(os.getcwd())

    print('-------------------------------------- new instance ',k,'-------')
    fLog= open('log.txt','a')
    fLog.write('instance: ' + str(k)+ '\n')
    fLog.close()

    #GCCif (TRACKOUT == 1)
    ntime_SIPout = 0
    #GCCendif  /* (TRACKOUT == 1) */
    #GCCif (INFLUX_TOP == 1)
    ntime_SIPin = 0
    #GCCendif /* (INFLUX_TOP == 1)*/

    i_MOMsave=1
    i_GVplot=1

    #GCCif (SPEED_VECTOR == 0)
    random.seed(32+(k+nr_ins_start)*123433)
    #GCCendif /* (SPEED_VECTOR == 0) */
    #GCCif (SPEED_VECTOR == 1)
    np.random.seed(32+(k+nr_ins_start)*123433)
    #GCCendif /* (SPEED_VECTOR == 1) */

    count_colls=0
    #GCCif (COUNT_COLLS == 1)
    count_colls=np.zeros(20,dtype=np.uint64)
    #GCCendif  /* (COUNT_COLLS == 1) */

    ###########################################################################
    #
    # Initialise SIP ensemble
    #
    ###########################################################################

    nEK_sip_ins = np.zeros(nr_sip_max)
    mEK_sip_ins = np.zeros(nr_sip_max)
    zEK_sip_ins = np.zeros(nr_sip_max)+zNan
    #GCCif (PROCESS != 2)
    wEK_sip_ins = np.zeros(nr_sip_max)
    #GCCendif /* (PROCESS != 2)*/
    zGBsep = np.array(np.arange(nz+1))*dz
    #print('zGBsep',zGBsep)
    iSIP_GBsep = np.zeros(nz+1,dtype='int')
    nr_SIPs_GB = np.zeros(nz,dtype='int')
    count_save = 0
    # ColumnModel-Initialization
        #  ! i_init:
        #  !    1 = empty domain
        #  !    2 = top GB only
        #  !    3 = top half domain
        #  !    4 = linearly decaying over total column from 2*g_init to 0
        #  !    5 = linearly decaying over top half from 2*g_init to 0
        #  !    6 = total domain
        #  !    7 = linearly decaying over top quarter from 2*g_init to 0
        #  !    8 = sin()-hill over top half with 2*g_init max peak
        #  !    9 = sin()-hill over top quarter with 2*g_init max peak

    if (i_init_1D == 1): i_start = nz + 1
    if (i_init_1D == 2): i_start = nz
    if (i_init_1D == 3): i_start = nz / 2
    if (i_init_1D == 4): i_start = 1
    if (i_init_1D == 5): i_start = nz / 2
    if (i_init_1D == 6): i_start = 1
    if (i_init_1D == 7): i_start = 3 * nz / 4
    if (i_init_1D == 8): i_start = nz / 2
    if (i_init_1D == 9): i_start = 3 * nz / 4
    i_start = i_start -1 # above block analagous to block in Botts F77 programm. Be aware that python indices start at 0, not at 1 like in Fortran

    for iz in range(0,nz):
        #print('init row ', iz,i_init_1D,i_start)
        if (i_init_1D <= 3) or (i_init_1D == 6): N0_pick = N0
        if (i_init_1D == 4) or (i_init_1D == 5) or (i_init_1D == 7) or (i_init_1D == 8) or (i_init_1D == 9):
            izp=iz+1
            if (i_init_1D == 4): gew = min(1.0,(1.0*izp/nz))
            if (i_init_1D == 5): gew = max(0.0,min(1.0,(2.0*(1.0*izp/nz-0.5))))
            if (i_init_1D == 7): gew = max(0.0,min(1.0,(4.0*(1.0*izp/nz-0.75))))
            if (i_init_1D == 8): gew = math.sin(max(0.0,min(1.0,(2.0*(1.0*(izp - 0.5 )/nz-0.50))))*3.14)
            if (i_init_1D == 9): gew = math.sin(max(0.0,min(1.0,(4.0*(1.0*(izp - 0.5 )/nz-0.75))))*3.14)
            #print('init ', iz, gew)
            N0_pick = N0 * 2 * gew
        if (iz >= i_start):
            [nr_SIPs,nEK_sip_tmp,mEK_sip_tmp,EK_krit,mdelta,mfix]=SI.InitSIP_ExpVert_singleSIP_WS(imlow,n10,r10,min10,eta_nu,xf0,N0_pick,dV,nr_sip_max,dV_skal=dV_skal)
            print(iz,nz,'iz,nz')
            # print('gew: ', gew)
            #[nr_SIPs,nEK_sip_tmp,mEK_sip_tmp]=SI.InitSIP_uniform(n10,r10,min10,const_mass2rad, gew=gew)
            nr_SIPs_GB[iz] = nr_SIPs
            if (count_save+nr_SIPs > nr_sip_max):
                print(count_save+nr_SIPs,nr_sip_max,iz)
                sys.exit("aktuelle SIP-Anzahl groesser als nr_sip_max")
            nEK_sip_ins[count_save:count_save+nr_SIPs]=nEK_sip_tmp[0:nr_SIPs]
            mEK_sip_ins[count_save:count_save+nr_SIPs]=mEK_sip_tmp[0:nr_SIPs]
            #print('count_save,nr_SIPs',count_save,nr_SIPs)
            for j in range(0,nr_SIPs):
                zEK_sip_ins[count_save+j]=(iz+random.random())*dz
                #zEK_sip_ins[count_save+j]=(i+0.5)*dz
        else:
            nr_SIPs = 0

        iSIP_GBsep[iz]=count_save
        count_save += nr_SIPs

    # total SIP number after end of column initialization
    nr_SIPs_tot=count_save
    #print('final nr_SIPs',count_save,nr_SIPs)
    iSIP_GBsep[nz]=count_save
    #iz=nz
    #print(iz,iSIP_GBsep[iz],zGBsep[iz])
    nr_SIPs_GB_max=nr_SIPs_GB.max()

    rSIP_tmp=FK.m2r(mEK_sip_ins[:nr_SIPs_tot],const_mass2rad)
    #GCCif (PROCESS != 2)
    wEK_sip_ins[:nr_SIPs_tot]=SD.Fallg(rSIP_tmp)
    #GCCendif /*(PROCESS != 2)*/

    #print('r',rSIP_tmp)
    #print('w',wEK_sip_ins[:nr_SIPs_tot])


#output at initialization

    for iz in range(0,nz):
        nr_SIPs = nr_SIPs_GB[iz]
        ia=iSIP_GBsep[iz]
        ie=iSIP_GBsep[iz+1]
        #Output
        if (nr_SIPs > 0): MOMsave[k,0,iz,:] = FK.Moments_k0_3(nEK_sip_ins[ia:ie],mEK_sip_ins[ia:ie])
        if (iPM >= 1): print('initial moments: ', MOMsave[k,0,iz,:])
        if (iPM >= 1): print('nr_SIPs: ', nr_SIPs)
        np.savetxt(fMom, MOMsave[k,0,iz,:].reshape(1,4), fmt=outp_format_long) # reshape necessary to write all 4 values in one row
        np.savetxt(fGV, np.array([nr_SIPs,iz]).reshape(1,2), fmt='%6i')
        if (nr_SIPs > 0):
            np.savetxt(fGV, nEK_sip_ins[ia:ie].reshape(1,nr_SIPs), fmt=outp_format)
            np.savetxt(fGV, mEK_sip_ins[ia:ie].reshape(1,nr_SIPs), fmt=outp_format)
            np.savetxt(fGV, zEK_sip_ins[ia:ie].reshape(1,nr_SIPs), fmt=outp_format)

    fGVMeta = open('SIP_meta.dat', 'a')
    fGVMeta.write("\n{:7d} {:5d}\n".format(nr_SIPs_GB.sum(),nr_SIPs_GB.max()))
    fGVMeta.close()


#prepare SIP processing

    #sort SIPs by z-position and find index range of each grid box
    #GCCif (INFLUX_TOP != 1)
    perm_sort=np.argsort(zEK_sip_ins[:nr_SIPs_tot+1]) #keep only SIPs in column + exactly one additional SIP with zNan (this extra SIPs is necessary for keeping get_isep routine simple)
    #GCCelse
    # do not shorten SIP arrays, if non-zero influx is chosen
    perm_sort=np.argsort(zEK_sip_ins)  # this statement keeps all nr_sip_max SIPs even if the SIP list is shorter than that
    #GCCendif /* (INFLUX_TOP != 1)*/

    nEK_sip_ins=nEK_sip_ins[perm_sort]
    mEK_sip_ins=mEK_sip_ins[perm_sort]
    zEK_sip_ins=zEK_sip_ins[perm_sort]
    #GCCif (PROCESS != 2)
    wEK_sip_ins=wEK_sip_ins[perm_sort]
    #GCCendif /*(PROCESS != 2)*/

    [iSIP_GBsep,nr_SIPs_GB] = FK.get_isep(zEK_sip_ins,zGBsep,nz)

    ausw_SIPs_active = np.flatnonzero(zEK_sip_ins <= zCol)
    nr_SIPs_tot=len(ausw_SIPs_active)
    print('nr_SIPs_tot',nr_SIPs_tot,len(perm_sort),zEK_sip_ins.min(),zEK_sip_ins.max(),zCol)

    #GCCif (TRACKCENTER == 1)
# track centers of ICs and SIPs
    tmp,nom,denom=FK.trackcenters(nEK_sip_ins,mEK_sip_ins,zEK_sip_ins,iSIP_GBsep,nr_SIPs_GB,nz,dz)
        #returns tmp: array of size (2,nz,nr_MOMs)
        #        nom: array of size (2,nr_MOMs)
        #        denom: array of size (2,nr_MOMs)
    zCenterIC_full [k,0,:,:,0]=tmp[0,:,:]
    zCenterSIP_full[k,0,:,:,0]=tmp[1,:,:]
    zCenterIC_full [k,0,:,:,1]=tmp[0,:,:]
    zCenterSIP_full[k,0,:,:,1]=tmp[1,:,:]
    zCenterIC_full [k,0,:,:,2]=tmp[0,:,:]
    zCenterSIP_full[k,0,:,:,2]=tmp[1,:,:]
    np.savetxt(fCenterIC_full,  zCenterIC_full[k,0,:,:,:].reshape(nz,4*3), fmt='%1.4f')
    np.savetxt(fCenterSIP_full, zCenterSIP_full[k,0,:,:,:].reshape(nz,4*3), fmt='%1.4f')
    nr_SIPs_GB_save[k,0,:]=nr_SIPs_GB

    zCenterIC_av[0,:,0]+=nom[0,:]
    zCenterIC_av_denom[0,:,0]+=denom[0,:]
    zCenterSIP_av[0,:,0]+=nom[1,:]
    zCenterSIP_av_denom[0,:,0]+=denom[1,:]
    zCenterIC_av[0,:,1]+=nom[0,:]
    zCenterIC_av_denom[0,:,1]+=denom[0,:]
    zCenterSIP_av[0,:,1]+=nom[1,:]
    zCenterSIP_av_denom[0,:,1]+=denom[1,:]
    zCenterIC_av[0,:,2]+=nom[0,:]
    zCenterIC_av_denom[0,:,2]+=denom[0,:]
    zCenterSIP_av[0,:,2]+=nom[1,:]
    zCenterSIP_av_denom[0,:,2]+=denom[1,:]
    #GCCendif  /* (TRACKCENTER == 1) */

    ibreak = 0
    iibreak = 0
    FluxOutAcctmp=np.zeros(5)
    FluxInAcctmp=np.zeros(5)
    ausw_zeroweights=np.flatnonzero(nEK_sip_ins <= 1e-50)
    ntmp_before=len(ausw_zeroweights)
    # Zeitschleife
    #iend=1
    for it in range(0,iend):
        t = it * dt
        #GCCif (INFLUX_TOP == 1)
        #let new SIPs fall into the domain, probabilistic implementation of an influx BC
        #generate a full SIP ensemble of a prescribed inflow SD
        [nr_SIPs_BC,nEK_sip_BC,mEK_sip_BC,EK_krit,mdelta,mfix]=SI.InitSIP_ExpVert_singleSIP_WS(imlow,n10,r10,min10,eta_nu,xf0,N0_pick,dV,nr_sip_max,dV_skal=dV_skal,silent=1)
        #compute how far each SIPs travels into the domain
        zdiff = SD.Fallg(FK.m2r(mEK_sip_BC,const_mass2rad))*dt
        # probabilities of SIPs being generated, i.e. fractional number of a SIP
        nrSIP_per_bin_frac = zdiff/dz*nr_SIPs_skal
        #print('nrSIP_per_bin_frac \n',nrSIP_per_bin_frac[mEK_sip_BC > 2.1361e-09])
        nrSIP_per_bin_floored = np.floor(nrSIP_per_bin_frac).astype(int) # at least nrSIP_per_bin_floored SIPs from a bin grid are generated
        prob_SIP = nrSIP_per_bin_frac - nrSIP_per_bin_floored # probability that nrSIP_per_bin_floored + 1 SIPs are generated from a bin
        p_vec = np.random.random(nr_SIPs_BC)
        nr_SIPs_save = nr_SIPs_tot
        for iSIP in range(nr_SIPs_BC):
            nr_SIPs_new = nrSIP_per_bin_floored[iSIP] + int(p_vec[iSIP]<prob_SIP[iSIP])
            for iSIP_per_bin in range(nr_SIPs_new):
                zEK_sip_ins[nr_SIPs_tot] = (nz * dz) - random.random() * zdiff[iSIP]
                nEK_sip_ins[nr_SIPs_tot] = nEK_sip_BC[iSIP]/nr_SIPs_skal
                mEK_sip_ins[nr_SIPs_tot] = mEK_sip_BC[iSIP]
                nr_SIPs_tot += 1
        FluxInAcctmp[0:4] += FK.Moments_k0_3(nEK_sip_ins[nr_SIPs_save:nr_SIPs_tot],mEK_sip_ins[nr_SIPs_save:nr_SIPs_tot])
        nr_SIPs_in = nr_SIPs_tot - nr_SIPs_save
        FluxInAcctmp[4] += nr_SIPs_in

        #blockwise output of all infalling SIPs
        if (nr_SIPs_in > 0):
            ntime_SIPin += 1
            np.savetxt(fGVin, np.array([int(t), nr_SIPs_in]).reshape(1,2), fmt='%6i')
            np.savetxt(fGVin, nEK_sip_ins[nr_SIPs_save:nr_SIPs_tot].reshape(1,nr_SIPs_in), fmt=outp_format)
            np.savetxt(fGVin, mEK_sip_ins[nr_SIPs_save:nr_SIPs_tot].reshape(1,nr_SIPs_in), fmt=outp_format)


        #GCCif (PROCESS == 1)
        # compute w of all new SIPs, only for onlySedi-case
        wEK_sip_ins[nr_SIPs_save:nr_SIPs_tot]=SD.Fallg(FK.m2r(mEK_sip_ins[nr_SIPs_save:nr_SIPs_tot],const_mass2rad))
        #GCCendif  /* (PROCESS == 1) */

        #Additional sort and SIP-GB-attribution required
        perm_sort=np.argsort(zEK_sip_ins)
        nEK_sip_ins=nEK_sip_ins[perm_sort]
        mEK_sip_ins=mEK_sip_ins[perm_sort]
        zEK_sip_ins=zEK_sip_ins[perm_sort]
        #GCCif (PROCESS != 2)
        wEK_sip_ins=wEK_sip_ins[perm_sort]
        #GCCendif /*(PROCESS != 2)*/
        [iSIP_GBsep,nr_SIPs_GB] = FK.get_isep(zEK_sip_ins,zGBsep,nz)
        #print('iSIP_GBsep, nr_SIPs_GB', iSIP_GBsep,nr_SIPs_GB)
        ausw_SIPs_active = np.flatnonzero(zEK_sip_ins <= zCol)
        nr_SIPs_tot=len(ausw_SIPs_active)

        #GCCendif /* (INFLUX_TOP == 1)*/

        #GCCif (TRACKCENTER == 1)
        if (it%imod_MOMsave ==0) & (it != 0):
            tmp,nom,denom=FK.trackcenters(nEK_sip_ins,mEK_sip_ins,zEK_sip_ins,iSIP_GBsep,nr_SIPs_GB,nz,dz) #returns array of size (2,nz,nr_MOMs)
            zCenterIC_full [k,i_MOMsave,:,:,0]=tmp[0,:,:]
            zCenterSIP_full[k,i_MOMsave,:,:,0]=tmp[1,:,:]
            nr_SIPs_GB_save[k,i_MOMsave,:]=nr_SIPs_GB
            zCenterIC_av[i_MOMsave,:,0]+=nom[0,:]
            zCenterIC_av_denom[i_MOMsave,:,0]+=denom[0,:]
            zCenterSIP_av[i_MOMsave,:,0]+=nom[1,:]
            zCenterSIP_av_denom[i_MOMsave,:,0]+=denom[1,:]
        #GCCendif  /* (TRACKCENTER == 1) */

        #GCCif (PROCESS != 2)
            #GCCif (PROCESS == 0)
        #update fall speeds of SIPs, not necessary when aggregation is switchd off
        rSIP_tmp=FK.m2r(mEK_sip_ins[ausw_SIPs_active],const_mass2rad)
        wEK_tmp=SD.Fallg(rSIP_tmp)
            #GCCendif /* (PROCESS == 0) */
            #GCCif (PROCESS == 1)
            #always take fall speeds computed during initialization
        wEK_tmp=wEK_sip_ins[ausw_SIPs_active]
        #print(wEK_tmp.shape,min(wEK_tmp),max(wEK_tmp))
            #GCCendif /* (PROCESS == 1) */
            #GCCif (WELLMIXED > 0)
        zEK_new=np.zeros(nr_SIPs_tot)
        zEK_new[ausw_SIPs_active]=zEK_sip_ins[ausw_SIPs_active] - wEK_tmp*dt
            #GCCendif /*(WELLMIXED > 0)*/
        #GCCendif /* (PROCESS != 2)*/

        #GCCif (PROCESS != 1)
            #GCCif (WELLMIXED <= 1) /* consider aggregation in each grid box separately */
        for iz in range(0,nz):
            if (nr_SIPs_GB[iz] > 1):
                ia=iSIP_GBsep[iz]
                ie=iSIP_GBsep[iz+1]
                #GCCif (WELLMIXED == 0)
                ibreak = AON.Aggregation(nr_SIPs_GB[iz],nEK_sip_ins[ia:ie],mEK_sip_ins[ia:ie],count_colls,cck,m_low,eta_indize,m_kernel,
                                         fLog_currColls=fLog_currColls, fLog_accColls=fLog_accColls, fLog_Combs=fLog_Combs, fLog_p=fLog_p)
                #GCCendif /* (WELLMIXED == 0)*/
                #GCCif (WELLMIXED == 1)
                ibreak = AON.Aggregation(nr_SIPs_GB[iz],nEK_sip_ins[ia:ie],mEK_sip_ins[ia:ie],zEK_sip_ins[ia:ie],zEK_new[ia:ie], count_colls,cck,m_low,eta_indize,m_kernel)
                #GCCendif /* (WELLMIXED == 1)*/
            #GCCendif /* (WELLMIXED <= 1) */
            #GCCif (WELLMIXED >= 2) /* consider aggregation in full column */
        ia=iSIP_GBsep[0]
        ie=iSIP_GBsep[nz]
        ibreak = AON.Aggregation(nr_SIPs_tot,nEK_sip_ins[ia:ie],mEK_sip_ins[ia:ie],zEK_sip_ins[ia:ie],zEK_new[ia:ie], count_colls,cck,m_low,eta_indize,m_kernel,
                                fLog_currColls=fLog_currColls, fLog_accColls=fLog_accColls, fLog_currColls_WM2D=fLog_currColls_WM2D, fLog_Combs=fLog_Combs, fLog_p=fLog_p)
            #GCCendif /* (WELLMIXED >= 2) */
        #GCCendif /* (PROCESS != 1)*/

        #GCCif (TRACKCENTER == 1)
        if (it%imod_MOMsave ==0) & (it != 0):
            tmp,nom,denom=FK.trackcenters(nEK_sip_ins,mEK_sip_ins,zEK_sip_ins,iSIP_GBsep,nr_SIPs_GB,nz,dz) #returns array of size (2,nz,nr_MOMs)
            zCenterIC_full [k,i_MOMsave,:,:,1]=tmp[0,:,:]
            zCenterSIP_full[k,i_MOMsave,:,:,1]=tmp[1,:,:]
            zCenterIC_av[i_MOMsave,:,1]+=nom[0,:]
            zCenterIC_av_denom[i_MOMsave,:,1]+=denom[0,:]
            zCenterSIP_av[i_MOMsave,:,1]+=nom[1,:]
            zCenterSIP_av_denom[i_MOMsave,:,1]+=denom[1,:]
        #GCCendif  /* (TRACKCENTER == 1) */

        #GCCif (PROCESS != 2)
        zEK_sip_ins[ausw_SIPs_active]=zEK_sip_ins[ausw_SIPs_active] - wEK_tmp*dt
        #GCCendif /* (PROCESS != 2)*/

        #GCCif (RANDOMZ == 1)
        zEK_sip_ins[ausw_SIPs_active] = np.random.random(nr_SIPs_tot)*nz*dz
        #GCCendif /* (RANDOMZ == 1)*/

        #GCCif (TRACKCENTER == 1)
        if (it%imod_MOMsave ==0) & (it != 0):
            tmp,nom,denom=FK.trackcenters(nEK_sip_ins,mEK_sip_ins,zEK_sip_ins,iSIP_GBsep,nr_SIPs_GB,nz,dz) #returns array of size (2,nz,nr_MOMs)
            zCenterIC_full [k,i_MOMsave,:,:,2]=tmp[0,:,:]
            zCenterSIP_full[k,i_MOMsave,:,:,2]=tmp[1,:,:]
            np.savetxt(fCenterIC_full,  zCenterIC_full[k,i_MOMsave,:,:,:].reshape(nz,4*3), fmt='%1.4f')
            np.savetxt(fCenterSIP_full, zCenterSIP_full[k,i_MOMsave,:,:,:].reshape(nz,4*3), fmt='%1.4f')
            zCenterIC_av[i_MOMsave,:,2]+=nom[0,:]
            zCenterIC_av_denom[i_MOMsave,:,2]+=denom[0,:]
            zCenterSIP_av[i_MOMsave,:,2]+=nom[1,:]
            zCenterSIP_av_denom[i_MOMsave,:,2]+=denom[1,:]
        #GCCendif  /* (TRACKCENTER == 1) */

        #--------update SIP list---------
            #-------remove zero weight SIPs--------
        ausw_zeroweights=np.flatnonzero(nEK_sip_ins <= 1e-50)
        ntmp=len(ausw_zeroweights)
        if (ntmp > ntmp_before):
            print("it = ",it,"ntmp = ", ntmp)
            print(nEK_sip_ins[ausw_zeroweights])
            print(zEK_sip_ins[ausw_zeroweights])
        zEK_sip_ins[ausw_zeroweights] = zNan
        #zEK_sip_ins[0:nr_SIPs_tot][nEK_sip_ins[ausw_SIPs_active] <= 1e-50] = zNan
        ntmp_before=ntmp
            #---------treatment of lower boundary condition-------------
        ausw_outf = np.flatnonzero(zEK_sip_ins <= 0.0)
        ntmp=len(ausw_outf)
        if (ntmp > 0):
            FluxOutAcctmp[0:4]=FluxOutAcctmp[0:4]+FK.Moments_k0_3(nEK_sip_ins[ausw_outf],mEK_sip_ins[ausw_outf])
            FluxOutAcctmp[4]=FluxOutAcctmp[4]+ntmp
            #GCCif (TRACKOUT == 1)
            ntime_SIPout += 1
            np.savetxt(fGVout, np.array([int(t),ntmp]).reshape(1,2), fmt='%6i')
            np.savetxt(fGVout, nEK_sip_ins[ausw_outf].reshape(1,ntmp), fmt=outp_format)
            np.savetxt(fGVout, mEK_sip_ins[ausw_outf].reshape(1,ntmp), fmt=outp_format)
            #GCCendif  /* (TRACKOUT == 1) */

            #GCCif (INFLUX_TOP == 2)
                #---periodic boundary conditions; insert outfalling SIPs at top boundary-------
            FluxInAcctmp=FluxOutAcctmp    #in rare cases of very heavy SIPs they might cross the lower BC twice during a single timestep, i.e. modulo operation substracts 2*zCol. Flux evaluation does not account for this rare case.
            zEK_sip_ins[ausw_outf]=zEK_sip_ins[ausw_outf]%zCol
            #GCCelse
                #---------remove outfalling SIPs--------------
            zEK_sip_ins[ausw_outf]=zNan
            #GCCendif /* (INFLUX_TOP == 2)*/

        #---------sort SIPs by z-position and find index range of each grid box----------------
        perm_sort=np.argsort(zEK_sip_ins)
        nEK_sip_ins=nEK_sip_ins[perm_sort]
        mEK_sip_ins=mEK_sip_ins[perm_sort]
        zEK_sip_ins=zEK_sip_ins[perm_sort]
        #GCCif (PROCESS != 2)
        wEK_sip_ins=wEK_sip_ins[perm_sort]
        #GCCendif /*(PROCESS != 2)*/

        #--------attribute SIPs to grid boxes----------
        [iSIP_GBsep,nr_SIPs_GB] = FK.get_isep(zEK_sip_ins,zGBsep,nz)

        ausw_SIPs_active = np.flatnonzero(zEK_sip_ins <= zCol)
        nr_SIPs_tot = len(ausw_SIPs_active)
        #print('nr_SIPs_tot ', nr_SIPs_tot,len(perm_sort))
        #>>>>>>>>>>save and output moment data>>>>>>>>>>>
        if (it%imod_MOMsave ==0) & (it != 0):
            for iz in range(0,nz):
                nr_SIPs = nr_SIPs_GB[iz]
                ia=iSIP_GBsep[iz]
                ie=iSIP_GBsep[iz+1]
                if (nr_SIPs > 0): MOMsave[k,i_MOMsave,iz,:] = FK.Moments_k0_3(nEK_sip_ins[ia:ie],mEK_sip_ins[ia:ie])
                np.savetxt(fMom, MOMsave[k,i_MOMsave,iz,:].reshape(1,4), fmt=outp_format_long)# reshape necessary to write all 4 values in one row
            i_MOMsave = i_MOMsave+1
        #<<<<<<<<<<save moment data<<<<<<<<<<<<

        #>>>>>>>>>> output SIP & save and output flux data>>>>>>>>>>-
        if (it%imod_GVplot == 0) & (it != 0):
            #---------output SIP data>>>>>>>>>>---------
            print('it,it*dt,nr_SIPs_tot',it,it*dt,nr_SIPs_tot)
            #print('iSIP_GBsep', iSIP_GBsep)
            #print('n[:5] ', nEK_sip_ins[:5])
            if (iPM >= 2): print('SIP-Daten rausschreiben, #,Iter-schritt, Zeit:', i_GVplot, it, it*dt)
            for iz in range(0,nz):
                nr_SIPs = nr_SIPs_GB[iz]
                ia=iSIP_GBsep[iz]
                ie=iSIP_GBsep[iz+1]
                np.savetxt(fGV, np.array([nr_SIPs,iz]).reshape(1,2), fmt='%6i')
                if (nr_SIPs > 0):
                    np.savetxt(fGV, nEK_sip_ins[ia:ie].reshape(1,nr_SIPs), fmt=outp_format)
                    np.savetxt(fGV, mEK_sip_ins[ia:ie].reshape(1,nr_SIPs), fmt=outp_format)
                    np.savetxt(fGV, zEK_sip_ins[ia:ie].reshape(1,nr_SIPs), fmt=outp_format)

            fGVMeta = open('SIP_meta.dat', 'a')
            fGVMeta.write("{:7d} {:5d}\n".format(nr_SIPs_GB.sum(),nr_SIPs_GB.max()))
            fGVMeta.close()

            #------output of influx/outflux----------
            skal=np.array([dAi,dAi,dAi,dAi,1])
            FluxInAcc[k,i_GVplot-1,:] = FluxInAcctmp*skal
            FluxOutAcc[k,i_GVplot-1,:]= FluxOutAcctmp*skal
            np.savetxt(fFluxInAcc, FluxInAcctmp.reshape(1,5)*skal, fmt=outp_format_flux)
            np.savetxt(fFluxOutAcc, FluxOutAcctmp.reshape(1,5)*skal, fmt=outp_format_flux)

            i_GVplot = i_GVplot+1

    #>>>>>>end of time iteration of a single instance>>>>>>>>>

    ## check influx SIP condition
    ## see "Sims1D/Tests/Test_11"
    #ia=iSIP_GBsep[0]
    #ie=iSIP_GBsep[nz]
    #ibin_eachSIP = np.ceil(np.log10(mEK_sip_ins[ia:ie]/mfix[0])*n10).astype(int)
    #nbin_allSIPs = np.zeros(50, dtype='int')
    #for iii in ibin_eachSIP:
        #nbin_allSIPs[iii] += 1
    #print('ibin_eachSIP: ', ibin_eachSIP)
    #print('nbin_allSIPs: ', nbin_allSIPs)
    #zmaxbin_allSIPs = np.zeros(50)
    #for iii,nbin_tmp in enumerate(nbin_allSIPs):
        #if (nbin_tmp > 0):
            #lowestGB=np.floor(np.min(zEK_sip_ins[ia:ie][ibin_eachSIP == iii]/dz)).astype(int)
            #occupiedGBs=nz-lowestGB
            #print("{:3d} {:5d} {:5d} {:5d} {:.2}".format(iii, nbin_tmp, lowestGB,occupiedGBs,((nbin_tmp/occupiedGBs)*100)//1/100))

    #GCCif (TRACKOUT == 1)
    fGVoutMeta  = open('SIPout_meta.dat', 'a')
    fGVoutMeta.write("{}\n".format(ntime_SIPout))
    fGVoutMeta.close()
    #GCCendif  /* (TRACKOUT == 1) */
    #GCCif (INFLUX_TOP == 1)
    fGVinMeta  = open('SIPin_meta.dat', 'a')
    fGVinMeta.write("{}\n".format(ntime_SIPin))
    fGVinMeta.close()
    #GCCendif /* (INFLUX_TOP == 1)*/


    #>>>>>>>>>>analyse computing time>>>>>>>>>>>>>>>>>>>>
    currenttime = time.time()
    currenttime_str = time.asctime( time.localtime(currenttime))
    endtime_expected=starttime+ (nr_inst/(k+1)) * (currenttime-starttime)
    endtime_expected_str = time.asctime( time.localtime(endtime_expected))
    if (iPM >= 2):
        print("Instance {} of {} finished".format(k+1,nr_inst))
        print("Start time/Current time/Expected end time: ")
    print(localtime, ' --- ', currenttime_str, ' --- ', endtime_expected_str)
    fLog= open('log.txt','a')
    fLog.write('total computing time in sec: '+ str(int(currenttime-starttime)) + '\n')
    fLog.write(currenttime_str+ ' --- '+ endtime_expected_str+ '\n')

#GCCif (COUNT_COLLS == 1)
    cc_sum=count_colls.sum()
    fLog.write('a '+"{}".format(cc_sum)+'\n')
    fLog.write('b '+" ".join("{}".format(x) for x in count_colls)+'\n')
    cc_frac=count_colls/cc_sum
    fLog.write('b '+" ".join("{:.4}".format(x) for x in cc_frac)+'\n')
#GCCendif  /* (COUNT_COLLS == 1) */

    fLog.close()
#<<<<<<<<<<<<<<<<<end of loop of all realisations<<<<<<<<<<<<<<<<<
#---------------------------------------------------------------------------------------------

#GCCif (TRACKCENTER == 1)
#Compute normalised IC and SIP centers, averaged over column and instances
zCenterIC_av=zCenterIC_av/zCenterIC_av_denom
zCenterSIP_av=zCenterSIP_av/zCenterSIP_av_denom
np.savetxt(fCenterIC_av,  zCenterIC_av.reshape(nr_MOMsave,4*3), fmt='%1.4f')
np.savetxt(fCenterSIP_av, zCenterSIP_av.reshape(nr_MOMsave,4*3), fmt='%1.4f')
fCenterIC_full.close()
fCenterSIP_full.close()
fCenterIC_av.close()
fCenterSIP_av.close()
#GCCendif  /* (TRACKCENTER == 1) */

#Diagnose "instantaneous" fluxes from accumulated fluxes
for i_inst in range(nr_inst):
    fo_tm1=np.zeros(5) # flux out at "t minus 1"
    fi_tm1=np.zeros(5) # flux in  at "t minus 1"
    skal=np.array([t_intervall_GVplot]*5)
    skal[4]=1
    for i_GVplot in range(nr_GVplot-1):
        fo=FluxOutAcc[i_inst,i_GVplot,:]
        fi=FluxInAcc[i_inst,i_GVplot,:]
        FluxOut[i_inst,i_GVplot,:]= (fo - fo_tm1)/skal
        FluxIn[i_inst,i_GVplot,:] = (fi - fi_tm1)/skal
        fo_tm1=fo
        fi_tm1=fi
        np.savetxt(fFluxIn, FluxIn[i_inst,i_GVplot,:].reshape(1,5), fmt=outp_format_flux)
        np.savetxt(fFluxOut, FluxOut[i_inst,i_GVplot,:].reshape(1,5), fmt=outp_format_flux)

fFluxInAcc.close()
fFluxOutAcc.close()
fFluxIn.close()
fFluxOut.close()

fMom.close()
fGV.close()

#GCCif (TRACKOUT == 1)
fGVout.close()
#GCCendif  /* (TRACKOUT == 1) */
#GCCif (INFLUX_TOP == 1)
fGVin.close()
#GCCendif /* (INFLUX_TOP == 1)*/

localtime = time.asctime( time.localtime(time.time()) )
print("End time computation:", localtime)

fLog= open('log.txt','a')
currenttime = time.time()
fLog.write('total computing time in sec: '+ str(int(currenttime-starttime)) + '\n')
fLog.write('finalised')
fLog.close()

#GCCif (MOM_meanTD >= 1)
    #Generate output file with mean moments, currently the average is taken over all instances and grid boxes
    #plot mean moments
MOMmean=FK.CIO_MOMmean(data_in=MOMsave,skal_m=skal_m,dV=dV,fp_out='')
    #GCCif (MOM_meanTD == 2)
MOMmean = MOMmean * dz * nz
    #GCCendif /* (MOM_meanTD == 2) */
PS.PlotMoments(MOMmean,t_vec_MOMsave,iMean=1)
#GCCendif /* (MOM_meanTD >= 1) */

#GCCif (GV_meanTD == 1 || RZ_scatter == 1)
    #Process SIP data, read data from file
    #read all SIPs of one instance in a single long list (i.e neglect the vertical information)
        #read SIP meta data
fGVMeta  = open('SIP_meta.dat', 'r')
nr_inst   = int(fGVMeta.readline())
nr_GVplot = int(fGVMeta.readline())
t_vec_GVplot= np.array(fGVMeta.readline().split(), dtype='float' )
nr_SIPs_plot=np.zeros([nr_inst,nr_GVplot],dtype='int')
nr_SIPs_max=0 #total SIP number over full column, maximum over time and instances
nr_SIPs_GB_maxmax=0 # nr_SIPs_GB_max: highest SIP number in a single grid box; nr_SIPs_GB_maxmax: take maximum over time and instances
for i_inst in range(nr_inst):
    fGVMeta.readline()
    for i_GVplot in range(nr_GVplot):
        #print(i_inst,i_GVplot)
        [nr_SIPs,nr_SIPs_GB_max]=np.array(fGVMeta.readline().split(), dtype='int')
        #print(nr_SIPs,nr_SIPs_GB_max)
        nr_SIPs_plot[i_inst,i_GVplot ]= nr_SIPs
        nr_SIPs_max=max([nr_SIPs,nr_SIPs_max])
        nr_SIPs_GB_maxmax=max([nr_SIPs_GB_max,nr_SIPs_GB_maxmax]) # currently not used in the following evaluation
print('nr_SIPs_max,nr_SIPs_GB_maxmax',nr_SIPs_max,nr_SIPs_GB_maxmax)
fGVMeta.close()
        #read actual SIP data
nEK_sip_plot=np.zeros([nr_inst,nr_GVplot,nr_SIPs_max])
mEK_sip_plot=np.zeros([nr_inst,nr_GVplot,nr_SIPs_max])
zEK_sip_plot=np.zeros([nr_inst,nr_GVplot,nr_SIPs_max])
fGV  = open('SIP.dat', 'r')
for i_inst in range(nr_inst):
    for i_GVplot in range(nr_GVplot):
        iSIP=0
        for k in range(nz):
            [nr_SIPs,iz]=np.array(fGV.readline().split(),dtype='int')
            #print('i_GVplot,nr_SIPs,iz',nr_SIPs,iz)
            if (nr_SIPs > 0):
                ia=iSIP
                ie=ia+nr_SIPs
                #print('ia,ie', ia,ie)
                nEK_sip_plot[i_inst,i_GVplot,ia:ie]=np.array(fGV.readline().split())
                mEK_sip_plot[i_inst,i_GVplot,ia:ie]=np.array(fGV.readline().split())
                #fGV.readline()
                zEK_sip_plot[i_inst,i_GVplot,ia:ie]=np.array(fGV.readline().split())
                iSIP=iSIP+nr_SIPs
        #print('-----------------')
        #print(nEK_sip_plot[i_inst,i_GVplot,:])wEK_sip_ins[:nr_SIPs_tot]
fGV.close()

#GCCif (GV_meanTD == 1)
PS.PlotGV(nEK_sip_plot,mEK_sip_plot,nr_SIPs_plot,t_vec_GVplot)
#GCCendif
#GCCif (RZ_scatter == 1 )
PS.Scatter_r_z(zEK_sip_plot,mEK_sip_plot,nr_SIPs_plot,t_vec_GVplot)
#GCCendif
#GCCendif  /* (GV_meanTD == 1 || RZ_scatter == 1)  */

#GCCif (MOM_prof == 1)
    #Generate output file with mean moments, currently the average is taken over all instances and grid boxes
    #plot mean moments
MOMmean_prof=FK.CIO_MOMmean(data_in=MOMsave,skal_m=skal_m,dV=dV,ikeep1D=1)
Times=[int(i) for i in np.arange(0,Tsim+1,600)]
PS.PlotMomentsProf(MOMmean_prof,t_vec_MOMsave,iMean=1,iTimes_select=Times)
#GCCendif /* (MOM_prof == 1) */

#GCCif (FLUX_time == 1)
PS.PlotFluxesTime(FluxIn,FluxOut,FluxInAcc,FluxOutAcc,t_vec_GVplot[1:],iplot_onlyMean=1)
#GCCendif /* (FLUX_time == 1) */

#GCCif (TRACKCENTER == 1)
#PS.PlotCenter(zCenterIC_full,nr_SIPs_GB_save,t_vec_MOMsave, 'CenterIC')
#PS.PlotCenter(zCenterSIP_full,nr_SIPs_GB_save,t_vec_MOMsave, 'CenterSIP')
nr_SIPs_acc=nr_SIPs_GB_save.sum(axis=2).sum(axis=0)
PS.PlotCenter(zCenterIC_av,nr_SIPs_acc,t_vec_MOMsave, 'CenterIC', av=1)
PS.PlotCenter(zCenterSIP_av,nr_SIPs_acc,t_vec_MOMsave, 'CenterSIP',av=1)
#GCCendif  /* (TRACKCENTER == 1) */



# end of column model section
#=============================== COLUMN-MODEL =============================================
#GCCendif /* (COLUMN == 1) */

