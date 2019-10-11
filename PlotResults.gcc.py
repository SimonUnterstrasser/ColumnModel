import sys
import os
import math
import numpy as np
import InputOutput as IO
import PlotSim as PS
#GCCif (BIN == 1)
import Referenzloesung as REF
#GCCendif /* (BIN == 1) */

text=None
label=None
title=None
ilastGVonly=None
iDmean=0
#GCCinclude "params_plot.txt"

#GCCif (((IREF == 3) && (KERNEL != 3)) || ((IREF == 2) && (KERNEL != 2)))
print('check preprocessor settings')
print('wrong combination of IREF and KERNEL')
exit(1)
#GCCendif
#GCCif (((IREF == 3) && (DISCRETE == 0)) || ((IREF == 2) && (DISCRETE == 0)))
print('check preprocessor settings')
print('wrong combination of IREF and DISCRETE')
exit(1)
#GCCendif

#GCCif ((MOM_meanTD == 2) && (COLUMN == 0)) 
print('check preprocessor settings')
print('wrong combination of MOM_meanTD and COLUMN')
exit(1)
#GCCendif


#GCCif (IMODE == 1)
def plotSingleSimulationGV(fp_in):
    #Einlesen der SIP-Daten
    nEK_sip_plot, mEK_sip_plot, zEK_sip_plot, nr_SIPs_plot, nr_SIPs_prof, t_vec_GVplot, V, skal_m = IO.readSingleSimulationSIPdata(fp_in)
    #[nEK_sip_plot, mEK_sip_plot, nr_SIPs_plot, t_vec_GVplot] = IO.readSingleSimulationSIPdataMARION(fp_in)
    #print(nr_SIPs_plot)
    nz=1
    #GCCif (COLUMN == 1)
    nz=IO.get_nz(fp_in)
    #GCCendif /*(COLUMN == 1)*/
    #Erstelle GV-Plot
    PS.PlotGV(nEK_sip_plot/nz, mEK_sip_plot, nr_SIPs_plot, t_vec_GVplot, ilastGVonly = ilastGVonly)
#<<<<<<<< end routine plotSingleSimulationGV(fp)

    #GCCif (DISCRETE  >= 1)
def plotSingleSimulationRStD(fp_in):
#Einlesen der SIP-Daten
    nEK_sip_plot, mEK_sip_plot,zDummy, nr_SIPs_plot, nr_SIPs_prof, t_vec_GVplot, V, skal_m = IO.readSingleSimulationSIPdata(fp_in)
#Erstelle Plot von relativer Standardabweichung
    PS.PlotRelStdDev(mEK_sip_plot, t_vec_GVplot)
#<<<<<<<< end routine plotSingleSimulationRStD(fp)
    #GCCendif /* (DISCRETE  >= 1) */

def plotSingleSimulationMOM(fp_in,iprof=0,itot=1):
# reads Moment data of box model or column model data
# for box model the design of the plot is straighforward, i.e. temporal evolution of moments
# for the colum model several options exists:
#    1. temporal evolution of column integrated moments (itot = 1)
#    2. profiles of moments at selected times (iprof = 1)

    MOMsave, t_vec_MOMsave, skal_m, dV=IO.readSingleSimulationMoments(fp_in)

    iSIPplot = 0
    nr_SIPs_prof = None
    t_vec_GVplot = None
    #GCCif (addSIPplot > 0)
    iSIPplot = 1
    nEK_sip_plot, mEK_sip_plot,zDummy, nr_SIPs_plot, nr_SIPs_prof, t_vec_GVplot, V, skal_m = IO.readSingleSimulationSIPdata(fp_in)
    nr_SIPs_mean = np.mean(nr_SIPs_plot, axis=0)
    nr_SIPs_prof_mean = np.mean(nr_SIPs_prof, axis=0)
    #GCCendif /* (addSIPplot > 0) */

    #GCCif (STD == 0)
    iplotStdDev = 0
    #GCCendif /* (STD == 0) */
    #GCCif (STD == 1)
    iplotStdDev = 1
    #GCCendif /* (STD == 1) */
    #GCCif (STD == 2)
    iplotStdDev = 2
    #GCCendif /* (STD == 2) */

    nz = 1  # for COLUMN = 0

    #GCCif (COLUMN == 1)
    #in column model version: 1. call profile plot 2. integrate over column 3. call total moments plot

    # >>>> 1. call profile plot
    dz=np.zeros(1)
    dz[0] = IO.get_dz(fp_in)
    nz=np.zeros(1, dtype=int)
    nz[0] = IO.get_nz(fp_in)

    if (iprof == 1):
        if (iSIPplot == 1):
            #GCCif (addSIPplot == 1)
            nr_SIPs_prof_mean = nr_SIPs_prof_mean*100/dz[0]
            #GCCendif /* (addSIPplot == 1) */
            #GCCif (addSIPplot == 2)
            nr_SIPs_prof_mean = nr_SIPs_prof_mean
            #GCCendif /* (addSIPplot == 2) */
        for Times_elem in Times_list:
            PS.PlotMomentsProf(MOMsave, t_vec_MOMsave, iplot_onlyMOMmean=1, 
                               nz_vec=nz, dz_vec=dz, iplotStdDev=iplotStdDev, iDmean=iDmean,
                               iTimes_select=Times_elem, 
                               iSIPplot=iSIPplot, nr_SIPs_prof_mean=nr_SIPs_prof_mean, t_vec_GVplot=t_vec_GVplot)

    # >>>> 2. integrate over column
    #nz=IO.get_nz(fp_in)
    MOMsave=np.mean(MOMsave,axis=2)
    #GCCendif /*(COLUMN == 1)*/

    #GCCif (MOM_meanTD == 2)
    MOMsave = MOMsave * dz[0] * nz[0]
    #GCCendif /* (MOM_meanTD == 2) */
    #GCCif (MOM_meanTD == 3)
    MOMsave = MOMsave * dV * nz[0]
    #GCCendif /* (MOM_meanTD == 3) */

    #Erstelle Momenten-Plot
    if (itot == 1): 
        #GCCif (MOM_skal_m == 0)
        PS.PlotMoments(MOMsave,t_vec_MOMsave,iplot_onlyMOMmean=1,iplotStdDev=iplotStdDev,iDmean=iDmean,
                       iSIPplot=iSIPplot,nr_SIPs_mean=nr_SIPs_mean,t_vec_GVplot=t_vec_GVplot)
        #GCCelse
        PS.PlotMoments(MOMsave,t_vec_MOMsave,iplot_onlyMOMmean=1,iplotStdDev=iplotStdDev,iDmean=iDmean,
                       skal_m=skal_m,
                       iSIPplot=iSIPplot,nr_SIPs_mean=nr_SIPs_mean,t_vec_GVplot=t_vec_GVplot)
        #GCCendif /* else (MOM_skal_m == 0) */
#<<<<<<<< end routine plotSingleSimulationMOM(fp)

def plotSingleSimulationFLUX(fp_in):
    #GCCif (STD == 0)
    iplotStdDev = 0
    #GCCendif /* (STD == 0) */
    #GCCif (STD == 1)
    iplotStdDev = 1
    #GCCendif /* (STD == 1) */
    #GCCif (STD == 2)
    iplotStdDev = 2
    #GCCendif /* (STD == 2) */
    #read flux data
    FluxIn,FluxOut,FluxInAcc,FluxOutAcc,t_vec_GVplot=IO.readSingleSimulationFluxes(fp_in)
    #generate flux plot
    PS.PlotFluxesTime(FluxIn,FluxOut,FluxInAcc,FluxOutAcc,t_vec_GVplot[1:],iplot_onlyMean=1,iplotStdDev=iplotStdDev)
#<<<<<<<< end routine plotSingleSimulationFLUX(fp_in)

def plotSingleSimulationGVout(fp_in,infalling=0):
    # Einlesen der SIP-Daten
    # returns a list of SIPs of length nr_SIPs_out with nu and m and the time tEK_sip_out each SIP crossed the lower BC
    # A is the basal area of the column
    nEK_sip_out, mEK_sip_out, tEK_sip_out, nr_SIPs_out, Tsim, A = IO.readSingleSimulationSIP_outfalling_data(fp_in,infalling)

    PS.PlotGV(nEK_sip_out/A, mEK_sip_out, nr_SIPs_out, [Tsim], outfallingGV = 1 + infalling)
#<<<<<<<< end routine plotSingleSimulationGV(fp)
#GCCendif /* (IMODE == 1) */

#GCCif (IMODE == 2)
import Misc as FK
import PlotSim as PS

def plotMultipleSimulationsMOM(filepath_input,iTime=0,iProf=0):
    nr_sims=len(filepath_input)
    print('call of plotMultipleSimulationMOM')
    print('nr_sims:',nr_sims)
    #if (label is not None) and (len(label) != nr_sims+iaddsimREF): 
        #print('wrong number of labels are given')
        #print('label',label)
        #print('nr_sims, iaddsimREF: ', nr_sims, iaddsimREF)
    iplot_mode=np.zeros(nr_sims,dtype=int) +2
    iplot_mode[0]=1
    iplot_mode[-1]=3

    #GCCif (STD == 0)
    iplotStdDev = 0
    MOM_StdDev = None
    MOMprof_StdDev = None
    #GCCendif /* (STD == 0) */
    #GCCif (STD == 1)
    iplotStdDev = 1
    #GCCendif /* (STD == 1) */
    #GCCif (STD == 2)
    iplotStdDev = 2
    #GCCendif /* (STD == 2) */

    nt_vec      = np.zeros(nr_sims,dtype='int')
    nr_inst_vec = np.zeros(nr_sims,dtype='int')
    dV_vec      = np.zeros(nr_sims)
    nz_vec      = np.zeros(nr_sims,dtype='int')

    #GCCif (COLUMN == 1)
    dz_vec      =np.zeros(nr_sims)
    #find out first the maximum dimensions of the profile data
    for i in range(nr_sims):
        fp=filepath_input[i]
        print('i, fp', i, fp)
        #GCCif (BIN != 1)
        nz,dz,Tsim,nr_inst,LWC,r0,xf0,ikernel,i_init_1D,i_process=IO.get_MetaData(fp)
        nz_vec[i]      = nz
        dz_vec[i]      = dz
        nr_inst_vec[i] = nr_inst
        #####GCCelse
        ####ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF,rgridREF,mgridREF=REF.get_RefMetaData(fp=fp)
        ####nz_vec[i]      = nzREF
        ####dz_vec[i]      = dzREF
        ####nt_vec[i]      = ntREF
        ####nr_inst_vec[i] = 1
        #GCCendif /* (BIN != 1) */
    if (iProf == 1):
        nz_max      = nz_vec.max()
        dz_max      = dz_vec.max()
        nr_inst_max = nr_inst_vec.max()

        nt_max=300
        nr_Mom = 4
        MOMsave=np.zeros([nr_sims,nt_max,nz_max,nr_Mom])
        t_vec_MOM_vec=np.zeros([nr_sims,nt_max])
        if (iplotStdDev == 1):
            MOMprof_StdDev = np.zeros([1,nr_sims,nt_max,nz_max,nr_Mom])
        if (iplotStdDev == 2):
            MOMprof_StdDev = np.zeros([2,nr_sims,nt_max,nz_max,nr_Mom])

    iSIPplot = 0
    nr_SIPs_prof_mean = None
    nr_SIPs_mean = None
    t_vec_GVplot_vec = None
    t_vec_GVplotsingleSim = None
    ntGV_vec = None
    #GCCif (addSIPplot > 0)
    iSIPplot = 1
    nz_max      = nz_vec.max()
    nr_SIPs_mean, nr_SIPs_prof_mean, t_vec_GVplot_vec, ntGV_vec = plotMultipleSimulationGV(filepath_input, igetSIPdataonly=1, nz_max=nz_max, dz_vec=dz_vec) # uebergibt SIP-Daten aller Sims
    t_vec_GVplotsingleSim = np.squeeze(t_vec_GVplot_vec[0,:])
    #GCCendif /* (addSIPplot > 0) */

    #GCCendif /*(COLUMN == 1)*/

    for i in range (nr_sims):
        fp=filepath_input[i]
        print('process simulation:',fp)
        #GCCif (BIN != 1)
        MOMsingleSim, t_vec_MOMsingleSim, skal_m, dV =IO.readSingleSimulationMoments(fp)
        dV_vec[i]      = dV
        #GCCelse
        #nzREF = nz_vec[i]
        #ntREF = nt_vec[i]
        ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF,rgridREF,mgridREF=REF.get_RefMetaData(fp=fp)

        MOMsingleSim = REF.get_RefProfileData(nzREF, ntREF, fp, iswap=1, fp=fp)
        t_vec_MOMsingleSim = TsimREF/(ntREF-1)*np.arange(ntREF)
        dV_vec[i]      = 1.0
        skal_m = 1
        #GCCendif /* (BIN != 1) */

        #GCCif (COLUMN == 0)
        nz_vec[i]      = 1
        #GCCendif /* (COLUMN == 0) */

        if (iProf == 1):
            nt_vec[i] = t_vec_MOMsingleSim.size
            print('AA',nt_vec[i])
            if (nt_vec[i] > nt_max): 
                print ('nt_vec[i], nt_max', nt_vec[i], nt_max)
                sys.exit('increase nt_max!!')
            #GCCif (STD > 0)
            MOMsave[i,:nt_vec[i],:nz_vec[i],:], MOMprof_StdDev[:,i,:nt_vec[i],:nz_vec[i],:] = \
                FK.CIO_MOMmean(data_in=MOMsingleSim,ikeep1D=1,igetSTD=iplotStdDev)
            #GCCelse
            MOMsave[i,:nt_vec[i],:nz_vec[i],:]=FK.CIO_MOMmean(data_in=MOMsingleSim,ikeep1D=1)
            #GCCendif /* else (STD > 0) */

            t_vec_MOM_vec[i,:nt_vec[i]]=t_vec_MOMsingleSim

        if (iTime == 1):
            #produce plot MeanMoments over time
            #individual call for each simulation

            iexcludetop = 0
            #GCCif (STD > 0)
            MOMmean, MOM_StdDev = FK.CIO_MOMmean(data_in=MOMsingleSim,igetSTD=iplotStdDev,iexcludetop=iexcludetop)
            #GCCelse
            print(MOMsingleSim.shape)
            MOMmean=FK.CIO_MOMmean(data_in=MOMsingleSim,iexcludetop=iexcludetop)
            print(MOMmean.shape)
            #GCCendif /* else (STD > 0) */

            print('sim: ', i)
            print('initial moments ', MOMmean[0])
            print('skal_m: ', skal_m)

            #GCCif (MOM_meanTD == 2)
            MOMmean = MOMmean * dz_vec[i] * nz_vec[i]
            #GCCif (STD > 0)
            MOM_StdDev = MOM_StdDev * dz_vec[i] * nz_vec[i]
            #GCCendif /* (STD > 0) */
            #GCCendif /* (MOM_meanTD == 2) */

            #GCCif (MOM_meanTD == 3)
            MOMmean = MOMmean * dV_vec[i] * nz_vec[i]
            #GCCif (STD > 0)
            MOM_StdDev = MOM_StdDev * dV_vec[i] * nz_vec[i]
            #GCCendif /* (STD > 0) */
            #GCCendif /* (MOM_meanTD == 3) */

            nr_SIPs_mean_pick = None
            #GCCif (addSIPplot > 0)
            nr_SIPs_mean_pick = nr_SIPs_mean[i,:]
            #GCCendif /* (addSIPplot > 0) */

            #GCCif (MOM_skal_m == 0)
            PS.PlotMoments(MOMmean,t_vec_MOMsingleSim,iMean=1,iplot_mode=iplot_mode[i],
                           label=label,title=title,nr_sims=nr_sims,text=text,
                           iplotStdDev=iplotStdDev, MOM_StdDev=MOM_StdDev, iDmean=iDmean,
                           iSIPplot=iSIPplot,nr_SIPs_mean=nr_SIPs_mean_pick,t_vec_GVplot=t_vec_GVplotsingleSim,
                           iexcludetop=iexcludetop)
            #GCCelse
            PS.PlotMoments(MOMmean,t_vec_MOMsingleSim,iMean=1,iplot_mode=iplot_mode[i],
                           label=label,title=title,nr_sims=nr_sims,text=text,
                           iplotStdDev=iplotStdDev, MOM_StdDev=MOM_StdDev, iDmean=iDmean, skal_m=skal_m,
                           iSIPplot=iSIPplot,nr_SIPs_mean=nr_SIPs_mean_pick,t_vec_GVplot=t_vec_GVplotsingleSim,
                           iexcludetop=iexcludetop)
            #GCCendif /* else (MOM_skal_m == 0) */

#produce plot MomentProfiles
    if (iProf == 1):
        for Times_elem in Times_list: 
            PS.PlotMomentsProf(MOMsave, t_vec_MOM_vec, iMultipleSims=1, iMean=1, 
                               nr_inst_vec=nr_inst_vec, nz_vec=nz_vec, nt_vec=nt_vec, dz_vec = dz_vec,
                               label=label,title=title,text=text,
                               iplot_onlyMOMmean=1,iTimes_select=Times_elem,
                               iplotStdDev=iplotStdDev, MOMprof_StdDev=MOMprof_StdDev, iDmean=iDmean, 
                               iSIPplot=iSIPplot, nr_SIPs_prof_mean=nr_SIPs_prof_mean, t_vec_GVplot=t_vec_GVplot_vec, ntGV_vec=ntGV_vec
                               )
#<<<<<<<< end routine  plotMultipleSimulationsMOM(filepath_input)

def plotMultipleSimulationGV(filepath_input, igetSIPdataonly=0, nz_max=0, dz_vec=0):
    nr_sims=len(filepath_input)
    print('call of plotMultipleSimulationGV')
    print('nr_sims:',nr_sims)

    nEK_sip_plot=np.zeros([nr_sims,nr_inst_max,nr_GVplot_max,nr_sip_max])
    mEK_sip_plot=np.zeros([nr_sims,nr_inst_max,nr_GVplot_max,nr_sip_max])
    nr_SIPs_plot=np.zeros([nr_sims,nr_inst_max,nr_GVplot_max],dtype='int')

    #GCCif (addSIPplot > 0)
    if (igetSIPdataonly == 1):
        #nz_max = 1000
        nr_SIPs_av=np.zeros([nr_sims,nr_GVplot_max])
        nr_SIPs_prof_av=np.zeros([nr_sims,nr_GVplot_max,nz_max])
        t_vec_GVplot_vec=np.zeros([nr_sims,nr_GVplot_max])
        ntGV_vec=np.zeros([nr_sims],dtype='int')
    #GCCendif /* (addSIPplot > 0) */

    nr_inst_vec=np.zeros(nr_sims,dtype='int')
    Vtot_vec=np.zeros(nr_sims)
    skal_m_vec=np.zeros(nr_sims)

    for i_sim in range(0,nr_sims):
        fp=filepath_input[i_sim]
        print('process simulation:',fp)
#Einlesen SIP-Daten
        nEK_sip_plot_sim, mEK_sip_plot_sim,zDummy, nr_SIPs_plot_sim, nr_SIPs_prof_sim, t_vec_GVplot, V, skal_m = IO.readSingleSimulationSIPdata(fp,nr_sip_max=nr_sip_max)
        nr_inst,nr_GVplot,nr_SIPs = nEK_sip_plot_sim.shape
        Vtot_vec[i_sim] = V
        skal_m_vec[i_sim] = skal_m

        print('nr_inst,nr_GVplot,nr_SIPs',nr_inst,nr_GVplot,nr_SIPs)
        print(nr_SIPs_plot_sim.shape)
        nr_inst_vec[i_sim] = nr_inst
        if (nr_inst > nr_inst_max): 
            print ('increase nr_inst_max: ', nr_inst, nr_inst_max )
        nr_SIPs_plot[i_sim,0:nr_inst,:nr_GVplot] = nr_SIPs_plot_sim
        nEK_sip_plot[i_sim,0:nr_inst,:nr_GVplot,:nr_SIPs]= nEK_sip_plot_sim
        mEK_sip_plot[i_sim,0:nr_inst,:nr_GVplot,:nr_SIPs]= mEK_sip_plot_sim

        #GCCif (addSIPplot > 0)
        if (igetSIPdataonly == 1):
            nr_SIPs_av[i_sim,:nr_GVplot]=np.mean(nr_SIPs_plot_sim,axis=0)
            prof_dim=nr_SIPs_prof_sim.shape # nr_inst,nr_GVplot, nz
            nz = prof_dim[2]
            #GCCif (addSIPplot == 1)
            nr_SIPs_prof_av[i_sim,:nr_GVplot,:nz]=np.mean(nr_SIPs_prof_sim,axis=0)*100/dz_vec[i_sim]
            #GCCendif /* (addSIPplot == 1) */
            #GCCif (addSIPplot == 2)
            nr_SIPs_prof_av[i_sim,:nr_GVplot,:nz]=np.mean(nr_SIPs_prof_sim,axis=0)
            #GCCendif /* (addSIPplot == 2) */
            ntGV_vec[i_sim] = nr_GVplot
            t_vec_GVplot_vec[i_sim,:nr_GVplot]= t_vec_GVplot
        #GCCendif /* (addSIPplot > 0) */

    if (igetSIPdataonly == 0):
    #Erstelle GV-Plot
        for Times_elem in Times_list: 
            #GCCif (MOM_skal_m == 0)
            PS.PlotGV(nEK_sip_plot,mEK_sip_plot,nr_SIPs_plot,t_vec_GVplot,
                    iMultipleSims=1,nr_inst_vec=nr_inst_vec,
                    label=label,title=title,ilastGVonly=ilastGVonly,iTimes_select=Times_elem,
                    V=Vtot_vec)
            #GCCelse
            PS.PlotGV(nEK_sip_plot,mEK_sip_plot,nr_SIPs_plot,t_vec_GVplot,
                    iMultipleSims=1,nr_inst_vec=nr_inst_vec,
                    label=label,title=title,ilastGVonly=ilastGVonly,iTimes_select=Times_elem,
                    V=Vtot_vec, skal_m=skal_m_vec)
            #GCCendif /* else (MOM_skal_m == 0) */
    else:
        return nr_SIPs_av, nr_SIPs_prof_av, t_vec_GVplot_vec, ntGV_vec


#<<<<<<<< end routine plotMultipleSimulationGV(filepath_input)

    #GCCif (DISCRETE >= 1)
def plotMultipleSimulationRStD(filepath_input,label=None,title=None):
    nr_sims=len(filepath_input)
    print('call of plotMultipleSimulationRStD')
    print('nr_sims:',nr_sims)

    mEK_sip_plot=np.zeros([nr_sims,nr_inst_max,nr_GVplot_max,nr_sip_max])
    nr_inst_vec=np.zeros(nr_sims,dtype='int')

    for i_sim in range(0,nr_sims):
        fp=filepath_input[i_sim]
        print('process simulation:',fp)
#Einlesen SIP-Daten
        nEK_sip_plot_sim, mEK_sip_plot_sim,zDummy, nr_SIPs_plot_sim, nr_SIPs_prof_sim, t_vec_GVplot, V, skal_m = IO.readSingleSimulationSIPdata(fp,nr_sip_max=nr_sip_max)
        nr_inst,nr_GVplot = nr_SIPs_plot_sim.shape
        nr_inst_vec[i_sim] = nr_inst
        #nr_SIPs_plot[i_sim,0:nr_inst,:] = nr_SIPs_plot_sim
        mEK_sip_plot[i_sim,0:nr_inst,:nr_GVplot,:]= mEK_sip_plot_sim

#Erstelle GV-Plot
    PS.PlotRelStdDev(mEK_sip_plot,t_vec_GVplot,iMultipleSims=1,nr_inst_vec=nr_inst_vec,label=label,title=title)
#<<<<<<<< end routine plotMultipleSimulationRStD(filepath_input,label=None,title=None)

    #GCCendif /* (DISCRETE >= 1) */

def plotMultipleSimulationFLUX(filepath_input):
    nr_sims=len(filepath_input)
    print('call of plotMultipleSimulationFLUX')
    print('nr_sims:',nr_sims)
    nt_max=80

    #GCCif (STD == 1)
    iplotStdDev = 1
    FluxIn_STD= np.zeros([1,nr_sims,nt_max,5])
    FluxOut_STD = np.zeros([1,nr_sims,nt_max,5])
    FluxInAcc_STD = np.zeros([1,nr_sims,nt_max,5])
    FluxOutAcc_STD = np.zeros([1,nr_sims,nt_max,5])
    #GCCendif /* (STD == 1) */
    #GCCif (STD == 2)
    iplotStdDev = 2
    FluxIn_STD= np.zeros([2,nr_sims,nt_max,5])
    FluxOut_STD = np.zeros([2,nr_sims,nt_max,5])
    FluxInAcc_STD = np.zeros([2,nr_sims,nt_max,5])
    FluxOutAcc_STD = np.zeros([2,nr_sims,nt_max,5])
    #GCCendif /* (STD == 2) */
    #GCCif (STD == 0)
    iplotStdDev = 0
    FluxIn_STD= None
    FluxOut_STD = None
    FluxInAcc_STD = None
    FluxOutAcc_STD = None
    #GCCendif /* (STD == 0) */

    FluxIn= np.zeros([nr_sims,nt_max,5])
    FluxOut = np.zeros([nr_sims,nt_max,5])
    FluxInAcc = np.zeros([nr_sims,nt_max,5])
    FluxOutAcc = np.zeros([nr_sims,nt_max,5])
    print('FluxIn.shape,FluxOut.shape,FluxInAcc.shape,FluxOutAcc.shape')
    print(FluxIn.shape,FluxOut.shape,FluxInAcc.shape,FluxOutAcc.shape)

    nt_vec      =np.zeros(nr_sims,dtype='int')
    t_vec_Flux_vec=np.zeros([nr_sims,nt_max])
    #read flux data
    for i in range(nr_sims):
        fp_in=filepath_input[i]
        print('process simulation:',fp_in)
        FluxIn_tmp,FluxOut_tmp,FluxInAcc_tmp,FluxOutAcc_tmp,t_vec_GVplot=IO.readSingleSimulationFluxes(fp_in)

        nt_vec[i] = t_vec_GVplot.size-1
        t_vec_Flux_vec[i,:nt_vec[i]] = t_vec_GVplot[1:]
        FluxIn[i,:nt_vec[i],:]=np.mean(FluxIn_tmp,axis=0)
        FluxOut[i,:nt_vec[i],:]=np.mean(FluxOut_tmp,axis=0)
        FluxInAcc[i,:nt_vec[i],:]=np.mean(FluxInAcc_tmp,axis=0)
        FluxOutAcc[i,:nt_vec[i],:]=np.mean(FluxOutAcc_tmp,axis=0)

        #GCCif (STD == 1)
        FluxIn_STD[0,i,:nt_vec[i],:]=np.std(FluxIn_tmp,axis=0)
        FluxOut_STD[0,i,:nt_vec[i],:]=np.std(FluxOut_tmp,axis=0)
        FluxInAcc_STD[0,i,:nt_vec[i],:]=np.std(FluxInAcc_tmp,axis=0)
        FluxOutAcc_STD[0,i,:nt_vec[i],:]=np.std(FluxOutAcc_tmp,axis=0)
        #GCCendif /* (STD == 1) */
        #GCCif (STD == 2)
        FluxIn_STD[:,i,:nt_vec[i],:]=np.abs(np.percentile(FluxIn_tmp,[10,90],axis=0)-np.squeeze(FluxIn[i,:nt_vec[i],:]))
        FluxOut_STD[:,i,:nt_vec[i],:]=np.abs(np.percentile(FluxOut_tmp,[10,90],axis=0)-np.squeeze(FluxOut[i,:nt_vec[i],:]))
        FluxInAcc_STD[:,i,:nt_vec[i],:]=np.abs(np.percentile(FluxInAcc_tmp,[10,90],axis=0)-np.squeeze(FluxInAcc[i,:nt_vec[i],:]))
        FluxOutAcc_STD[:,i,:nt_vec[i],:]=np.abs(np.percentile(FluxOutAcc_tmp,[10,90],axis=0)-np.squeeze(FluxOutAcc[i,:nt_vec[i],:]))
        #GCCendif /* (STD == 2) */

    #GCCif (STD > 0)
    print('FluxIn_STD.shape,FluxOut_STD.shape,FluxInAcc_STD.shape,FluxOutAcc_STD.shape')
    print(FluxIn_STD.shape,FluxOut_STD.shape,FluxInAcc_STD.shape,FluxOutAcc_STD.shape)
    #GCCendif /* (STD > 0) */

    #generate flux plot
    PS.PlotFluxesTime(FluxIn,FluxOut,FluxInAcc,FluxOutAcc,t_vec_Flux_vec,
                      iMean=1,iMultipleSims=1,io_sep=1,nt_vec=nt_vec,
                      label=label,title=title,text=text,
                      iplotStdDev = iplotStdDev, FluxIn_STD = FluxIn_STD, FluxOut_STD = FluxOut_STD, FluxInAcc_STD = FluxInAcc_STD, FluxOutAcc_STD = FluxOutAcc_STD)

def plotMultipleSimulationGVout(filepath_input,infalling=0):
    nr_sims=len(filepath_input)
    print('call of plotMultipleSimulationGVout')
    print('nr_sims:',nr_sims)

    nEK_sip_out=np.zeros([nr_sims,nr_inst_max,nr_sip_max])
    mEK_sip_out=np.zeros([nr_sims,nr_inst_max,nr_sip_max])
    tEK_sip_out=np.zeros([nr_sims,nr_inst_max,nr_sip_max])
    nr_SIPs_out=np.zeros([nr_sims,nr_inst_max],dtype='int')
    nr_inst_vec=np.zeros(nr_sims,dtype='int')
    #A_vec=np.zeros(nr_sims)

    for i_sim in range(0,nr_sims):
        fp=filepath_input[i_sim]
        print('process simulation:',fp)
#Einlesen SIP-Daten
        nEK_sip_out_sim, mEK_sip_out_sim, tEK_sip_out_sim, nr_SIPs_out_sim, Tsim, A = IO.readSingleSimulationSIP_outfalling_data(fp,infalling)

        nr_inst,nr_SIPs = nEK_sip_out_sim.shape
        #A_vec[i_sim] = A_sim
        print('nr_SIPs_out_sim', nr_SIPs_out_sim)
        print('nr_inst,nr_SIPs', nr_inst, nr_SIPs)
        print(nr_SIPs_out_sim.shape)
        nr_inst_vec[i_sim] = nr_inst
        if (nr_inst > nr_inst_max): 
            print ('increase nr_inst_max: ', nr_inst, nr_inst_max )
        nr_SIPs_out[i_sim,0:nr_inst] = nr_SIPs_out_sim
        nEK_sip_out[i_sim,0:nr_inst,:nr_SIPs]= nEK_sip_out_sim/A
        mEK_sip_out[i_sim,0:nr_inst,:nr_SIPs]= mEK_sip_out_sim

#Erstelle GV-Plot
    PS.PlotGV(nEK_sip_out,mEK_sip_out,nr_SIPs_out,[3600],iMultipleSims=1,nr_inst_vec=nr_inst_vec,label=label,title=title, outfallingGV=1+infalling)
#<<<<<<<< end routine plotMultipleSimulationGVout(filepath_input)

#GCCendif /* (IMODE == 2) */


#--------------- end defintion of subroutines --------------------------------------------------------------------


#GCCif (IMODE == 1)
print('filepath',filepath)

    #GCCif (GV_meanTD == 1)
plotSingleSimulationGV(filepath)
    #GCCendif /*(GV_meanTD == 1)*/

    #GCCif (rSIGMA_time  >= 1)
plotSingleSimulationRStD(filepath)
    #GCCendif /* (rSIGMA_time >= 1) */

    #GCCif (MOM_meanTD >= 1)
plotSingleSimulationMOM(filepath)
    #GCCendif /* (MOM_meanTD >= 1) */

    #GCCif (COLUMN == 1)

        #GCCif (FLUX_time == 1)
plotSingleSimulationFLUX(filepath)
        #GCCendif /* (FLUX_time == 1) */

        #GCCif (MOM_prof == 1)
plotSingleSimulationMOM(filepath,iprof=1,itot=0)
        #GCCendif /* (MOM_prof == 1) */

        #GCCif (GV_out == 1)
plotSingleSimulationGVout(filepath)
        #GCCendif /*(GV_out == 1)*/

        #GCCif (GV_in == 1)
plotSingleSimulationGVout(filepath, infalling=1)
        #GCCendif /*(GV_in == 1)*/
    #GCCendif /* (COLUMN == 1) */

#GCCendif /* (IMODE == 1) */

#GCCif (IMODE == 2)

    #GCCif (GV_meanTD == 1)
plotMultipleSimulationGV(filepath_input)
    #GCCendif /*(GV_meanTD == 1)*/

    #GCCif (rSIGMA_time >= 1)
plotMultipleSimulationRStD(filepath_input)
    #GCCendif /* (rSIGMA_time >= 1) */

    #GCCif (MOM_meanTD >= 1 || MOM_prof == 1)
iProf=None
iTime=None
    #GCCif (MOM_meanTD >= 1)
iTime=1
    #GCCendif /* (MOM_meanTD >= 1) */
    #GCCif (MOM_prof == 1)
iProf=1
    #GCCendif /* (MOM_prof == 1) */
plotMultipleSimulationsMOM(filepath_input,iTime=iTime,iProf=iProf)
    #GCCendif /* (MOM_meanTD >= 1 || MOM_prof == 1) */

    #GCCif (FLUX_time == 1)  
plotMultipleSimulationFLUX(filepath_input)
    #GCCendif /* (FLUX_time == 1) */

    #GCCif (GV_out == 1)
plotMultipleSimulationGVout(filepath_input)
    #GCCendif /*(GV_out == 1)*/

    #GCCif (GV_in == 1)
plotMultipleSimulationGVout(filepath_input, infalling=1)
    #GCCendif /*(GV_in == 1)*/

#GCCendif /* (IMODE == 2 ) */
#print(os.path.dirname(os.path.abspath(__file__)))

