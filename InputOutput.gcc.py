import os
import numpy as np
import Misc as FK
#GCCif (COLUMN == 0)
def readSingleSimulationSIPdata(fp_in, nr_sip_max=500):
# >>>> box model version >>>>>>
# Read SIP data of a single simulation.
# In the default case, all SIP data are stored in the specified folder fp_in.
# If a simulation was split into several subs-simulations, then read SIP data from the various subfolders.

# default case, where all SIP data lie in folder fp_in
    i_singlefolder = 1
    dirs_subsim=['']
    # nr_inst_tot not yet known
    nr_subsims = 1

# Cases with sub-simulations
    # if a given simulation is divided into sub-simulations, then by definition file InstanzenMeta.dat exists

    if os.path.isfile(fp_in + 'InstanzenMeta.dat'):
        i_singlefolder = 0
        fu = open(fp_in + 'InstanzenMeta.dat','r')
        nr_subsims = int(fu.readline())
        nr_inst_tot = int(fu.readline())
        dirs_subsim = fu.readlines()
        for i,s in enumerate(dirs_subsim):  # remove trailing \n and other whitespace
            dirs_subsim[i] = s.strip()+'/'

        print('nr_subsims, nr_inst_tot: ', nr_subsims, nr_inst_tot)
        print(dirs_subsim)
        fu.close()

    itot_inst = 0
    for i_subsim,folder in enumerate(dirs_subsim):
    # read SIP data of all realisations
        fp = fp_in + folder
        fu = open(fp + 'SIP_meta.dat','r')
        nr_inst = int(fu.readline())
        nr_GVplot = int(fu.readline())
        t_vec_GVplot= np.array(fu.readline().split(), dtype='float' )
        #print(t_vec_GVplot)
        fu.close()
        fu = open(fp + 'Moments_meta.dat','r')
        dV = float(fu.readline())
        skal_m = float(fu.readline())
        fu.close()

        if (i_subsim == 0):
            if (i_singlefolder == 1): nr_inst_tot=nr_inst
            nEK_sip_plot=np.zeros([nr_inst_tot,nr_GVplot,nr_sip_max])
            mEK_sip_plot=np.zeros([nr_inst_tot,nr_GVplot,nr_sip_max])
            nr_SIPs_plot=np.zeros([nr_inst_tot,nr_GVplot],dtype='int')

        fu = FK.openfile(fp + 'SIP.dat',options='r')
        for i_inst in range(0,nr_inst):
            for i_plot in range(0,nr_GVplot):
                nr_SIPs = int(fu.readline())
                nr_SIPs_plot[itot_inst,i_plot] = nr_SIPs
                #print(i_inst,iplot,nr_SIPs)
                strline=fu.readline().split()
                nEK_sip_plot[itot_inst,i_plot,0:nr_SIPs]= np.array( strline, dtype='float' ) #, sep=','
                strline=fu.readline().split()
                mEK_sip_plot[itot_inst,i_plot,0:nr_SIPs]= np.array( strline, dtype='float' )
            itot_inst += 1
        fu.close()

    return nEK_sip_plot, mEK_sip_plot, None, nr_SIPs_plot, None, t_vec_GVplot, dV, skal_m
#<<<<<<<<< box model version <<<<<<<<<<<<<
#<<<<<<<< end subroutine readSingleSimulationSIPdata

def readSingleSimulationMoments(fp_in):
# >>>> box model version >>>>>>
# box model version
# read moments data of a single simulation
# In the default case, all moment data are stored in the specified folder fp_in.
# If a simulation was split into several subs-simulations, then read moment data from the various subfolders.

# default case, where all SIP data lie in folder fp_in
    i_singlefolder = 1
    dirs_subsim=['']
    # nr_inst_tot not yet known
    nr_subsims = 1

# Cases with sub-simulations
    # if a given simulation is divided into sub-simulations, then by definition file InstanzenMeta.dat exists

    if os.path.isfile(fp_in + 'InstanzenMeta.dat'):
        i_singlefolder = 0
        fu = open(fp_in + 'InstanzenMeta.dat','r')
        nr_subsims = int(fu.readline())
        nr_inst_tot = int(fu.readline())
        dirs_subsim = fu.readlines()
        for i,s in enumerate(dirs_subsim):  # remove trailing \n and other whitespace
            dirs_subsim[i] = s.strip()+'/'

        print('nr_subsims, nr_inst_tot: ', nr_subsims, nr_inst_tot)
        print(dirs_subsim)
        fu.close()

    itot_inst = 0
    for i_subsim,folder in enumerate(dirs_subsim):
        fp = fp_in + folder

        # read moments data
        fu = open(fp + 'Moments_meta.dat','r')
        dV = float(fu.readline())
        skal_m = float(fu.readline())
        nr_inst = int(fu.readline())
        nr_MOMsave = int(fu.readline())
        t_vec_MOMsave= np.array(fu.readline().split(), dtype='float' )
        #print(t_vec_GVplot)
        fu.close()

        if (i_subsim == 0):
            if (i_singlefolder == 1): nr_inst_tot=nr_inst
            MOMsave=np.zeros([nr_inst_tot,nr_MOMsave,4])
        fu = FK.openfile(fp + 'Moments.dat',options='rb')
        MOMsave[itot_inst:itot_inst+nr_inst,:,:]=np.loadtxt(fu).reshape((nr_inst,nr_MOMsave,4))
        itot_inst += nr_inst

    for i in range (4): MOMsave[:,:,i] = MOMsave[:,:,i]*skal_m**i/dV

    return MOMsave, t_vec_MOMsave, skal_m, dV
#<<<<<<<<< box model version <<<<<<<<<<<<<
#<<<<<<<< end subroutine readSingleSimulationMoments
#GCCendif /* (COLUMN == 0) */
#<<<<<<<<<<<<<<<<<<<<<<<<< box model <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


#>>>>>>>>>>>>>>>>>>>>>>>>> column model >>>>>>>>>>>>>>>>>>>>>>>>>>>>
#GCCif (COLUMN == 1)

def get_SubSimsData(fp_in):
    fu = open(fp_in + 'InstanzenMeta.dat','r')
    nr_subsims = int(fu.readline())
    nr_inst_tot = int(fu.readline())
    dirs_subsim = fu.readlines()
    for i,s in enumerate(dirs_subsim):  # remove trailing \n and other whitespace
        dirs_subsim[i] = s.strip()+'/'

    print('nr_subsims, nr_inst_tot: ', nr_subsims, nr_inst_tot)
    print(dirs_subsim)
    fu.close()
    return 0,nr_subsims,nr_inst_tot,dirs_subsim

def get_nz(fp):
    return get_MetaData(fp)[0]

def get_dz(fp):
    return get_MetaData(fp)[1]

def get_MetaData(fp):
    #get meta data of single column model simulation
    dirs_subsim=['']
    #check existence of file InstanzenMeta.dat. If it exists, then read Meta.dat of first sub-simulation.
    if os.path.isfile(fp + 'InstanzenMeta.dat'):
        i_singlefolder,nr_subsims,nr_inst_tot,dirs_subsim = get_SubSimsData(fp)

    fMeta = open(fp +dirs_subsim[0] +'Meta.dat','r')
    str_array = fMeta.readline().split()
    nz=int(str_array[0])
    dz=float(str_array[1])
    Tsim=float(str_array[2])
    nr_inst=int(str_array[3])
    str_array = fMeta.readline().split()
    LWC=float(str_array[0])
    r0=float(str_array[1])
    xf0=float(str_array[2])
    ikernel=int(str_array[3])
    i_init_1D=int(str_array[4])
    i_process=int(str_array[5])

    return nz,dz,Tsim,nr_inst,LWC,r0,xf0,ikernel,i_init_1D,i_process
#<<<<<<<< end routine get_MetaData

def readSingleSimulationSIPdata(fp_in, nr_sip_max=1):
# >>>>colum model version >>>>>>
# Read SIP data of a single simulation.
# In the default case, all SIP data are stored in the specified folder fp_in.
# If a simulation was split into several subs-simulations, then read SIP data from the various subfolders.# nr_sip_max is dummy argument, introduced for consistency with COLUMN =0 version of this routine

# default case, where all SIP data lie in folder fp_in
    i_singlefolder = 1
    dirs_subsim=['']
    # nr_inst_tot not yet known
    nr_subsims = 1

# Cases with sub-simulations
    # if a given simulation is divided into sub-simulations, then by definition file InstanzenMeta.dat exists

    if os.path.isfile(fp_in + 'InstanzenMeta.dat'):
        i_singlefolder,nr_subsims,nr_inst_tot,dirs_subsim = get_SubSimsData(fp_in)

    nr_inst=np.zeros(nr_subsims,dtype='int')

    for i_subsim,folder in enumerate(dirs_subsim):
    # Read SIP meta data of all realisations
    # get maximum number of SIPs of full column and of a single grid box over all realisations.
        fp = fp_in + folder
        fGVMeta = open(fp + 'SIP_meta.dat','r')
        nr_inst[i_subsim] = int(fGVMeta.readline())
        nr_GVplot = int(fGVMeta.readline())
        t_vec_GVplot= np.array(fGVMeta.readline().split(), dtype='float' )
        if (i_subsim == 0):
            if (i_singlefolder == 1): nr_inst_tot=nr_inst[i_subsim]
            nr_SIPs_max=0 #total SIP number over full column, maximum over time and instances
            nr_SIPs_GB_maxmax=0 # nr_SIPs_GB_max: highest SIP number in a single grid box; nr_SIPs_GB_maxmax: take maximum over time and instances
        for i_inst in range(nr_inst[i_subsim]):
            fGVMeta.readline()
            for i_GVplot in range(nr_GVplot):
                [nr_SIPs,nr_SIPs_GB_max]=np.array(fGVMeta.readline().split(), dtype='int')
                nr_SIPs_max=max([nr_SIPs,nr_SIPs_max])
                nr_SIPs_GB_maxmax=max([nr_SIPs_GB_max,nr_SIPs_GB_maxmax]) # currently not used in the following evaluation
        fGVMeta.close()

    itot_inst = 0
    for i_subsim,folder in enumerate(dirs_subsim):
    # read SIP data of all realisations
        if (i_subsim == 0):
            if (i_singlefolder == 1): nr_inst_tot=nr_inst[0]
            nEK_sip_plot=np.zeros([nr_inst_tot,nr_GVplot,nr_SIPs_max])
            mEK_sip_plot=np.zeros([nr_inst_tot,nr_GVplot,nr_SIPs_max])
            zEK_sip_plot=np.zeros([nr_inst_tot,nr_GVplot,nr_SIPs_max])
            nz=get_nz(fp)
            nr_SIPs_prof=np.zeros([nr_inst_tot,nr_GVplot,nz],dtype='int')
            nr_SIPs_plot=np.zeros([nr_inst_tot,nr_GVplot],dtype='int')
            fu = open(fp + 'Moments_meta.dat','r')
            dV = float(fu.readline())
            skal_m = float(fu.readline())
            fu.close()
            V = nz * dV
        fGV = FK.openfile(fp + 'SIP.dat',options='r')
        print(folder)
        for i_inst in range(0,nr_inst[i_subsim]):
            #print(i_subsim,i_inst,itot_inst)
            for i_plot in range(0,nr_GVplot):
                iSIP=0

                for k in range(nz):
                    [nr_SIPs,iz]=np.array(fGV.readline().split(),dtype='int')
                    nr_SIPs_prof[itot_inst,i_plot,k]= nr_SIPs
                    if (nr_SIPs > 0):
                        ia=iSIP
                        ie=ia+nr_SIPs
                        nEK_sip_plot[itot_inst,i_plot,ia:ie]=np.array(fGV.readline().split())
                        mEK_sip_plot[itot_inst,i_plot,ia:ie]=np.array(fGV.readline().split())
                        zEK_sip_plot[itot_inst,i_plot,ia:ie]=np.array(fGV.readline().split())
                        iSIP=iSIP+nr_SIPs
                nr_SIPs_plot[itot_inst,i_plot]= iSIP
            itot_inst += 1
        fGV.close()

    return nEK_sip_plot, mEK_sip_plot, zEK_sip_plot, nr_SIPs_plot, nr_SIPs_prof, t_vec_GVplot, V, skal_m
#<<<<<<<<<colum model version <<<<<<<<<<<<<
#<<<<<<< end subroutine readSingleSimulationSIPdata

def readSingleSimulationMoments(fp_in):
# >>>>colum model version >>>>>>
# read Moments.dat in folder fp_in
# if simulation is divided in several sub-simulations, then read according files from all sub-folders

# Default Case, where simulation data lie in specified folder fp_in
    i_singlefolder = 1
    dirs_subsim=['']
    # nr_inst_tot not yet known
    nr_subsims = 1

#Case, where multiple sub-simulations exist and data lie in sub-folders
    #if file InstanzenMeta.dat exists, then sub-simulations exist (the sub-folders are listed in InstanzenMeta.dat)

    if os.path.isfile(fp_in + 'InstanzenMeta.dat'):
        i_singlefolder,nr_subsims,nr_inst_tot,dirs_subsim = get_SubSimsData(fp_in)

    itot_inst = 0
    for i_subsim,folder in enumerate(dirs_subsim):
        fp = fp_in + folder

        #Read meta file
        fu = open(fp + 'Moments_meta.dat','r')
        dV = float(fu.readline())
        skal_m = float(fu.readline())
        #print(dV,skal_m)
        nr_inst = int(fu.readline())
        print('i_subsim, nr_inst', i_subsim, nr_inst,itot_inst)
        nr_MOMsave = int(fu.readline())
        t_vec_MOMsave= np.array(fu.readline().split(), dtype='float' )
        #print(t_vec_GVplot)
        fu.close()

        if (i_subsim == 0):
            if (i_singlefolder == 1): nr_inst_tot=nr_inst
            nz=get_nz(fp)
            MOMsave=np.zeros([nr_inst_tot,nr_MOMsave,nz,4])
            print('MOMsave.shape', MOMsave.shape)

        fu = FK.openfile(fp + 'Moments.dat',options='rb')
        MOMsave[itot_inst:itot_inst+nr_inst,:,:,:]=np.loadtxt(fu).reshape((nr_inst,nr_MOMsave,nz,4))
        itot_inst += nr_inst
    for i in range (4): MOMsave[:,:,:,i] = MOMsave[:,:,:,i]*(skal_m**i)/dV
    return MOMsave, t_vec_MOMsave, skal_m, dV
#<<<<<<<<<colum model version <<<<<<<<<<<<<
#<<<<<<<< end subroutine readSingleSimulationMoments

def readSingleSimulationMomentsMeanFinal(fp_in, iPlotextras=0, iTimesCross=0, thresh=None):
# >>>>colum model version >>>>>>
# read MomentsMean.dat in folder fp_in
# if simulation is divided in several sub-simulations then read according files from all sub-folders

# Default Case, where simulation data lie in specified folder fp_in
    i_singlefolder = 1
    dirs_subsim=['']
    # nr_inst_tot not yet known
    nr_subsims = 1

#Case, where multiple sub-simulations exist and data lie in sub-folders
    #if file InstanzenMeta.dat exists, then sub-simulations exist (the sub-folders are listed in InstanzenMeta.dat)

    if os.path.isfile(fp_in + 'InstanzenMeta.dat'):
        i_singlefolder,nr_subsims,nr_inst_tot,dirs_subsim = get_SubSimsData(fp_in)

    for i_subsim,folder in enumerate(dirs_subsim):
        fp = fp_in + folder

        #Read meta file
        fu = open(fp + 'Moments_meta.dat','r')
        dV = float(fu.readline())
        skal_m = float(fu.readline())
        #print(dV,skal_m)
        nr_inst = int(fu.readline())
        nr_MOMsave = int(fu.readline())
        t_vec_MOMsave= np.array(fu.readline().split(), dtype='float' )
        #print(t_vec_GVplot)
        fu.close()

        if (i_subsim == 0):
            if (i_singlefolder == 1): nr_inst_tot=nr_inst
            nz=get_nz(fp)
            if (iTimesCross == 0):
                MOMsave=np.zeros([4])
            else:
                MOMsave = np.zeros([nr_MOMsave,4])
                indices = np.zeros([2,3],dtype='int')

        fu = FK.openfile(fp + 'MomentsMean.dat',options='rb')
        if (iTimesCross == 0):
            if (fu is None):
                return 0, 0, 0, 0
            nt_time_pick = nr_MOMsave-1
            if (iPlotextras == 10):
                nt_time_pick = 72
            MOMsave[:] += nr_inst*np.loadtxt(fu).reshape((nr_MOMsave,4))[nt_time_pick,:]
        else:
            if (fu is None):
                return 0, 0, 0, 0
            MOMsave[:,:] += nr_inst*np.loadtxt(fu).reshape((nr_MOMsave,4))

    # Averaging over final data
    MOMsave = MOMsave/dV/nr_inst_tot

    if (iTimesCross == 0):
        for i in range (4):
            MOMsave[i] = MOMsave[i]*(skal_m**i)
        return MOMsave, skal_m, dV, 1
    else:
        iMom = 0
        iiMom = 0
        for ithr in range(3):
            print('iMom,ithr',iMom,ithr)
            indices[iiMom,ithr] = np.where(MOMsave[:,iMom].reshape((nr_MOMsave)) < thresh[iiMom,ithr])[0][0]
        iMom = 2
        iiMom = 1
        for ithr in range(3):
            print('iMom,ithr',iMom,ithr)
            indices[iiMom,ithr] = np.where(MOMsave[:,iMom].reshape((nr_MOMsave)) > thresh[iiMom,ithr])[0][0]

        return t_vec_MOMsave[indices]


#<<<<<<<<<colum model version <<<<<<<<<<<<<
#<<<<<<<< end subroutine readSingleSimulationMomentsMean


def readSingleSimulationFluxes(fp_in):
# >>>>colum model version >>>>>>
# read flux data of a single simulation
# if simulation is divided in several sub-simulations then read according files from all sub-folders

# Default Case, where simulation data lie in specified folder fp_in
    i_singlefolder = 1
    dirs_subsim=['']
    # nr_inst_tot not yet known
    nr_subsims = 1

#Case, where multiple sub-simulations exist and data lie in sub-folders
    #if file InstanzenMeta.dat exists, then sub-simulations exist (the sub-folders are listed in InstanzenMeta.dat)

    if os.path.isfile(fp_in + 'InstanzenMeta.dat'):
        i_singlefolder,nr_subsims,nr_inst_tot,dirs_subsim = get_SubSimsData(fp_in)

    #nz,dz,Tsim,nr_inst,LWC,r0,xf0,ikernel,i_init_1D,i_process=get_MetaData(fp_in)

    itot_inst = 0
    for i_subsim,folder in enumerate(dirs_subsim):
        fp = fp_in + folder

        fGVMeta = open(fp + 'SIP_meta.dat','r')
        nr_inst = int(fGVMeta.readline())
        nr_GVplot = int(fGVMeta.readline())
        t_vec_GVplot= np.array(fGVMeta.readline().split(), dtype='float' )
        fGVMeta.close()

        if (i_subsim == 0):
            if (i_singlefolder == 1): nr_inst_tot=nr_inst
            FluxIn=np.zeros([nr_inst_tot,nr_GVplot-1,5])
            FluxOut=np.zeros([nr_inst_tot,nr_GVplot-1,5])
            FluxInAcc=np.zeros([nr_inst_tot,nr_GVplot-1,5])
            FluxOutAcc=np.zeros([nr_inst_tot,nr_GVplot-1,5])

        #FluxIn,FluxOut,FluxInAcc,FluxOutAcc
        FluxIn[itot_inst:itot_inst+nr_inst,:,:]=np.loadtxt(fp+'Fluxes_in.dat').reshape(nr_inst,nr_GVplot-1,5)
        FluxInAcc[itot_inst:itot_inst+nr_inst,:,:]=np.loadtxt(fp+'Fluxes_in_acc.dat').reshape(nr_inst,nr_GVplot-1,5)
        FluxOut[itot_inst:itot_inst+nr_inst,:,:]=np.loadtxt(fp+'Fluxes_out.dat').reshape(nr_inst,nr_GVplot-1,5)
        FluxOutAcc[itot_inst:itot_inst+nr_inst,:,:]=np.loadtxt(fp+'Fluxes_out_acc.dat').reshape(nr_inst,nr_GVplot-1,5)
        itot_inst += nr_inst

    return FluxIn,FluxOut,FluxInAcc,FluxOutAcc,t_vec_GVplot
#<<<<<<<<<colum model version <<<<<<<<<<<<<
#<<<<<<<< end subroutine readSingleSimulationFluxes


def readSingleSimulationSIP_outfalling_data(fp_in,infalling=0):
# if simulation is divided in several sub-simulations then read according files from all sub-folders
    #fu_log = open(fp_in + 'LOG_out.txt','w')
# Default Case, where simulation data lie in specified folder fp_in
    i_singlefolder = 1
    dirs_subsim = ['']
    # nr_inst_tot not yet known
    nr_subsims = 1
    suffix=['out','in'][infalling]

#Case, where multiple sub-simulations exist and data lie in sub-folders
    #if file InstanzenMeta.dat exists, then sub-simulations exist (the sub-folders are listed in InstanzenMeta.dat)

    if os.path.isfile(fp_in + 'InstanzenMeta.dat'):
        i_singlefolder, nr_subsims, nr_inst_tot, dirs_subsim = get_SubSimsData(fp_in)

    nr_inst = np.zeros(nr_subsims, dtype='int')

    #get number of instances of each subsim
    for i_subsim, folder in enumerate(dirs_subsim):
        fp = fp_in + folder
        fGVMeta = open(fp + 'SIP_meta.dat','r')
        nr_inst[i_subsim] = int(fGVMeta.readline())
        nr_GVplot = int(fGVMeta.readline())
        t_vec_GVplot = np.array(fGVMeta.readline().split(), dtype='float' )
        fGVMeta.close()
        if (i_subsim == 0):
            if (i_singlefolder == 1): nr_inst_tot = nr_inst[i_subsim]

    print('nr_inst: ', nr_inst)
    # get number of outfalling SIPs per instance
    nr_SIPs_out = np.zeros(nr_inst_tot, dtype='int')
    nr_times_out = np.zeros(nr_inst_tot, dtype='int')
    itot_inst = 0
    for i_subsim, folder in enumerate(dirs_subsim):
        fp = fp_in + folder
        #print(fp)
        #tmpread = np.loadtxt(fp+'Fluxes_out_acc.dat').reshape(nr_inst[i_subsim], nr_GVplot-1, 5)
        #print(itot_inst, nr_inst[i_subsim], tmpread.shape, tmpread[:, nr_GVplot-2, 4].shape)
        #print('II', itot_inst,itot_inst+nr_inst[i_subsim])
        #print(nr_SIPs_out[itot_inst:itot_inst+nr_inst[i_subsim]].shape)
        nr_SIPs_out[itot_inst:itot_inst+nr_inst[i_subsim]] = \
           np.loadtxt(fp+'Fluxes_'+suffix+'_acc.dat').reshape(nr_inst[i_subsim], nr_GVplot-1, 5)[:, nr_GVplot-2, 4]
        #print('nrSIPs_out: ', nr_SIPs_out[itot_inst:itot_inst+nr_inst[i_subsim]])
        nr_times_out[itot_inst:itot_inst+nr_inst[i_subsim]] = \
            np.loadtxt(fp+'SIP'+suffix+'_meta.dat')
        #print('nrTimes_out: ', nr_times_out[itot_inst:itot_inst+nr_inst[i_subsim]])
        itot_inst += nr_inst[i_subsim]
    # end meta data processing

    nr_SIPs_out_max = nr_SIPs_out.max()
    itot_inst = 0
    for i_subsim, folder in enumerate(dirs_subsim):
        fp = fp_in + folder
        #Einlesen SIP-Daten aller Instanzen
        if (i_subsim == 0):
            nEK_sip_out = np.zeros([nr_inst_tot, nr_SIPs_out_max])
            mEK_sip_out = np.zeros([nr_inst_tot, nr_SIPs_out_max])
            tEK_sip_out = np.zeros([nr_inst_tot, nr_SIPs_out_max])
            dz=get_dz(fp)
            fu = open(fp + 'Moments_meta.dat','r')
            dV = float(fu.readline())
            fu.close()
            A = dV / dz
        fGVout = FK.openfile(fp + 'SIP'+suffix+'.dat',options='r')

        for i_inst in range(nr_inst[i_subsim]):
            iSIP = 0
            itimes_tot=0
            #print('BB', i_subsim, i_inst, nr_times_out[itot_inst])
            for itimes in range(nr_times_out[itot_inst]):
                [time, nr_SIPs]=np.array(fGVout.readline().split(), dtype='int')
                #fu_log.write("{:3d} {:3d} {:3d} {:5d} {:5d} {:5d} {:7d} {:7d} {:7d}\n".format(i_subsim, i_inst, itot_inst, itimes,itimes_tot, nr_times_out[itot_inst], time, iSIP, nr_SIPs))
                ia = iSIP
                ie = ia + nr_SIPs
                nEK_sip_out[itot_inst, ia:ie] = np.array(fGVout.readline().split())
                mEK_sip_out[itot_inst, ia:ie] = np.array(fGVout.readline().split())
                tEK_sip_out[itot_inst, ia:ie] = time
                iSIP += nr_SIPs
                itimes_tot += 1
            itot_inst += 1

        fGVout.close()

    return nEK_sip_out, mEK_sip_out, tEK_sip_out, nr_SIPs_out, t_vec_GVplot[nr_GVplot-1], A
#<<<<<<<<<colum model version <<<<<<<<<<<<<
#<<<<<<< end subroutine readSingleSimulationSIP_outfalling_data


#GCCendif /* (COLUMN == 1) */