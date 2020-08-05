import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import sys
import numpy as np
import math
import Misc as FK
import Referenzloesung as REF
from cycler import cycler
import string

print('matplotlib ', matplotlib.__version__)
import platform
print(platform.python_version())
iExtraLegend = 0
leg_extra = None
xticks = None
xlabel_Mom= [' ', 't / s', 't / min']
lsREF =':'
leg_ncol=1
iMom_select_const = [0,2] # None
labelMoments_dict={}
iPanelAnnotate=None
i_inline=0
iBin_noLegend = 0
nrPanels=1
iPanel=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
iPanelREF=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
iEmptydom = 0
iHalfDom  = 0
iBoxM_LWC_N_r_var = 0
iPlotextras = 0
iClip = 0
isavedata_filename = None
#GCCif (IREF == 9)
nrsimREF = 1
ntREFmax = 10000
#GCCendif /* (IREF == 9) */

#GCCif (DISCRETE == 0)
ytitle = '$g_{\mathrm{ln} r}$ / (g / cm$^3$)'
xscale="log"
    #GCCif (KERNEL == 1 || KERNEL == 2)
xlimGV=(1e0,5e3)
ylimGV=(1e-4,1e1)
    #GCCendif /* (KERNEL == 1 || KERNEL == 2) */
    #GCCif (KERNEL == 0)
xlimGV=(1e0,2e4)
ylimGV=(1e-4,1e1)
    #GCCendif /* (KERNEL == 0) */
#GCCendif /* (DISCRETE == 0) */

#GCCif (DISCRETE >= 1)
#GCCif (IREF == 2)
ytitle = 'droplet mass conc. / (g / cm$^3$)'
xscale="linear"
xlimGV=(15,60)
xticks=np.arange(15, 60, 5)

#xlimGV=(10,60)
#xticks=np.arange(10, 60, 5)
#xlimGV=(10,130)
#xticks=np.arange(10, 130, 10)
ylimGV=(1e-9,1e-6)
#GCCendif /* (IREF == 2) */
#GCCif (IREF == 3)
ytitle = 'droplet mass conc. / (g / cm$^3$)'
xscale="linear"
xlimGV=(10,70)
xticks=np.arange(10, 71, 10)
ylimGV=(1e-9,1e-6)
#GCCendif /* (IREF == 3) */
#GCCendif /* (DISCRETE >= 1) */

#GCCifdef COMP
#the following statement includes a parameter file via a preprocessor directive
#GCCinclude "params.txt"
#GCCdefine ICOMPPLOT 1
Time_cycler='ColorReg'
Sim_cycler='Line'
#GCCendif

#GCCifdef PLOT
#the following statement includes a parameter file via a preprocessor directive
#GCCinclude "params_plot.txt"
#GCCdefine ICOMPPLOT 2
#GCCendif

plt.rc('lines',markersize=4)

def PlotMoments(MOMsave,t_vec_MOMsave,fp_out='',
                iMean=0,iplot_onlyMOMmean=0,iplot_mode=0,
                label=[None],title=None,nr_sims=None,text=None,
                iDmean=0, iplotStdDev=0, MOM_StdDev=None, skal_m=None,
                iSIPplot=0, nr_SIPs_mean=None, t_vec_GVplot=None,
                iexcludetop=0
                ):
    print('iMom_select:',iMom_select_const,'end')
    if iMom_select_const is None:
        iMom_select = [0,1,2,3]
        nrMom_select = 4
    else:
        iMom_select = iMom_select_const.copy()
        nrMom_select= len(iMom_select)

    #generates plots with the temporal evolution of selected moments
    #MOMsave contains moment data, separately for each realisation (iMean=0)
    #MOMsave contains moment data, average over all realisations (iMean=1)
    #fp_out: output path
    #iplot_onlyMOMmean: 0 plot each realisation separately; =1 plot only ensemble average
    # the PlotMoments routine is called multiple times
    # if the routine is called with iplot_mode = 0, then the plot is initialised
    # if the routine is called with iplot_mode > 1 then curves are added one by one (the first call is with iplot_mode = 1, the last call with iplot_mode = 3, all intermediate calls with 2)

#GCCif (IREF == 9)
    if (iplot_mode == 0) or (iplot_mode == 3):
        ntREF_vec=np.zeros(nrsimREF, dtype='int')
        print(nrsimREF)
        if (nrsimREF > 1):
            MOMsaveREF = np.zeros([nrsimREF,4,ntREFmax])
            DmeanREF   = np.zeros([nrsimREF,ntREFmax])
        for isimREF in range(nrsimREF):
            ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF,rgridREF,mgridREF=REF.get_RefMetaData(isimREF=isimREF)
            print('REF meta data', ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF)
            ntREF_vec[isimREF] = ntREF
            if (ntREF > ntREFmax):
                print('increase ntREFmax', ntREF, ntREFmax)

            t_vecREF = TsimREF/(ntREF-1)*np.arange(ntREF)
            if (itime_legend == 2): t_vecREF = t_vecREF/60
            #print(t_vecREF)
            MOMtmp = REF.get_RefProfileData(nzREF,ntREF,isimREF=isimREF)  #function output : np.zeros([4,nt,nz])
            #MOMtmp = np.mean(MOMtmp, axis=2)
            if (iexcludetop == 0):
                MOMtmp = np.mean(MOMtmp, axis=2)
            else:
                nz_smallDom=int(nzREF/40*iexcludetop)
                print('nz_smallDom, nz: ', nz_smallDom, nzREF)
                MOMtmp = np.mean(MOMtmp[:,:,:nz_smallDom], axis=2)

            #GCCif (MOM_meanTD == 2)
            MOMtmp = MOMtmp * dzREF * nzREF
            #GCCendif /* (MOM_meanTD == 2) */
            #GCCif (MOM_meanTD == 3)
            MOMtmp = MOMtmp * dVREF * nzREF   # in bin model dV is not a relevant parameter and hence not specified. Specify dVREF in param section
            #GCCendif /* (MOM_meanTD == 3) */
            if (nrsimREF == 1):
                MOMsaveREF = np.expand_dims(MOMtmp, axis=0)
            else:
                MOMsaveREF[isimREF,:,:ntREF] = MOMtmp
            if (iDmean == 1):
                a = MOMtmp[1,:]
                b = MOMtmp[0,:]
                c = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
                Dmeantmp = FK.m2r(c,const_mass2rad)*2e6 # diameter in um
                if (nrsimREF == 1):
                    DmeanREF = np.expand_dims(Dmeantmp, axis=0)
                else:
                    DmeanREF[isimREF,:ntREF] = Dmeantmp
#GCCendif /* (IREF == 9) */

    global i_simMOM, ax_vec, fig0

    if (iMean == 1):
        MOMmean=MOMsave.copy()
        iplot_onlyMOMmean = 1
        MOMsave = np.expand_dims(MOMsave, axis=0) # add a dummy first dimension with length 1

    [nr_inst,nt_time,nr_Mom]=MOMsave.shape

    if (iMean == 0):
        #MOMmean=FK.CIO_MOMmean(data_in=MOMsave,skal_m=skal_m,dV=dV)
        MOMmean=FK.CIO_MOMmean(data_in=MOMsave,fp_out='')

    if (skal_m is not None):
        for i in range(4):
            MOMmean[:,i] = MOMmean[:,i] / skal_m**i
            MOMsave[:,:,i] = MOMsave[:,:,i] / skal_m**i
            if (MOM_StdDev is not None):
                print(MOM_StdDev.shape)
                MOM_StdDev[:,:,i] = MOM_StdDev[:,:,i] / skal_m**i

    if (iplotStdDev > 0 and MOM_StdDev is None):
        if (iMean == 1):
            print("if only Mean Values are provided, then StdDev/Percentiles can not be computed and must be provided as well")
            iplotStdDev = 0
        if (iMean == 0):
            if (iplotStdDev == 1):
                print("Compute standard deviation")
                MOM_StdDev = np.std(MOMsave, axis=0)
                MOM_StdDev = np.expand_dims(MOM_StdDev, axis=0)
            if (iplotStdDev == 2):
                print("Compute 10, 90 percentile",iMean)
                MOM_StdDev = np.abs(np.percentile(MOMsave, [10,90], axis=0)-MOMmean) # given as positive deviations from the mean value

    print('Plot moment data: ', nr_inst,nt_time,nr_Mom)

    print('iplot_mode', iplot_mode)
    if (iplot_mode == 0) or (iplot_mode == 1):
        i_simMOM = -1
        Simcycler_dict  = cycler_dict(nr_sims)
        #print(Simcycler_dict)
        plt.rc('axes', prop_cycle=Simcycler_dict[Sim_cycler])
        #print(Simcycler_dict,Sim_cycler)
        if (iplot_onlyMOMmean == 0):
            custom_cycler   = cycler(color=[plt.cm.cool(i) for i in np.linspace(0, 1, nr_inst)])
            plt.rc('axes', prop_cycle=custom_cycler) #this call together with preceding cycler definition must be placed before plt.figure call!
        #GCCif (MOM_Panel2Dsims == 0)
        yfigsize = (nrMom_select+iDmean+iSIPplot) * 1.5
        xfigsize = 4 + 3.5*(nrPanels-1)
        fig0 = plt.figure(figsize=(xfigsize,yfigsize), dpi=300)
        ax_vec = fig0.subplots(nrMom_select+iDmean+iSIPplot,nrPanels,squeeze=False)
        if (nrMom_select+iDmean+iSIPplot > 1 or nrPanels > 1):
            if (iPlotextras == 6) or (iBoxM_LWC_N_r_var == 1):
                fig0.subplots_adjust(left=0.22, bottom=0.15, right=0.85, top=0.91, wspace=0.3, hspace=0.3)
            elif ((iHalfDom == 1) and (iExtraLegend != 5)) or ((iEmptydom == 1) and (iExtraLegend != 5)):
                fig0.subplots_adjust(left=0.22, bottom=0.15, right=0.85, top=0.94, wspace=0.3, hspace=0.3)
            else:
                fig0.subplots_adjust(left=0.22, bottom=0.15, right=0.85, top=0.9, wspace=0.1, hspace=0.3)
        if (title is not None):
            if (i_inline == 0):
                plt.suptitle(title)
            else:
                for ax_tmp in ax_vec.flat:
                    ax_tmp.text(0.4,0.75,title, fontsize= 10, transform=ax_tmp.transAxes,  weight='bold',bbox=dict(facecolor='mediumorchid', edgecolor='black', boxstyle = 'square, pad=0.2'))
        #GCCendif /* (MOM_Panel2Dsims == 0) */
        #GCCif (MOM_Panel2Dsims == 1)
        yfigsize = nr_rows * 1.5
        xfigsize = 4 + 3.5*(nr_cols-1)

        #Default settings for iClip = 0
        colwidth = [1]*nr_cols # default: all columns have equal width
            # xfigsize and yfigsize are not changed
        tmin = [0]*nr_cols
        wspace = 0.1
        top = 0.91

        if (iClip == 1):
            colwidth = [1,1]
            xfigsize = 4.0
            tmin = [30,30]
            wspace = 0.2
            top = 0.92
        if (iClip == 2):
            colwidth = [1]
            xfigsize = 2.0
            tmin = [30,30,30]
            top = 0.92
        if (iClip == 3):
            colwidth = [2,1]
            xfigsize = 6.0
            tmin = [0,30]
            wspace = 0.15
        fig0 = plt.figure(figsize=(xfigsize,yfigsize), dpi=300)

        ax_vec = fig0.subplots(nr_rows,nr_cols,squeeze=False,gridspec_kw={'width_ratios': colwidth})
        fig0.subplots_adjust(left=0.22, bottom=0.15, right=0.85, top=top, wspace=wspace, hspace=0.3)
        iDmean = 0
        if (titleMAIN is not None):
            if (iPlotextras == 2 or \
                iPlotextras == 3 or \
                iPlotextras == 4 or \
                iPlotextras == 5):
                plt.suptitle(titleMAIN,va='bottom')  #va='bottom'
            else:
                plt.suptitle(titleMAIN)
        #GCCendif /* (MOM_Panel2Dsims == 1) */
        print(ax_vec.shape)

        plt.set_cmap('Reds')
        #ylabel_vec=['$\lambda_0$','$\lambda_1$','$\lambda_2$','$\lambda_3$']

        str_unit = 'kg'
        if (skal_m is not None):
            str_unit = r'$m_e$'

        #GCCif (MOM_meanTD <= 1)
        ylabel_vec = ['$\lambda_0$ / m$^{-3}$',
                      '$\lambda_1$ / (' + str_unit + ' m$^{-3}$)',
                      '$\lambda_2$ / (' + str_unit + '$^{2}$ m$^{-3}$)',
                      '$\lambda_3$ / (' + str_unit + '$^{3}$ m$^{-3}$)']
        #GCCendif /* (MOM_meanTD == 1) */

        #GCCif (MOM_meanTD == 2)
        ylabel_vec=['$\lambda_0$ / m$^{-2}$',
                    '$\lambda_1$ / (' + str_unit + ' m$^{-2}$)',
                    '$\lambda_2$ / (' + str_unit + '$^{2}$ m$^{-2}$)',
                    '$\lambda_3$ / (' + str_unit + '$^{3}$ m$^{-2}$)']
        #GCCendif /* (MOM_meanTD == 2) */

        #GCCif (MOM_meanTD == 3)
        ylabel_vec=['$\lambda_0$',
                    '$\lambda_1$ / ' + str_unit + '',
                    '$\lambda_2$ / ' + str_unit + '$^{2}$',
                    '$\lambda_3$ / ' + str_unit + '$^{3}$']
        #GCCendif /* (MOM_meanTD == 3) */

    #GCCif (MOM_NORMmass == 1)
    ylabel_vec[1] = '$\lambda_1$ / $\lambda_1(t=0)$'
    #GCCendif /* (MOM_NORMmass == 1) */

    if (iDmean == 1):
        #compute Dmean, must be done before first moment is normalised
        a = MOMmean[:,1]
        if (skal_m is not None):
            a = a * skal_m
        b = MOMmean[:,0]
        c = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
        Dmean = FK.m2r(c, const_mass2rad)*2e6 # diameter in um

    #GCCif (MOM_NORMmass == 1)
    for k in range(nr_inst):
        MOMsave[k,:,1] = MOMsave[k,:,1]/MOMsave[k,0,1] # normalise first moment (= mass)
    #GCCendif /* (MOM_NORMmass == 1) */

    if (itime_legend == 2):
        t_vec_MOMsave = t_vec_MOMsave/60
        if (iSIPplot == 1):
            t_vec_GVplot  = t_vec_GVplot/60

    i_simMOM += 1

    #GCCif (MOM_Panel2Dsims == 0)
    ipc=iPanel[i_simMOM]
    print('-------label--------',label)
    nr_labels = len(label)
    #GCCendif /* (MOM_Panel2Dsims == 0) */

    #GCCif (MOM_Panel2Dsims == 1)
    ipc = iPanel[i_simMOM] % nr_cols
    ipr = iPanel[i_simMOM] // nr_cols
    print('nr_rows,nr_cols',nr_rows,nr_cols)
    print('iPanel[i_simMOM],i_simMOM,ipr, ipc', iPanel[i_simMOM], i_simMOM, ipr, ipc)
    label = label_vec[iPanel[i_simMOM]]
    if label is None:
        nr_labels = 0
        print('no legend')
    else:
        nr_labels = len(label)
        print('-------label--------',label)
    #GCCendif /* (MOM_Panel2Dsims == 1) */

    if (iDmean == 1):
        #additional plot with Dmean-evolution
        #print('Dmean',Dmean)
        if (i_simMOM < nr_labels):
            print('set label:',i_simMOM, label[i_simMOM])
            #ax_vec[0].set_label(label[i_simMOM])
            ax_vec[0,ipc].plot(t_vec_MOMsave,Dmean, markevery=10,label=label[i_simMOM])
        else:
            ax_vec[0,ipc].plot(t_vec_MOMsave,Dmean, markevery=10)
        ax_vec[0,ipc].get_xaxis().set_ticklabels([])
        ax_vec[0,ipc].set_ylabel('$D_\mathrm{mean}$ / $\mu$m')
        ax_vec[0,ipc].set_xlim([0,math.ceil(max(t_vec_MOMsave))])
        ax_vec[0,ipc].set_ylim([0,100])
        #ax_vec[0,ipc].minorticks_on()

    #GCCif (addSIPplot > 0)
    if (iSIPplot == 1):
        ax_vec[iDmean,ipc].plot(t_vec_GVplot,nr_SIPs_mean[:], markevery=10)
        print('nr_SIPs',nr_SIPs_mean[:])
        ax_vec[iDmean,ipc].get_xaxis().set_ticklabels([])
        ax_vec[iDmean,ipc].set_ylabel('$n_\mathrm{SIP}$')
        ax_vec[iDmean,ipc].set_xlim([0,math.ceil(max(t_vec_GVplot))])
        #ax_vec[iDmean,ipc].set_ylim([0,150])
        ax_vec[iDmean,ipc].minorticks_on()
    #GCCendif /* (addSIPplot > 0) */

    for iMom in range(nrMom_select):
        iiMom = iMom_select[iMom]
        print('*****************************',iiMom)

        #print('t_vec_MOMsave', t_vec_MOMsave)
        #print('MOMsave[k,:,iiMom]', MOMsave[k,:,iiMom])
        #GCCif (MOM_Panel2Dsims == 0)
        ipr=iMom+iDmean+iSIPplot
        if (iplot_mode == 0) or (iplot_mode == 1):
            ax_vec[ipr,ipc].set_ylabel(ylabel_vec[iiMom])
        #GCCendif /* (MOM_Panel2Dsims == 0) */

        if (iplot_onlyMOMmean == 0):
            for k in range(0,nr_inst):
                ax_vec[ipr,ipc].plot(t_vec_MOMsave,MOMsave[k,:,iiMom],'r:')
            if (iplotStdDev == 0):
                line,=ax_vec[ipr,ipc].plot(t_vec_MOMsave,MOMmean[:,iiMom],'k:',linewidth=2.0)
            else:
                #print('AA', MOMmean[:,iiMom].shape,MOM_StdDev[:,iiMom].shape)
                ax_vec[ipr,ipc].errorbar(t_vec_MOMsave,MOMmean[:,iiMom], fmt='k:',linewidth=2.0, markevery=10, yerr=np.squeeze(MOM_StdDev[:,:,iiMom]), errorevery=10)
        else:
            if (iplotStdDev == 0):
                #print('MOMmean[:,iiMom]: ', MOMmean[:,iiMom])
                line,=ax_vec[ipr,ipc].plot(t_vec_MOMsave,MOMmean[:,iiMom], markevery=10)
            else:
                #print(ip,i_simMOM)
                #print(MOM_StdDev[:,iiMom]/MOMmean[:,iiMom])
                line,=ax_vec[ipr,ipc].errorbar(t_vec_MOMsave,MOMmean[:,iiMom],  markevery=10, yerr=np.squeeze(MOM_StdDev[:,:,iiMom]), errorevery=10+i_simMOM)
            if (i_simMOM < nr_labels):
                ax_vec[ipr,ipc].set_label(label[i_simMOM])

        #GCCif (MOM_Panel2Dsims == 0)
        if (iplot_mode == 0) or (iplot_mode == 3):
            for ipii in range(nrPanels):
                #ax_vec[ipr,ipii].minorticks_on()
                if (iEmptydom == 1):
                    print('---------------------------------------------------------------')
                    print('choose EmptyDom plot settings', iiMom)
                    if (iiMom == 0):
                        ax_vec[ipr,ipii].set_ylim([0,1.1e6])
                        #ax_vec[ipr,ipii].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
                        ax_vec[ipr,ipii].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                        #ax_vec[ipr,ipii].set_ylim([0,.1e5])
                    if (iiMom == 2 or iiMom == 3):
                        ax_vec[ipr,ipii].semilogy()
                        if (iiMom == 2):
                            print('YES!!!!!!!!!!')
                            ax_vec[ipr,ipii].set_ylim([1e-16,1e-10])
                    ax_vec[ipr,ipii].set_xlim([0,math.ceil(max(t_vec_MOMsave))])
                    ax_vec[ipr,ipii].xaxis.set_ticks([0,20,40,60,80,100,120])
                    ax_vec[ipr,ipii].minorticks_on()
                if (iHalfDom == 1):
                    print('choose HalfDom plot settings')
                    if (iDmean == 1):
                        ax_vec[0,ipii].set_ylim([0,50])
                        ax_vec[0,ipii].yaxis.set_ticks([0,10,20,30,40,50])
                        ax_vec[0,ipii].minorticks_on()
                        ax_vec[0,ipii].grid(b=True)

                    if (iiMom != 1):
                        ax_vec[ipr,ipii].semilogy()
                    else:
                        ax_vec[ipr,ipii].set_ylim([0,None])

                    if (iiMom == 0):
                        ax_vec[ipr,ipii].set_ylim([1e6,1e9])
                        ax_vec[ipr,ipii].yaxis.set_ticks([1e6,1e7,1e8,1e9])
                    #if (iiMom == 1):
                        #ax_vec[ipr,ipii].set_ylim([0,3e-4])

                    if (iiMom == 2):
                        ax_vec[ipr,ipii].yaxis.set_ticks([1e-15,1e-13,1e-11,1e-9,1e-7])

                    ax_vec[ipr,ipii].set_xlim([0,60])
                    ax_vec[ipr,ipii].xaxis.set_ticks([0,10,20,30,40,50,60])
                    ax_vec[ipr,ipii].minorticks_on()
                    ax_vec[ipr,ipii].grid(b=True)

                if (iBoxM_LWC_N_r_var == 1):
                    print('choose BoxModel LWC_N_r_var plot settings')
                    if (iDmean == 1):
                        print('set Dmean options')
                        ax_vec[0,ipii].set_ylim([0,500]) # Dmean-plot
                        ax_vec[0,ipii].yaxis.set_ticks([0,100,200,300,400,500])
                        ax_vec[0,ipii].set_xlim([0,100])
                        ax_vec[0,ipii].xaxis.set_ticks([0,20,40,60,80,100])
                        ax_vec[0,ipii].minorticks_on()
                        ax_vec[0,ipii].grid(b=True)

                    #ax_vec[ipr,ipii].set_xlim([0,60])
                    #ax_vec[ipr,ipii].xaxis.set_ticks([0,10,20,30,40,50,60])
                    ax_vec[ipr,ipii].set_xlim([0,100])
                    ax_vec[ipr,ipii].xaxis.set_ticks([0,20,40,60,80,100])
                    ax_vec[ipr,ipii].semilogy()
                    if (iiMom == 0):
                        ax_vec[ipr,ipii].set_ylim([1e3,1e9])
                        ax_vec[ipr,ipii].yaxis.set_ticks([1e3,1e5,1e7,1e9])
                    if (iiMom == 2):
                        ax_vec[ipr,ipii].set_ylim([1e-15,1e-6])
                        ax_vec[ipr,ipii].yaxis.set_ticks([1e-15,1e-12,1e-9,1e-6])
                    ax_vec[ipr,ipii].minorticks_on()
                    ax_vec[ipr,ipii].grid(b=True)


                if ((1-iEmptydom)*(1-iHalfDom)*(1-iBoxM_LWC_N_r_var) == 1):
                    print('choose generic settings')

                    ax_vec[ipr,ipii].minorticks_on()
                    ax_vec[ipr,ipii].grid(b=True)
                    if (iiMom != 1):
                        ax_vec[ipr,ipii].semilogy()
                    if (iiMom == 0):
                        ax_vec[ipr,ipii].set_ylim([1e5,4e8])
                        ax_vec[ipr,ipii].yaxis.set_ticks([1e5,1e6,1e7,1e8])
                    if (iiMom == 1):
                        ax_vec[ipr,ipii].set_ylim([0,3e-4])

                    ax_vec[ipr,ipii].set_xlim([0,60])
                    ax_vec[ipr,ipc].xaxis.set_ticks([0,10,20,30,40,50,60])

                    ax_vec[ipr,ipii].minorticks_on()
                    ax_vec[ipr,ipii].yaxis.set_ticks_position('both')
                    # ax_vec[ipr,ipii].tick_params(labelright = True)
                print('max t:  ',math.ceil(max(t_vec_MOMsave)))
                if (ipii > 0):  # include y-axis only in left most plot
                    ax_vec[ipr,ipii].get_yaxis().set_ticklabels([])

                if (iMom == nrMom_select-1):  # include x-axis only in lower most plot
                    ax_vec[ipr,ipii].set_xlabel(xlabel_Mom[itime_legend])
                else:
                    ax_vec[ipr,ipii].get_xaxis().set_ticklabels([])
    # end iteration over moments
        #GCCendif /* (MOM_Panel2Dsims == 0) */

    #GCCif (MOM_Panel2Dsims == 1)
    for ipr in range(nr_rows):
        for ipc in range(nr_cols):
            iiPanel = ipr * nr_cols + ipc
            if (iplot_mode == 0) or (iplot_mode == 1):
                print('================================================================')
                if (iPlotextras == 4):
                    print('iPlotextras 4')
                    if (iiMom == 0):
                        ax_vec[ipr,ipc].set_ylim([4e7,3e8])
                elif (iPlotextras == 5):
                    print('iPlotextras 5')
                    ax_vec[ipr,ipc].semilogy()
                    if (iiMom == 0):
                        ax_vec[ipr,ipc].set_ylim([1e4,4e8])
                        ax_vec[ipr,ipc].yaxis.set_ticks([1e4,1e5,1e6,1e7,1e8])
                else:
                    ax_vec[ipr,ipc].semilogy()
                    if (iiMom == 0):
                        ax_vec[ipr,ipc].set_ylim([1e5,4e8])
                        ax_vec[ipr,ipc].yaxis.set_ticks([1e5,1e6,1e7,1e8])

                ax_vec[ipr,ipc].minorticks_on()
                ax_vec[ipr,ipc].grid(b=True)

                ax_vec[ipr,ipc].set_xlim([tmin[ipc],60])
                ax_vec[ipr,ipc].xaxis.set_ticks(np.arange(tmin[ipc],61,10))

                if (ipc == 0):  # include y-axis only in left most plot
                    ax_vec[ipr,ipc].set_ylabel(ylabel_vec[iiMom])
                else:
                    ax_vec[ipr,ipc].get_yaxis().set_ticklabels([])
                    ax_vec[ipr,ipc].get_yaxis().set_ticklabels([],minor=True)

                if (ipr == nr_rows-1):  # include x-axis only in lower most plot
                    ax_vec[ipr,ipc].set_xlabel(xlabel_Mom[itime_legend])
                else:
                    ax_vec[ipr,ipc].get_xaxis().set_ticklabels([])
            if (i_inline == 0):
                ax_vec[ipr,ipc].set_title(title_vec[iiPanel])
            else:
                ax_vec[ipr,ipc].text(xpos[iiPanel], zpos[iiPanel], title_vec[iiPanel], fontsize= fontsize_inline, transform=ax_vec[ipr,ipc].transAxes,  weight='bold', bbox=dict(facecolor='mediumorchid', edgecolor='black', boxstyle = 'square, pad=0.2'))

            if (iplot_mode == 0) or (iplot_mode == 3):
                print('================================================================')
                text = text_vec[iiPanel]
                print('iiPanel,text: ', iiPanel, text)
                if text is not None:
                    print('len(text),type(text)',len(text),type(text))
                    if type(text) is str:
                        ax_vec[ipr,ipc].text(0.05, 0.08, text, horizontalalignment='left',verticalalignment='bottom', transform=ax_vec[ipr,ipc].transAxes,bbox=dict(facecolor='white',pad=2),fontsize=6)
                    else:
                        ha='left'
                        va='bottom'
                        if (len(text) >= 4):
                            ha=text[3]
                        if (len(text) == 5):
                            va=text[4]
                        print('ha,va',ha,va)
                        ax_vec[ipr,ipc].text(text[1],text[2], text[0], horizontalalignment=ha,verticalalignment=va, transform=ax_vec[ipr,ipc].transAxes,bbox=dict(facecolor='white',pad=2),fontsize=6)
                if iPanelAnnotate is not None:
                    delh = 0.0
                    if (iClip == 3) and (ipc == 1):
                        delh = 0.05
                    ax_vec[ipr,ipc].text(iPanelAnnotate[0]-delh,iPanelAnnotate[1], string.ascii_lowercase[iiPanel]+')', transform=ax_vec[ipr,ipc].transAxes,  weight='bold')
            #GCCif (IREF == 1)
                # read moment data of reference solution
                [timeREF, MomentsREF] = REF.defineMoments_0D()
                if (itime_legend == 2): timeREF = timeREF/60
                iiMom = iMom_select[iMom]
                line,=ax_vec[ipr,ipc].plot(timeREF,MomentsREF[:,iiMom], color='k', label='BIN')
                print('aa',line.get_c(),line.get_ls())
            #GCCendif /* (IREF == 1) */

            #GCCif (IREF == 9)
                color_tmp = 'k'
                label_tmp = 'BIN'
                isimREF = 0
                line,=ax_vec[ipr,ipc].plot(t_vecREF[:ntREF_vec[isimREF]],MOMsaveREF[isimREF,iiMom,:ntREF_vec[isimREF]], color=color_tmp, label=label_tmp, linestyle=lsREF)
            #GCCendif /* (IREF == 9) */

                if (iExtraLegend == 1):
                    leg_extra = None
                    if (ipr == 1) and (ipc == 0):
                        label_extra = ['24', '49', '98', '197', '296', '497', '993']
                        fs=labelMoments_dict.get('fs',5)
                        leg_extra = ax_vec[ipr,ipc].legend(label_extra, fontsize=fs, loc='lower center', title='$N_\mathrm{SIP,GB}$',title_fontsize = fs*1.3,ncol=2) #,title_fontsize=fs
                        bb = leg_extra.get_bbox_to_anchor().inverse_transformed(ax_vec[ipr,ipc].transAxes)
                        # Change to location of the legend.
                        xOffset = 0.1
                        bb.x0 += xOffset
                        bb.x1 += xOffset
                        leg_extra.set_bbox_to_anchor(bb, transform = ax_vec[ipr,ipc].transAxes)

                    if (ipr == 3) and (ipc == 0):
                        #label_extra = ['24','124', '248','197','993','1987']
                        label_extra = ['98', '496', '993', '497', '2484', '4968']
                        fs=labelMoments_dict.get('fs',5)
                        leg_extra = ax_vec[ipr,ipc].legend(label_extra, fontsize=fs, loc='lower center', title='$N_\mathrm{SIP,GB}$',title_fontsize = fs*1.3,ncol=2) #,title_fontsize=fs
                        # Get the bounding box of the original legend
                        bb = leg_extra.get_bbox_to_anchor().inverse_transformed(ax_vec[ipr,ipc].transAxes)
                        # Change to location of the legend.
                        xOffset = 0.15
                        bb.x0 += xOffset
                        bb.x1 += xOffset
                        leg_extra.set_bbox_to_anchor(bb, transform = ax_vec[ipr,ipc].transAxes)
                if (iExtraLegend == 2):
                    leg_extra = None
                    if (ipr == 0) and (ipc == 0):
                        label_extra = ['24','49', '98','197','296','497','993']
                        fs=labelMoments_dict.get('fs',5)
                        leg_extra = ax_vec[ipr,ipc].legend(label_extra, fontsize=fs, loc='lower left', title='$N_\mathrm{SIP,GB}$',ncol=2)
                        bb = leg_extra.get_bbox_to_anchor().inverse_transformed(ax_vec[ipr,ipc].transAxes)
                        # Change to location of the legend.label_vec
                        xOffset = 0.2
                        bb.x0 += xOffset
                        bb.x1 += xOffset
                        leg_extra.set_bbox_to_anchor(bb, transform = ax_vec[ipr,ipc].transAxes)

                if (iBin_noLegend == 0):
                    print('mit Bin curve in legend')
                    label=label_vec[iiPanel]+['BIN']
                else:
                    print('ohne Bin curve in legend')
                    label=label_vec[iiPanel]

                if (label[0] is not None):
                    print('label[0]', label[0])
                    lc=labelMoments_dict.get('loc','lower left')
                    fs=labelMoments_dict.get('fs', 5)
                    leg_ncol=labelMoments_dict.get('ncol',1)

                    ipr_leg=ipr
                    ipc_leg=ipc
                    #ipr_leg = 2
                    #lc ='lower right'

                    if (leg_ncol == 1):
                        print('A labels:', label[:nr_sims+iaddsimREF])
                        lab_sc = None
                        if (iiPanel == 2) & (iPlotextras == 1):
                            leg_ncol = 2
                        if (iiPanel == 6) & (iPlotextras == 1):
                            lab_sc = 0.2

                        leg_tmp = ax_vec[ipr_leg,ipc_leg].legend(label[:nr_sims+iaddsimREF], fontsize=fs, loc=lc, ncol=leg_ncol,labelspacing = lab_sc)
                        #legend = plt.legend(handles=[one, two, three], title="title", loc=4, fontsize='small', fancybox=True)
                    else:
                        handles, labels = ax_vec[0,ipc_leg].get_legend_handles_labels()
                        print('handles', handles)
                        print('labels', labels)
                        ax_vec[ipr_leg,ipc_leg].legend(flip(handles,leg_ncol),flip(labels,leg_ncol), fontsize=fs, loc=lc, ncol=leg_ncol)
                if (iExtraLegend > 0):
                    if (leg_extra is not None):
                        ax_vec[ipr,ipc].add_artist(leg_extra)

    #GCCendif /* (MOM_Panel2Dsims == 1) */

    if (iplot_mode == 0) or (iplot_mode == 3):
        #GCCif (MOM_Panel2Dsims == 0)
        if text is not None:
            print('len(text),type(text)',len(text),type(text))
            if type(text) is str:
                ax_vec[0,ipc].text(0.96, 0.92, text, horizontalalignment='right',verticalalignment='top', transform=ax_vec[0,ipc].transAxes,bbox=dict(facecolor='none',pad=2),fontsize=6)
            else:
                ha='right'
                if (len(text) == 4):
                    ha=text[3]
                ax_vec[0,ipc].text(text[1],text[2], text[0], horizontalalignment=ha,verticalalignment='top', transform=ax_vec[0,ipc].transAxes,bbox=dict(facecolor='none',pad=2),fontsize=6)

        #GCCif ((IREF == 1) || (IREF == 3))
        # read moments of reference solution
        [timeREF, MomentsREF] = REF.defineMoments_0D()
        print('timeREF A',timeREF)
        if (itime_legend == 2): timeREF = timeREF/60
        print('timeREF B',timeREF)
        for iMom in range(nrMom_select):
            iiMom = iMom_select[iMom]
            ipr=iMom+iDmean+iSIPplot
            print('CAB',iMom,nrMom_select,iiMom,ipr,ipc,iDmean,iSIPplot)
            #GCCif (IREF == 3)
            if (iiMom == 2):
                M2_init=skal_m*skal_m*nplot*dVi
                tau=1./(C_original*M2_init)
                ax_vec[ipr,ipc].plot([tau,tau],[min(MomentsREF[:,iiMom]),max(MomentsREF[:,iiMom])],'k',label='Analytical')
            if (iiMom == 1):
                continue
            #GCCendif
            ax_vec[ipr,ipc].plot(timeREF,MomentsREF[:,iiMom],'k',label='BIN')
        #GCCendif
        #GCCif (IREF == 9)
        if (nrsimREF == 1):
            color_tmp = 'k'
            label_tmp = 'BIN'
            label.append('BIN')
        else:
            color_tmp = None
            label_tmp = None

        for isimREF in range(nrsimREF):
            print('nrsimREF: ', nrsimREF)
            ipc=iPanelREF[isimREF]
            for iMom in range(nrMom_select):
                iiMom = iMom_select[iMom]
                ipr=iMom+iDmean+iSIPplot
                if (iiMom != 1):
                    line,=ax_vec[ipr,ipc].plot(t_vecREF[:ntREF_vec[isimREF]],MOMsaveREF[isimREF,iiMom,:ntREF_vec[isimREF]], color=color_tmp, label=label_tmp, linestyle=lsREF)
                else:
                    #GCCif (MOM_NORMmass == 1)
                    Momtmp=MOMsaveREF[isimREF,iiMom,:ntREF_vec[isimREF]]/MOMsaveREF[isimREF,iiMom,0] # normalise first moment (= mass)
                    #GCCelse
                    Momtmp=MOMsaveREF[isimREF,iiMom,:ntREF_vec[isimREF]]
                    #GCCendif /* (MOM_NORMmass == 1) */
                    line,=ax_vec[ipr,ipc].plot(t_vecREF[:ntREF_vec[isimREF]],Momtmp, color=color_tmp, label=label_tmp, linestyle=lsREF)  #,col_list[ii]+':'
                print('label_tmp', label_tmp)
                ax_vec[ipr,ipc].set_label(label_tmp)
                print('aa',line.get_c(),line.get_ls())
            if (iDmean == 1):
                ax_vec[0,ipc].plot(t_vecREF[:ntREF_vec[isimREF]], DmeanREF[isimREF,:ntREF_vec[isimREF]], color=color_tmp, label=label_tmp, linestyle=lsREF)

        #GCCendif /* (IREF == 9) */

        if (label[0] is not None):
            if (iExtraLegend == 3):
                lc = 'upper right'
                fs = labelMoments_dict.get('fs',6)
                leg_ncol = labelMoments_dict.get('ncol',1)
                ipr_leg = 0
                ipc_leg = 0
                dummy_lines = []
                for b_idx, b in enumerate(['-','--', '-.',(0, (3, 10, 1, 10)), ':']):  # use linestyles order as in C4L5 cycler
                    dummy_lines.append(ax_vec[ipr_leg,ipc_leg].plot([],[], c="black", ls = b)[0])
                legend2 = ax_vec[ipr_leg,ipc_leg].legend([dummy_lines[i] for i in range(5)], ["reg", "WM2D", " nS", "LS", "BIN"], loc='upper left',fontsize=fs,handlelength=4)
                print('legend2 created', ipr_leg, ipc_leg)
                #ax_vec[ipr_leg,ipc_leg].add_artist(legend2) this call is not necessary as the additional legend and the original legend are added in different panels.

            lc=labelMoments_dict.get('loc','lower left')
            fs=labelMoments_dict.get('fs',6)
            leg_ncol=labelMoments_dict.get('ncol',1)
            ipr_leg=labelMoments_dict.get('ipr',0)
            ipc_leg=labelMoments_dict.get('ipc',0)

            kwargs = {}
            #if (iPlotextras == 20):
                #kwargs = {"bbox_to_anchor": (0.02, 0.2)}

            if (leg_ncol == 1):
                print('B labels:', label[:nr_sims+iaddsimREF], ipr_leg, ipc_leg)
                ax_vec[ipr_leg,ipc_leg].legend(label[:nr_sims+iaddsimREF], fontsize=fs, loc=lc,**kwargs)
            else:
                handles, labels = ax_vec[0,ipc_leg].get_legend_handles_labels()
                print('C handles', handles)
                print('C labels', labels)
                ax_vec[ipr_leg,ipc_leg].legend(flip(handles,leg_ncol),flip(labels,leg_ncol), fontsize=fs, loc=lc, ncol=leg_ncol)

        if (iExtraLegend == 5):
            for ipc in range(3):
                vert = 0.2 #0.92
                ax_vec[0,ipc].text(0.5, vert, ["AON-WM2D", "AON-regular", "BIN"][ipc], horizontalalignment='right',verticalalignment='top', transform=ax_vec[0,ipc].transAxes,bbox=dict(facecolor='none',pad=2),fontsize=labelMoments_dict.get('fs',6))

        #GCCendif /* (MOM_Panel2Dsims == 0) */
        plt.savefig(fp_out+'Mom.png', format='png', bbox_inches='tight', dpi=300)
        plt.savefig(fp_out+'Mom.pdf', format='pdf', bbox_inches='tight', dpi=300)

        plt.close()
    print('end of plot function')
# end function PlotMoments -----------------------------------------------------------------------------------

def PlotMomentsFinal(MOMFinal,fp_out='', skal_m=None):
    # MOMFinal=np.zeros([nr_sims,nr_Mom])

#GCCif (IREF == 1)
    #Wang reference data
    [timeREF, MomentsREF] = REF.defineMoments_0D()
    MomRefFinal = MomentsREF[-1,:]
#GCCendif /* (IREF == 1) */

#GCCif (IREF == 9)
    ntREF_vec=np.zeros(nrsimREF, dtype='int')
    print('nrsimREF', nrsimREF)
    MomRefFinal = np.zeros([nrsimREF,4])
    for isimREF in range(nrsimREF):
        ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF,rgridREF,mgridREF=REF.get_RefMetaData(isimREF=isimREF)
        #print('REF meta data', ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF)
        ntREF_vec[isimREF] = ntREF
        if (ntREF > ntREFmax):
            print('increase ntREFmax', ntREF, ntREFmax)

        t_vecREF = TsimREF/(ntREF-1)*np.arange(ntREF)

        MOMtmp = REF.get_RefProfileData(nzREF,ntREF,isimREF=isimREF)  #function output : np.zeros([4,nt,nz])
        MOMtmp = np.mean(MOMtmp, axis=2)

        nt_time_pick = -1
        if (iPlotextras == 10):
            nt_time_pick = [18,18][isimREF]
        MomRefFinal[isimREF,:] = MOMtmp[:,nt_time_pick]
#GCCendif /* (IREF == 9) */

    [nr_sims,nr_Mom]=MOMFinal.shape
    iMom_select = iMom_select_const.copy()
    nrMom_select= len(iMom_select)

    if (skal_m is not None):
        for i in range(4):
            MOMFinal[:,i] = MOMFinal[:,i] / skal_m**i

    str_unit = 'kg'
    if (skal_m is not None):
        str_unit = r'$m_e$'

    ylabel_vec = ['$\lambda_0(t=1\:$h) / m$^{-3}$',
                    '$\lambda_1(t=1\:$h)  / (' + str_unit + ' m$^{-3}$)',
                    '$\lambda_2(t=1\:$h)  / (' + str_unit + '$^{2}$ m$^{-3}$)',
                    '$\lambda_3(t=1\:$h)  / (' + str_unit + '$^{3}$ m$^{-3}$)']

    for iMom in iMom_select:
        if (isavedata_filename is not None):
            datafile = open(isavedata_filename+'_{0:01}.dat'.format(iMom),'w')
            datafile.write(titleMAIN+'\n')
            datafile.write(str(len(iperm)) + '   '+ str(nr_seriestypes)+'\n')

        yfigsize = 4.5
        xfigsize = 3.0
        Simcycler_dict  = cycler_dict(1)

        if type(Sim_cycler) is not list:
            plt.rc('axes', prop_cycle=Simcycler_dict[Sim_cycler])
            fig = plt.figure(figsize=(xfigsize,yfigsize), dpi=300)
            ax = fig.subplots(2,1)

        if type(Sim_cycler) is list:
            fig = plt.figure(figsize=(xfigsize,yfigsize), dpi=300)
            ax = fig.subplots(2,1)
            for iii,SC in enumerate(Sim_cycler):
                ax[iii].set_prop_cycle(Simcycler_dict[SC])

        if (titleMAIN is not None):
            plt.suptitle(titleMAIN,va='bottom')  #va='bottom'

        for iseries in iperm:
            ia = iseries_indices[iseries]
            ie = iseries_indices[iseries+1]
            nr = ie-ia
            xgrid = iseries_shift[iseries] + np.arange(nr)
            line,=ax[iseriestype[iseries]].plot(xgrid,MOMFinal[iseries_indices[iseries]:iseries_indices[iseries+1],iMom],label=text_vec[iseries],marker='*')
            #get current linestyle and colour!
            ls_current = line.get_ls()
            col_current = line.get_c()

            if (isavedata_filename is not None):
                datafile.write(" ".join("{}".format(x) for x in [iseries,iseriestype[iseries],ls_current,col_current])+'\n')
                datafile.write(text_vecFILE[iseries]+'\n')
                datafile.write(" ".join("{}".format(x) for x in xgrid)+'\n')
                #datafile.write(" ".join("'{}'".format(x) for x in series_vecs[iseriestype[iseries]])+'\n')
                datafile.write(" ".join("{}".format(x) for x in MOMFinal[iseries_indices[iseries]:iseries_indices[iseries+1],iMom])+'\n')

        for ist in range(nr_seriestypes):
            #GCCif (IREF == 1)
            if (isavedata_filename is not None):
                datafile.write(" ".join("{}".format(x) for x in [ist,'iref1',':','k'])+'\n')
                datafile.write(" ".join("{}".format(x) for x in [0,len(series_vecs[ist])-1])+'\n')
                datafile.write(" ".join("{}".format(x) for x in [MomRefFinal[iMom],MomRefFinal[iMom]])+'\n')
            ax[ist].plot([0,len(series_vecs[ist])-1], [MomRefFinal[iMom],MomRefFinal[iMom]],'k:')
            #GCCendif /* (IREF == 1) */
            #GCCif (IREF == 9)
            if (isavedata_filename is not None):
                datafile.write(" ".join("{}".format(x) for x in [ist,'iref9',':','k'])+'\n')
                datafile.write(" ".join("{}".format(x) for x in [0,len(series_vecs[ist])-1])+'\n')
                datafile.write(" ".join("{}".format(x) for x in [MomRefFinal[ist,iMom],MomRefFinal[ist,iMom]])+'\n')
            ax[ist].plot([0,len(series_vecs[ist])-1], [MomRefFinal[ist,iMom],MomRefFinal[ist,iMom]],'k:')
            print('REF value: ', MomRefFinal[ist,iMom])
            #GCCendif /* (IREF == 9) */
            if (isavedata_filename is not None):
                datafile.write(" ".join("{}".format(x) for x in range(len(series_vecsFILE[ist])))+'\n')
                datafile.write(" ".join("{}".format(x) for x in series_vecsFILE[ist])+'\n')
            ax[ist].xaxis.set_ticks(np.arange(len(series_vecs[ist])))
            ax[ist].set_xticklabels(series_vecs[ist])

            ax[ist].semilogy()
            ax[ist].set_ylabel(ylabel_vec[iMom])
            if (iMom == 0):
                ax[ist].set_ylim([1e6,1e8])
                ax[ist].yaxis.set_ticks([1e6,1e7,1e8])

                if (iPlotextras == 10):
                    if (ist == 0):
                        ax[ist].set_ylim([1e2,5e4])
                        ax[ist].yaxis.set_ticks([1e2,1e3,1e4])
                    if (ist == 1):
                        ax[ist].set_ylim([0.7e8,1.5e8])
                        ax[ist].set_yscale('linear')
                if (iPlotextras == 11):
                    ax[ist].set_ylim([1e5,5e7])
                    ax[ist].yaxis.set_ticks([1e5,1e6,1e7])
                    if (ist == 0):
                        ax[ist].text(-0.1,1.8e7,'$dt_\mathrm{Wang}$ in s:',fontsize= 8)
                        for item,text in enumerate(['10', '10', '10', '2', '1', '0.1']):
                            ax[ist].text(item,1e7,text,fontsize= 8,ha='center')

            ax[ist].legend(fontsize=6,handlelength=4)

        #change for top row only
        ax[0].tick_params(axis = 'x', bottom = False, top = True, labelbottom = False, labeltop = True)

        if (isavedata_filename is not None):
            datafile.close()

        plt.savefig(fp_out+"MomFinal_{0:01}.png".format(iMom), format='png', bbox_inches='tight')
        plt.savefig(fp_out+"MomFinal_{0:01}.pdf".format(iMom), format='pdf', bbox_inches='tight')
# end function PlotMomentsFinal -----------------------------------------------------------------------------------

def PlotMomentsTimeCross(TimesCross,fp_out=''):
    # TimesCross = np.zeros([nr_sims,2,3])
    [nr_sims,nr_MomTC,nr_thr]=TimesCross.shape
    print('nr_sims,nr_MomTC,nr_thr',nr_sims,nr_MomTC,nr_thr)
#GCCif (IREF == 9)
    ntREF_vec=np.zeros(nrsimREF, dtype='int')
    print('nrsimREF', nrsimREF)
    TimesCrossREF= np.zeros([nrsimREF,nr_MomTC,nr_thr])
    #MOMsaveREF = np.zeros([nrsimREF,4,ntREFmax])
    for isimREF in range(nrsimREF):
        ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF,rgridREF,mgridREF=REF.get_RefMetaData(isimREF=isimREF)
        #print('REF meta data', ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF)
        ntREF_vec[isimREF] = ntREF
        if (ntREF > ntREFmax):
            print('increase ntREFmax', ntREF, ntREFmax)

        t_vecREF = TsimREF/(ntREF-1)*np.arange(ntREF)

        MOMtmp = REF.get_RefProfileData(nzREF,ntREF,isimREF=isimREF)  #function output : np.zeros([4,nt,nz])
        MOMtmp = np.mean(MOMtmp, axis=2)
        indices = np.zeros([2,3],dtype='int')
        iMom = 0
        iiMom = 0
        for ithr in range(3):
            indices[iiMom,ithr] = np.where(MOMtmp[iMom,:].reshape((ntREF)) < thresholds_TimeCross[iiMom,ithr])[0][0]
        iMom = 2
        iiMom = 1
        for ithr in range(3):
            indices[iiMom,ithr] = np.where(MOMtmp[iMom,:].reshape((ntREF)) > thresholds_TimeCross[iiMom,ithr])[0][0]
        TimesCrossREF[isimREF,:,:]=t_vecREF[indices]
        print(isimREF,t_vecREF[indices])
#GCCendif /* (IREF == 9) */

    for iiMom in range(nr_MomTC):
        yfigsize = 3.0
        xfigsize = 3.5
        Simcycler_dict  = cycler_dict(1)
        plt.rc('axes', prop_cycle=Simcycler_dict[Sim_cycler])
        for ithr in range(3):
            fig = plt.figure(figsize=(xfigsize,yfigsize), dpi=300)
            ax = fig.subplots(1,1)
            for iseries in range(nr_sims//4):
                ia = iseries*4
                ie = iseries*4 + 4
                nr = ie-ia
                xgrid = np.arange(nr)
                label_text = None
                if (iseries < 4):
                    label_text = text_vec[iseries]
                ax.plot(xgrid,TimesCross[ia:ie,iiMom,ithr]/60,label=label_text)
                if (iseries % 4 == 3):
                    iseriesREF = iseries//4
                    ia = iseriesREF*4
                    ie = iseriesREF*4 + 4
                    nr = ie-ia
                    xgrid = np.arange(nr)
                    label_text = None
                    if (iseriesREF == 0):
                        label_text = text_vec[4]
                    ax.plot(xgrid,TimesCrossREF[ia:ie,iiMom,ithr]/60,label=label_text)

            ax.xaxis.set_ticks(xgrid)
            ax.set_xticklabels(xticks_vec)
            ax.set_ylabel('$T_\mathrm{cross}$ / min')
            legend_save = ax.legend(fontsize=6,loc='upper left')

            if (iPlotextras == 11):
                # Extra text lower left corner
                ax.text(-0.03, -0.04, 'LWCvar\nDNCvar', horizontalalignment='right',verticalalignment='top', transform=ax.transAxes) #fontsize=6
                # add extra legend for line styles

                lc = 'upper center'
                fs = 6
                dummy_lines = []
                for b_idx, b in enumerate(['-','--', ':']):  # use linestyles order as in C5L3 cycler
                    dummy_lines.append(ax.plot([],[], c="black", ls = b)[0])
                legend2 = ax.legend([dummy_lines[i] for i in range(3)], ["DNCvar @ LWC$_\mathrm{init}$ fix", "LWCvar @ DNC$_\mathrm{init}$ fix", "LWCvar @ $r_\mathrm{init}$ fix"], loc='lower left',fontsize=fs)
                print('legend2 created')
                ax.add_artist(legend_save)

            iMom=[0, 2][iiMom]
            filename = fp_out+"Tcross_Moment{0:01}_iTresh{1:01}".format(iMom,ithr)
            plt.savefig(filename+".png", format='png', bbox_inches='tight')
            plt.savefig(filename+".pdf", format='pdf', bbox_inches='tight')
# end function PlotMomentsTimeCross -----------------------------------------------------------------------------------


def PlotGV(nEK_sip_plot, mEK_sip_plot, nr_SIPs_plot, t_vec_GVplot, fp_out='', iMultipleSims=0, nr_inst_vec=None, label=[None], title=None, ilastGVonly=None, iTimes_select=None, V=None, outfallingGV=0, skal_m=None):

# Erstelle GV-Plot

#GCCif (IREF == 9)
    if (outfallingGV == 0):
        ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF,rgridREF,mgridREF=REF.get_RefMetaData()
        print('REF meta data', ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF)
        t_vecREF=TsimREF/(ntREF-1)*np.arange(ntREF)
        z_vecREF=dzREF*np.arange(nzREF)
        #print(t_vecREF)
        SizeDistrREF=REF.get_RefSDdata(ntREF,nzREF,nbinREF) #function output : np.zeros(nt,nz+1,nbin)
        SizeDistrREF=np.mean(SizeDistrREF[:,0:nzREF,:],axis=1) # SizeDistrREF contains averaged SD shape: nt,nbin
#GCCendif /* (IREF == 9) */
#GCCif (IREF == 1)
    if (outfallingGV == 0):
        #reference solution given every 600s
        rgridREF,SizeDistrREF,nbinREF = REF.read_SDdata_Wang()
#GCCendif /* (IREF == 1) */

    if (iMultipleSims == 0):
        nEK_sip_plot = np.expand_dims(nEK_sip_plot, axis=0) # add a dummy first dimension with length 1
        mEK_sip_plot = np.expand_dims(mEK_sip_plot, axis=0) # add a dummy first dimension with length 1
        nr_SIPs_plot = np.expand_dims(nr_SIPs_plot, axis=0) # add a dummy first dimension with length 1

    if (outfallingGV >= 1):
        nEK_sip_plot = np.expand_dims(nEK_sip_plot, axis=2) # add a dummy first dimension with length 1
        mEK_sip_plot = np.expand_dims(mEK_sip_plot, axis=2) # add a dummy first dimension with length 1
        nr_SIPs_plot = np.expand_dims(nr_SIPs_plot, axis=2) # add a dummy first dimension with length 1

    [nr_sims,nr_inst,nr_GVplot,nr_SIPs]= nEK_sip_plot.shape

#GCCif (ICOMPPLOT == 1)
    #GCCif (COLUMN  == 0)
    nz = 1
    #GCCendif /*(COLUMN == 0)*/
    V = np.array([nz * dV])
#GCCelse
    if (outfallingGV == 0):
        if (V.size != nr_sims):
            print('provide total volumina of each simulation')
    else:
        V = np.ones([nr_sims])
#GCCendif /*(ICOMPPLOT == 1)*/


    if (nr_inst_vec is None):
        nr_inst_vec=np.array([nr_inst])
        print('nr_inst_vec',nr_inst_vec)
    print('Plotte GV: ',nr_sims, nr_inst,nr_GVplot,nr_SIPs)
    #print('ytitle')
    #print(ytitle)
# determine at which times size distributions are plotted
    #iGV_select = range(0,nr_GVplot)
    iGV_select = list(np.arange(0,nr_GVplot,3))
    #print(iGV_select,ilastGVonly)
    if (ilastGVonly == 1): iGV_select=[nr_GVplot-1]
    if (iTimes_select) is not None: iGV_select = [int(iT/TintervalGV) for iT in iTimes_select]  # assumes that SIP data is in any case stored every 200s.
    if (outfallingGV >= 1):
        iGV_select = [0]

    iGV_select=FK.as_list(iGV_select) # make sure that iGV_select is a list, which is iterable
    nrGV_select=len(iGV_select)

    print('iGV_select', iGV_select)

    skal_g = 1e3  # in the default case, convert from kg/m3 to g/cm3
    if (outfallingGV >= 1):
        ytitle_select = '$g_{\mathrm{ln} r}$ / (g / cm$^2$)'
        ylimGV_select = (1e-1,5e1)
        xlimGV_select = (.8e2,2e3)
        skal_g = 1e1 # convert from kg/m2 to g/cm2
    else:
        ytitle_select = ytitle
        ylimGV_select = ylimGV
        xlimGV_select = xlimGV

    suffix = ['', 'out', 'in'][outfallingGV]
#define plot styles
    Timecycler_dict = cycler_dict(nrGV_select)
    #Simcycler_dict  = cycler_dict(nr_sims+iaddsimREF)
    Simcycler_dict  = cycler_dict(nr_sims)

    #choose plot style for specific plot
    if ((nr_sims+iaddsimREF) == 1):
        custom_cycler = ( Timecycler_dict[Time_cycler])
    else:
        if (nrGV_select == 1):
            custom_cycler = ( Simcycler_dict[Sim_cycler] )
        else:
            custom_cycler = ( Timecycler_dict[Time_cycler] * Simcycler_dict[Sim_cycler] )
    plt.rc('axes', prop_cycle=custom_cycler)  #this call together with preceding cycler definition must be placed before plt.figure call!

    plt.figure(figsize=(4,3), dpi=500)
    #fig, ax = plt.subplots(111)
    plt.xlabel('Radius / '+ r'$\mu$' +'m')
    plt.ylabel(ytitle_select)
    plt.xscale(xscale)
    plt.yscale("log")
    #plt.xlim((1e-1,2e3))
    #plt.ylim((1e-4,1e1))
    #xlimGV=(10,200)
    plt.xlim(xlimGV_select)
    plt.ylim(ylimGV_select)
    if (xticks is not None): plt.xticks(xticks)

    print(nr_SIPs_plot.shape)

    for i_plot in iGV_select:
        for i_sim in range(0,nr_sims):
            print('Plotte GV ',i_plot,' von  Sims ',  i_sim,' zum Zeitpunkt: ',t_vec_GVplot[i_plot] )
            m1D_ins_mean=np.zeros(nplot)
            g1D_ins_mean=np.zeros(nplot)
            nr_inst= nr_inst_vec[i_sim]
            for i_inst in range(0,nr_inst):
                #print(i_inst)
                nr_SIPs=nr_SIPs_plot[i_sim,i_inst,i_plot]
                #print("i_sim,i_inst,i_plot,nr_SIPs", i_sim,i_inst,i_plot,nr_SIPs)
    #GCCif (DISCRETE == 0)
                #print('call MapSIPBin, inst=: ', i_inst)
                print('nr_SIPs,i_plot,i_sim,i_inst ', nr_SIPs,i_plot,i_sim,i_inst)
                if (nr_SIPs > 0):
                    [nEK_bin_plot,mEK_bin_plot,summe2,mdelta_plot,z]=FK.MapSIPBin(nEK_sip_plot[i_sim,i_inst,i_plot,0:nr_SIPs],mEK_sip_plot[i_sim,i_inst,i_plot,0:nr_SIPs],nr_SIPs,n10_plot,r10_plot,min10_plot)

                    m1D_ins_mean = m1D_ins_mean + mEK_bin_plot * nEK_bin_plot
                    g1D_ins_mean = g1D_ins_mean + nEK_bin_plot
            # adding up all realisations is finished

            m1D_ins_mean = m1D_ins_mean / g1D_ins_mean
                # weighted mean: wuerde man bei m1d analog zu g1D verfahren gaebe es Probleme, weil Bins, die nur von manchen Instanzen belegt sind, keine vernuenftige mittlere Masse liefern wuerden.
                # bei unbelegten Bins wird nun durch 0 dividiert -> praktischerweise sorgen die nan-values in m1D fuer einen GV-Plot mit Luecken, so wie gewuenscht.
            g1D_ins_mean = g1D_ins_mean / (nr_inst * V[i_sim])  # units kg / m^3 or kg / m^2
            g1D_ins_mean=3*g1D_ins_mean*(m1D_ins_mean**2)*skal_g # convert to units g / cm^3 or g / cm^2
            m1D_ins_mean=FK.m2r(m1D_ins_mean,const_mass2rad)*1e6 # convert radius in um
            if (i_sim >= len(label)):
                label_tmp = None
            else:
                label_tmp = label[i_sim]
            plt.plot(m1D_ins_mean, g1D_ins_mean,markevery=5,label=label_tmp)
            #print(np.nanmin(m1D_ins_mean),np.nanmax(m1D_ins_mean),np.nanmin(g1D_ins_mean),np.nanmax(g1D_ins_mean))
        #<<<<<< loop over i_sim<<<<<<<<<<

    #GCCendif /* (DISCRETE == 0) */

    #GCCif (DISCRETE >= 1)
                if (nr_SIPs > 0):
                    nEK_bin_plot = FK.CountSIPs(nEK_sip_plot[i_sim,i_inst,i_plot,0:nr_SIPs],mEK_sip_plot[i_sim,i_inst,i_plot,0:nr_SIPs],nr_SIPs,nplot)
                g1D_ins_mean = g1D_ins_mean + nEK_bin_plot
            g1D_ins_mean = g1D_ins_mean * dVi / nr_inst # average concentration in 1/m^3
            #print('Summe N', sum(g1D_ins_mean)*dV)
            #print('n1D_ins_mean : ', g1D_ins_mean)
            g1D_ins_mean = g1D_ins_mean  * (np.arange(nplot)+1)*skal_m[i_sim]* 1e3* 1e-6 # mean droplet mass concentration in g/cm^3
            #print('Summe M', sum(g1D_ins_mean)/skal_m*1e-3)
            m1D_ins_mean=(np.arange(nplot)+1)**(1./3.)*FK.m2r(skal_m[i_sim],const_mass2rad)*1e6 # radius in um
            #m1D_ins_mean=np.arange(nplot)+1
            indexlist=g1D_ins_mean.nonzero()
            #print('g1D_ins_mean: ', g1D_ins_mean[indexlist])
            #print('m1D_ins_mean: ', m1D_ins_mean[indexlist])
            plt.plot(m1D_ins_mean[indexlist], g1D_ins_mean[indexlist],'o-',label=label[i_sim])
    #GCCendif /* (DISCRETE >= 1) */

    #GCCif (IREF == 9 || IREF == 1)
    if (outfallingGV == 0):
        custom_cycler=Timecycler_dict[Time_cycler]
        #print(custom_cycler)
        #plt.rc('axes', prop_cycle=custom_cycler)
        plt.gca().set_prop_cycle(custom_cycler)
        for i_plot in iGV_select:
            #GCCif (IREF == 9)
            plt.plot(rgridREF*1e6, SizeDistrREF[i_plot,:]*1e3,color='k',marker='',label='BIN')
            #print(min(rgridREF*1e6),max(rgridREF*1e6),SizeDistrREF[i_plot,:].min(),SizeDistrREF[i_plot,:].max())
            #GCCendif /* (IREF == 9) */
            #GCCif (IREF == 1)
            indREF= int(i_plot/3)
            #print(indREF,i_plot)
            plt.plot(rgridREF[indREF,:nbinREF[indREF]]*1e3,SizeDistrREF[indREF,:nbinREF[indREF]],color='k',marker='',label='BIN')
            #GCCendif /* (IREF == 1) */
    #GCCendif  /* (IREF == 9 || IREF == 1) */

    #GCCif (IREF == 2)
    if (outfallingGV == 0):
        REF.PlotSD_0D()
    #GCCendif /* (IREF == 2) */

    if (title is not None): plt.title(title)
    if (label[0] is not None):
        if (iaddsimREF == 0) or (outfallingGV >= 1):
            plt.legend(fontsize=8,ncol=leg_ncol)
        else:
            if (leg_ncol == 1):
                #plt.legend(label[:nr_sims+iaddsimREF],fontsize=10,ncol=leg_ncol)
                handles, labels = plt.gca().get_legend_handles_labels()
                print('labels: ', labels)
                handles = handles[:nr_sims] + [handles[nr_sims*nrGV_select]]
                labels  = labels[:nr_sims]  + [labels[nr_sims*nrGV_select]]
                plt.legend(handles,labels,fontsize=8,loc='upper center')
            else:
                handles, labels = plt.gca().get_legend_handles_labels()
                handles = handles[:nr_sims] + [handles[nr_sims*nrGV_select]]
                labels  = labels[:nr_sims]  + [labels[nr_sims*nrGV_select]]
                plt.legend(flip(handles,leg_ncol),flip(labels,leg_ncol),fontsize=8,ncol=leg_ncol,loc='upper center')

    if iTimes_select is None:
        plt.savefig(fp_out+'SD'+suffix+'.png', format='png', bbox_inches='tight')
        plt.savefig(fp_out+'SD'+suffix+'.pdf', format='pdf', bbox_inches='tight')
    else:
        list_of_str = ["{0:04}".format(i) for i in iTimes_select]
        #print('list_of_str',list_of_str)
        plt.savefig(fp_out+'SD'+suffix+ '_'.join(list_of_str)+'.png', format='png', bbox_inches='tight')
        plt.savefig(fp_out+'SD'+suffix+ '_'.join(list_of_str)+'.pdf', format='pdf', bbox_inches='tight')

    plt.close(1)
# end function PlotGV  -----------------------------------------------------------------------------------

def PlotRelStdDev(mEK_sip_plot, t_vec_GVplot,iMultipleSims=0,fp_out='',nr_inst_vec=None,label=None,title=None,ilastGVonly=None):
# Generate RStD-Plot.
# Optinoally, additional plots are generated:
# 1. histogram of mass of largest SIP
# 2. mean masses of largest and second largest SIP

    ihistogram = 0 # 0 oder 1
    imaxelement = 1 # 0 oder 1

    if (iMultipleSims == 0):
        mEK_sip_plot = np.expand_dims(mEK_sip_plot, axis=0) # add a dummy first dimension with length 1
    [nr_sims,nr_inst,nr_GVplot,nr_SIPs]= mEK_sip_plot.shape

    if (nr_inst_vec is None): nr_inst_vec=np.array([nr_inst])
    print('Plotte RStD: ',nr_sims, nr_inst,nr_GVplot,nr_SIPs)
    #print(plt.rcParams)
    #ls_vec=['-', ':','--', '-.']
    #cycl_col_reg = cycler(color=list('bgrcmy'))  # regular color definitions
    #cycl_col_1   = cycler(color=[plt.cm.cool(i) for i in np.linspace(0, 1, 7)])
    #cycl_ls = cycler(linestyle=ls_vec[0:nr_sims])
    #custom_cycler = ( cycl_col_reg*cycl_ls )

    Simcycler_dict  = cycler_dict(nr_sims)
    print(Simcycler_dict)
    plt.rc('axes', prop_cycle=Simcycler_dict[Sim_cycler])

    fig0 = plt.figure(figsize=(4,6), dpi=300)
    ax_vec = fig0.subplots(1+ihistogram+imaxelement,1)
    if ((ihistogram == 0) and  (imaxelement == 0)):
        ax_RStD=ax_vec
    else:
        ax_RStD=ax_vec[0]

    fig0.subplots_adjust(left=0.22, bottom=0.15, right=0.85, top=0.94, wspace=0.3, hspace=0.3)

    plt.xlabel('Time / s')
    plt.xscale("linear")
    plt.yscale("linear")
    #plt.xlim((1e-1,2e3))
    #plt.ylim((1e-4,1e1))
    #plt.xlim((0,2600))
    #plt.ylim((0,1))
    #if (xticks is not None): plt.xticks(xticks)
    plt.rc('lines.markersize')

    for i_sim in range(0,nr_sims):
        print('Plotte RStD von Sims ',  i_sim)
        nr_inst = nr_inst_vec[i_sim]
        maxM=np.amax(mEK_sip_plot[i_sim,0:nr_inst,:,:],axis=2)
        print('maxM',maxM.shape)
        StD=np.std(maxM,axis=0)
        MeanV=np.mean(maxM,axis=0)
        RStD=StD/MeanV
        print('RStD',RStD.shape)

        if (ihistogram == 1):
        #plot coloured histogram
            plot_ind=1
            ax_vec[plot_ind].plot(t_vec_GVplot,MeanV)
            ax_vec[plot_ind].plot(t_vec_GVplot,MeanV-StD,'')
            ax_vec[plot_ind].plot(t_vec_GVplot,MeanV+StD,'')
            for i_plot in  range(0,nr_GVplot):
                binedges=0.5+np.arange(41)
                bincenters=1+np.arange(40)
                hist1D,dummy=np.histogram(maxM[0:nr_inst,i_plot],bins=binedges,density=True)
                print(hist1D)
                ax_vec[plot_ind].scatter(np.zeros(40)+t_vec_GVplot[i_plot],bincenters,5, c= 1+(np.log(hist1D)/50.))

        if (imaxelement == 1):
        # plot time series of largest and second largest element
            plot_ind = 1+ ihistogram
            maxM_save = maxM
            maxM_mean, max2M_mean=FK.get_first_second_highest_SIPmass(mEK_sip_plot[i_sim,0:nr_inst,:,:]) # returns array of shape (nr_inst,nr_GVplot)
            print('M1a',np.mean(maxM,axis=0))
            print('M1b',maxM_mean)
            print('M2' ,max2M_mean)
            line, = ax_vec[plot_ind].plot(t_vec_GVplot,maxM_mean,'')
            col_tmp = line.get_c()
            #ax_vec[plot_ind].set_xlim([0,2500])
            ax_vec[plot_ind].plot(t_vec_GVplot,max2M_mean,'.',color=col_tmp)
            ax_vec[plot_ind].set_ylabel("$<m_{max}>, <m_{2max}>$")

        # create plot of Relative Standard Deviation
        ax_RStD.plot(t_vec_GVplot,RStD)
        #ax_RStD.set_xlim([0,2500])
        ax_RStD.set_ylabel("$\sigma(m_{max})/<m_{max}>$")

#!GCCif (IREF == 2)
    REF.PlotRStD_0D(ax_RStD)
#!GCCendif /* (IREF == 2) */

    plt.savefig(fp_out+'Hist_max_mass.png', format='png', bbox_inches='tight')
    plt.savefig(fp_out+'Hist_max_mass.pdf', format='pdf', bbox_inches='tight')
    plt.close(1)
# end function PlotRelStdDev  -----------------------------------------------------------------------------------

def PlotMomentsProf(MOMsave, t_vec_MOMsave, fp_out='', iMultipleSims=0, iMean=0, iplot_onlyMOMmean=0,
                    nr_inst_vec=None, nz_vec=None,  nt_vec=None, dz_vec = None,
                    label=[None], title=None,text=None,
                    iTimes_select=None,
                    iDmean=0, iplotStdDev=0, MOMprof_StdDev=None,
                    dz_plot = None,
                    iSIPplot=0, nr_SIPs_prof_mean=None, t_vec_GVplot=None, ntGV_vec=None
                    ):
    #print('iTimes_select', iTimes_select)
    #Plot profiles of moments at given times
    dz_plot = 100
    iyscal_km = 1  # =0 yscale uses unit m, =1 yscale uses unit km
   #MOMsave contains moment data, separately for each realisation (iMean=0)
    #MOMsave contains moment data, average over all realisations (iMean=1)
    #iplot_onlyMOMmean: 0 plot each realisation separately; =1 plot only ensemble average
    #dz is only defined "params.txt", not in "params_plot.txt" => if ICOMPPLOT = 2, then dz must be provided as an argument
    iaddsim1DREF=0
#GCCif (IREF == 9)
    ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF,rgridREF,mgridREF=REF.get_RefMetaData()
    print('REF meta data', ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF)
    if (TsimREF == 3601) or (TsimREF == 7201):
        TsimREF = TsimREF - 1
    t_vecREF=TsimREF/(ntREF-1)*np.arange(ntREF)
    IndexTimeDictREF={t_vecREF[i]: i for i in range(ntREF)}
    if (iyscal_km == 0):
        z_vecREF=dzREF*np.arange(nzREF)
    else:
        z_vecREF=dzREF*np.arange(nzREF)*1e-3
    #print(t_vecREF)
    MOMsaveREF=REF.get_RefProfileData(nzREF,ntREF)  #function output : np.zeros([4,nt,nz])
    iaddsim1DREF = 1
    if (iDmean == 1):
        DmeanREF = np.zeros([ntREF,nzREF])
        for iREF in range(ntREF):
            a = MOMsaveREF[1,iREF,:]
            b = MOMsaveREF[0,iREF,:]
            c = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
            DmeanREF[iREF,:] = FK.m2r(c,const_mass2rad)*2e6 # diameter in um
        DmeanREF[DmeanREF == 0] = np.inf
#GCCendif /* (IREF == 9) */

    if (iMultipleSims == 0):
        MOMsave = np.expand_dims(MOMsave, axis=0) # add a dummy first dimension with length 1
        nt_vec = np.array([t_vec_MOMsave.size],dtype='int')
        t_vec_MOMsave = np.expand_dims(t_vec_MOMsave, axis=0) # add a dummy first dimension with length 1

    if (iMean == 1):
        MOMmean=MOMsave.copy()
        iplot_onlyMOMmean = 1
        MOMsave = np.expand_dims(MOMsave, axis=1) # fuege als 1.Dimension Dimension mit Laenge 1 ein

    if (iMean == 0):
        MOMmean=FK.CIO_MOMmean(data_in=MOMsave,ikeep1D=1)

    [nr_sims,nr_inst,nt_time,nz,nr_Mom]=MOMsave.shape
    print('Plotte Moment-Profile: ',nr_sims, nr_inst,nt_time,nz,nr_Mom)

    if (iplotStdDev > 0 and MOMprof_StdDev is None):
        if (iMean == 1):
            print("if only Mean Values are provided, then StdDev can not be computed and must be provided as well")
            iplotStdDev = 0
        if (iMean == 0):
            if (iplotStdDev == 1):
                print("Compute standard deviation")
                MOMprof_StdDev = np.std(MOMsave, axis=1)
                MOMprof_StdDev = np.expand_dims(MOMprof_StdDev, axis=0)
            if (iplotStdDev == 2):
                print("Compute 10, 90 percentile")
                MOMprof_StdDev = np.abs(np.percentile(MOMsave, [10,90], axis=1)-MOMmean) # given as positive deviations from the mean value

    if (iMultipleSims == 1):
        if (nz_vec      is None): print('Multiple Sims, but no nz_vec      provided. Assume nz_vec      =',nz       ,' for all sims')
        if (nt_vec      is None): print('Multiple Sims, but no nt_vec      provided. Assume nt_vec      =',nt_time  ,' for all sims')
        if (iMean != 1): print('provide mean profiles, if iMultipleSims = 1')

    #GCCif (ICOMPPLOT  == 2)
    if (dz_vec is None): sys.exit('provide dz_vec!')
    if (nz_vec is None): sys.exit('provide nz_vec!')
    #if (iMultipleSims == 0):
        #dz_vec = np.array([dz_vec])
        ##nz_vec = np.array([nz])
    #GCCendif /* (ICOMPPLOT  == 2) */
    #GCCif (ICOMPPLOT  == 1)
    dz_vec = np.array([dz])
    nz_vec = np.array([nz])
    #GCCendif /* (ICOMPPLOT  == 1) */

    #determine the indices of the selected times
    if iTimes_select is None:
        iTimes_select=[int(i) for i in np.arange(0,t_vec_MOMsave.max()+1,600)]
    iTimes_select=FK.as_list(iTimes_select) # make sure that iTimes_select is a list, which is iterable
    nrTimes_select=len(iTimes_select)
    iind_select=np.zeros([nr_sims,nrTimes_select],dtype='int')

    for i_sim in range(0,nr_sims):
        #create temporary dictionary to find indices
        IndexTimeDict={t_vec_MOMsave[i_sim,i]: i for i in range(nt_vec[i_sim])}
        for ii,Time in enumerate(iTimes_select):
            iind_select[i_sim,ii]=IndexTimeDict.get(Time,-1)

    labelwithBin=label
    #GCCif (IREF == 9)
    iind_REF=[IndexTimeDictREF.get(Time) for Time in iTimes_select ]
    labelwithBin=label+['BIN']
    print('labelwithBin',labelwithBin)
    print('iTimes_select: ', iTimes_select)
    print(IndexTimeDictREF, 'BBB', IndexTimeDict)
    print(iind_REF,'AAA',iind_select)
    #GCCendif /* (IREF == 9) */

    #GCCif (addSIPplot > 0)
    if (iSIPplot == 1):
        if (iMultipleSims == 0):
            nr_SIPs_prof_mean = np.expand_dims(nr_SIPs_prof_mean, axis=0) # add a dummy first dimension with length 1
            ntGV_vec = np.array([t_vec_GVplot.size],dtype='int')
            t_vec_GVplot = np.expand_dims(t_vec_GVplot, axis=0) # add a dummy first dimension with length 1

        iindGV_select=np.zeros([nr_sims,nrTimes_select],dtype='int')
        for i_sim in range(0,nr_sims):
            #create temporary dictionary to find indices
            IndexTimeGVDict={t_vec_GVplot[i_sim,i]: i for i in range(ntGV_vec[i_sim])}
            for ii,Time in enumerate(iTimes_select):
                iindGV_select[i_sim,ii]=IndexTimeGVDict.get(Time,-1)

        [nrSIP_sims,ntSIP_time,nzSIP]=nr_SIPs_prof_mean.shape
    #GCCendif /* (addSIPplot > 0) */


#define plot styles
    Timecycler_dict = cycler_dict(nrTimes_select)
    Simcycler_dict  = cycler_dict(nr_sims+iaddsim1DREF)
    Simcycler_dictwithoutREF  = cycler_dict(nr_sims)
    #choose plot style for specific plot
    if ((nr_sims+iaddsim1DREF) == 1):
        custom_cycler = ( Timecycler_dict[Time_cycler] )
    else:
        if (nrTimes_select == 1):
            custom_cycler = Simcycler_dict[Sim_cycler]
            custom_cyclerwithoutREF = Simcycler_dictwithoutREF[Sim_cycler]
        else:
            custom_cycler = ( Timecycler_dict[Time_cycler] * Simcycler_dict[Sim_cycler] )
            custom_cyclerwithoutREF = ( Timecycler_dict[Time_cycler] * Simcycler_dictwithoutREF[Sim_cycler] )

    xlabel_vec=['$\lambda_0$ / m$^{-3}$',
                '$\lambda_1$ / (kg m$^{-3}$)',
                '$\lambda_2$ / (kg$^{2}$ m$^{-3}$)',
                '$\lambda_3$ / (kg$^{3}$ m$^{-3}$)',
                '$D_\mathrm{mean}$ / $\mu$m',
                '$N_\mathrm{SIPs}$ / 100m']
    #GCCif (addSIPplot == 2)
    xlabel_vec[5]='$N_\mathrm{SIP,GB}$'
    #GCCendif /* (addSIPplot == 2) */
    print('iMom_select_const before',iMom_select_const)
    if iMom_select_const is None:
        iMom_select=list(range(nr_Mom+iDmean+iSIPplot))
        if (iDmean == 0) and (iSIPplot == 1):
            iMom_select[nr_Mom]=nr_Mom+1
    else:
        iMom_select = iMom_select_const.copy()
        print('else branch: ', iMom_select, iMom_select_const)
        if (iDmean == 1): iMom_select.append(4)
        if (iSIPplot == 1): iMom_select.append(5)

    print('shape MOMmean', MOMmean.shape)
    print('iMom_select',iMom_select)
    print('iTimes_select',iTimes_select)
    print('iplot_onlyMOMmean',iplot_onlyMOMmean)

    if (dz_plot is not None):
        print('BEFORE nz_vec:', nz_vec)
        print('BEFORE dz_vec:', dz_vec)

        # Simulations with dz smaller than dz_plot are re-gridded with a resolution of dz_plot
        for i_sim in range(0,nr_sims):
            if (dz_vec[i_sim] < dz_plot):
                stride = dz_plot/dz_vec[i_sim]
                if (stride != int(stride)):
                    print('dz_plot is not a multiple of dz_vec', dz_plot, dz_vec)
                stride = int(stride)
                nz_new = nz_vec[i_sim]/stride
                if (nz_new != int(nz_new)):
                    print('Lz is not a multiple of dz_plot', nz_new, dz_vec)
                nz_new = int(nz_new)
                print('Regridding: i_sim, stride, nz, nz_new = ', i_sim, stride, nz_vec[i_sim], nz_new)
                tmp = np.mean(MOMsave[i_sim,:,:,:nz_vec[i_sim],:].reshape(nr_inst,nt_time,nz_new,stride,nr_Mom),axis=3)
                MOMsave[i_sim,:,:,:nz_new,:] = tmp
                tmp = np.mean(MOMmean[i_sim,:,:nz_vec[i_sim],:].reshape(nt_time,nz_new,stride,nr_Mom),axis=2)
                MOMmean[i_sim,:,:nz_new,:] = tmp

                if (iplotStdDev == 1):
                    tmp = np.mean(MOMprof_StdDev[:,i_sim,:,:nz_vec[i_sim],:].reshape(1,nt_time,nz_new,stride,nr_Mom),axis=3)
                    MOMprof_StdDev[:,i_sim,:,:nz_new,:]= tmp
                if (iplotStdDev == 2):
                    tmp = np.mean(MOMprof_StdDev[:,i_sim,:,:nz_vec[i_sim],:].reshape(2,nt_time,nz_new,stride,nr_Mom),axis=3)
                    MOMprof_StdDev[:,i_sim,:,:nz_new,:]= tmp
                #GCCif (addSIPplot > 0)
                if (iSIPplot == 1):
                    tmp = np.mean(nr_SIPs_prof_mean[i_sim,:,:nz_vec[i_sim]].reshape(ntSIP_time,nz_new,stride),axis=2)
                    nr_SIPs_prof_mean[i_sim,:,:nz_new] = tmp
                #GCCendif /* (addSIPplot > 0) */

                dz_vec[i_sim] = dz_plot
                nz_vec[i_sim] = nz_new
        print('NEW nz_vec:', nz_vec)
        print('NEW dz_vec:', dz_vec)

    #GCCif (Prof_Panel == 1)
    nr_panels = len(iMom_select)
    fig0 = plt.figure(figsize=(3 + 2 * (nr_panels-1),3), dpi=300)
    print('handlelength', plt.rcParams["legend.handlelength"])
    ax_vec = fig0.subplots(ncols=nr_panels,squeeze=False)
    #GCCendif /* (Prof_Panel == 1) */

    for iiMom,iMom in enumerate(iMom_select):

        #GCCif (Prof_Panel == 0)
        fig0 = plt.figure(figsize=(4,6), dpi=300)
        ax_vec = fig0.subplots(ncols=1,squeeze=False)
        print('a: ',ax_vec.shape)
        ip = 0
        #GCCelse
        ip = iiMom
        #GCCendif /* (Prof_Panel == 1) */

        print('iMom,iiMom,ip: ', iMom, iiMom, ip)
        #ax = fig0.add_subplot(111)
        #fig0.subplots_adjust(left=0.22, bottom=0.15, right=0.85, top=0.94, wspace=0.3, hspace=0.3)
        #if (title is not None): plt.suptitle(title)
        #ax.set_cmap('Reds')

        if (iMom != 5):
            #plt.rc('axes', prop_cycle=custom_cycler)
            ax_vec[0,ip].set_prop_cycle(custom_cycler)
        else:
            #plt.rc('axes', prop_cycle=custom_cyclerwithoutREF)
            ax_vec[0,ip].set_prop_cycle(custom_cyclerwithoutREF)

        ax_vec[0,ip].set_xlabel(xlabel_vec[iMom])

        if (iMom <= 3):
            prof_data = MOMsave[:,:,:,:,iMom]
            prof_data_mean = MOMmean[:,:,:,iMom]
            if (iplotStdDev > 0):
                prof_data_STD = MOMprof_StdDev[:,:,:,:,iMom]
        if (iMom == 4):
            a = MOMsave[:,:,:,:,1]
            b = MOMsave[:,:,:,:,0]
            c = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
            prof_data = FK.m2r(c,const_mass2rad)*2e6 # diameter in um
            a = MOMmean[:,:,:,1]
            b = MOMmean[:,:,:,0]
            c = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
            prof_data_mean = FK.m2r(c,const_mass2rad)*2e6 # diameter in um
            if (iplotStdDev > 0):
                prof_data_STD = np.zeros([1,nr_sims,nt_time,nz])
        if (iMom == 5):
            prof_data_mean = nr_SIPs_prof_mean
            prof_data = None

        #get x-axis range over all times
        xMax=0
        for ii,iTimes in enumerate(iTimes_select):
            for i_sim in range(0,nr_sims):
                if (iMom != 5):
                    iT=iind_select[i_sim,ii]
                    if (iT != -1):
                        xMax=max(prof_data[i_sim,:,iT,:nz_vec[i_sim]].max(),xMax)
                else:
                    iT=iindGV_select[i_sim,ii]
                    if (iT != -1):
                        xMax=max(prof_data_mean[i_sim,iT,:nz_vec[i_sim]].max(),xMax)

        print('xMax: ',xMax)

        if (iHalfDom == 1):
            print('use halfDomLinDec plot options')
            ilog_vec=[1,1,1,1,0,1]
            ilog=ilog_vec[iMom]
            if (ilog==1):
                ax_vec[0,ip].set_xscale("log")
                #skal = np.array([1e-8,1e-2,1e-8,1e-10,1e-2,1e4])
                skal = np.array([1e-8,1e-3,1.1e-8,1e-10,1e-2,1e4])
                #xMin_act = xMax*skal[iMom]
                #xMax_act = xMax*5
                if (iMom == 0):
                    xMax_act = 3e9
                    xMin_act = xMax_act * skal[iMom]
                    ax_vec[0,ip].set_xlim([xMin_act,xMax_act])
                    locmin = matplotlib.ticker.LogLocator(base=10.0,numticks=12)
                    ax_vec[0,ip].xaxis.set_minor_locator(locmin)
                    ax_vec[0,ip].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
                if (iMom == 1):
                    xMax_act = 1.01e-2
                    xMin_act = 1e-5 #xMax_act * skal[iMom]
                    #ax_vec[0,ip].minorticks_on()
                    locmin = matplotlib.ticker.LogLocator(base=10.0,numticks=12)
                    ax_vec[0,ip].xaxis.set_minor_locator(locmin)
                    ax_vec[0,ip].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
                if (iMom == 2):
                    xMax_act = 1e-7
                    xMin_act = xMax_act * skal[iMom]
                    locmin = matplotlib.ticker.LogLocator(base=10.0,numticks=12)
                    ax_vec[0,ip].xaxis.set_minor_locator(locmin)
                    ax_vec[0,ip].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
                if (iMom == 5):
                    xMin_act = 1
                    xMax_act = 1e3
                    locmin = matplotlib.ticker.LogLocator(base=10.0,numticks=12)
                    ax_vec[0,ip].xaxis.set_minor_locator(locmin)
                    ax_vec[0,ip].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
                print('iMom,ip, xMin_act, xMax_act: ', iMom,ip, xMin_act, xMax_act)
                ax_vec[0,ip].set_xlim([xMin_act,xMax_act])
            else:
                ax_vec[0,ip].set_yscale("linear")
                if (iMom == 5):
                    xMax = 200
                if (iMom == 4):
                    xMax = 800
                ax_vec[0,ip].set_xlim([0,xMax])

        if (iEmptydom == 1):
            print('use Emptydom plot options')
            ilog_vec=[1,1,1,1,0,1]
            ilog=ilog_vec[iMom]
            if (ilog==1):
                ax_vec[0,ip].set_xscale("log")
                #skal = np.array([1e-8,1e-2,1e-8,1e-10,1e-2,1e4])
                skal = np.array([1e-7,1e-3,1.0e-5,1e-10,1e-2,1e4])
                xMin_act = xMax*skal[iMom]
                xMax_act = xMax*5
                if (iMom == 0):
                    xMax_act = 1e8
                    xMin_act = xMax_act * skal[iMom]
                    ax_vec[0,ip].set_xlim([xMin_act,xMax_act])
                    locmin = matplotlib.ticker.LogLocator(base=10.0,numticks=12)
                    ax_vec[0,ip].xaxis.set_minor_locator(locmin)
                    ax_vec[0,ip].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
                if (iMom == 2):
                    xMax_act = 5e-11
                    xMin_act = xMax_act * skal[iMom]
                    locmin = matplotlib.ticker.LogLocator(base=10.0,numticks=12)
                    ax_vec[0,ip].xaxis.set_minor_locator(locmin)
                    ax_vec[0,ip].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
                if (iMom == 5):
                    xMin_act = 1e-1
                    xMax_act = 1e2

                    locmin = matplotlib.ticker.LogLocator(base=10.0,numticks=12)
                    ax_vec[0,ip].xaxis.set_minor_locator(locmin)
                    ax_vec[0,ip].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
                print('iMom,ip, xMin_act, xMax_act: ', iMom,ip, xMin_act, xMax_act)
                ax_vec[0,ip].set_xlim([xMin_act,xMax_act])
            else:
                ax_vec[0,ip].set_yscale("linear")
                if (iMom == 4):
                    ax_vec[0,ip].xaxis.set_ticks([0,200,400,600])
                    ax_vec[0,ip].set_xlim([0,740])

        if ((1-iEmptydom)*(1-iHalfDom) == 1):
            ilog_vec=[1,1,1,1,0,1]
            ilog=ilog_vec[iMom]
            if (ilog==1):
                ax_vec[0,ip].set_xscale("log")
                #skal = np.array([1e-8,1e-2,1e-8,1e-10,1e-2,1e4])
                skal = np.array([1e-8,1e-3,1e-6,1e-10,1e-2,1e4])
                xMin_act = xMax*skal[iMom]
                xMax_act = xMax*5
                if (iMom == 5):
                    xMin_act = .1
                    xMax_act = 1e2
                    #xMin_act = 1
                    #xMax_act = 1e3
                ax_vec[0,ip].set_xlim([xMin_act,xMax_act])
            else:
                ax_vec[0,ip].set_yscale("linear")
                if (iMom == 5):
                    xMax = 200
                if (iMom == 4):
                    xMax = 800
                ax_vec[0,ip].set_xlim([0,xMax])


        if (ip > 0):  # include y-axis only in left most plot
            print('suppress y ticks')
            ax_vec[0,ip].get_yaxis().set_ticklabels([]) # supress y-axis tick labels in all but the left most plot
        else:
            if (iyscal_km == 0):
                ax_vec[0,ip].set_ylabel('$z$ / m')
            else:
                ax_vec[0,ip].set_ylabel('$z$ / km')


        if (iMom == 4):
            #empty grid boxes have no well defind mean radius, data points are left out in profile plot
            print('max(Dmean): ', np.nanmax(prof_data))
            prof_data_mean[prof_data_mean == 0] = np.inf
            prof_data[prof_data == 0] = np.inf
            if (iEmptydom == 1):
                print('remove spikes in Dmean')
                prof_data_mean[prof_data_mean > 550] = np.inf
                prof_data[prof_data > 550] = np.inf

        #plt.ylim([700,1000])
        for ii,iTimes in enumerate(iTimes_select):
            for i_sim in range(0,nr_sims):
                if (iyscal_km == 0):
                    z_vec=(np.arange(nz_vec[i_sim])+0.5)*dz_vec[i_sim]
                else:
                    z_vec=(np.arange(nz_vec[i_sim])+0.5)*dz_vec[i_sim]*1e-3
                iT=iind_select[i_sim,ii]
                if (iMom == 5):
                    iT=iindGV_select[i_sim,ii]
                if ((iplot_onlyMOMmean == 0) and (iMom !=5)):
                    for k in range(0,nr_inst):
                        ax_vec[0,ip].plot(prof_data[i_sim,k,iT,:nz_vec[i_sim]],z_vec,'r:')
                    if (iplotStdDev == 0):
                        ax_vec[0,ip].plot(prof_data_mean[i_sim,iT,:nz_vec[i_sim]],z_vec,'k:')
                    else:
                        ax_vec[0,ip].errorbar(prof_data_mean[i_sim,iT,:nz_vec[i_sim]], z_vec, xerr=np.squeeze(prof_data_STD[:,i_sim,iT,:nz_vec[i_sim]]), fmt='k:',linewidth=2.0,capsize=4)
                else:
                    if (iplotStdDev == 0):
                        ax_vec[0,ip].plot(prof_data_mean[i_sim,iT,:nz_vec[i_sim]],z_vec)
                    else:
                        ax_vec[0,ip].errorbar(prof_data_mean[i_sim,iT,:nz_vec[i_sim]],z_vec, xerr=np.squeeze(prof_data_STD[:,i_sim,iT,:nz_vec[i_sim]]),capsize=4)

            #plt.savefig(fp_out+'ProfMom'+ str(iMom)+ '_t' + str(iTimes) + '.png', format='png', bbox_inches='tight',dpi=300)
            #plt.close()


#GCCif (IREF == 9)
            if (iMom <= 3):
                print('ip, iMom, ii, iind_REF[ii]', ip, iMom, ii, iind_REF[ii])
                ax_vec[0,ip].plot(MOMsaveREF[iMom,iind_REF[ii],:],z_vecREF,color='k')  #,col_list[ii]+':'
            if (iMom == 4):
                ax_vec[0,ip].plot(DmeanREF[iind_REF[ii],:],z_vecREF,color='k')  #,col_list[ii]+':'
#GCCendif /* (IREF == 9) */

        if (iHalfDom == 1):
            if (iMom == 0):
                #for label in ax.get_xticklabels()[::2]:
                    #label.set_visible(False)

                xticklabels = ax_vec[0,ip].get_xticklabels()
                print('before 2: ', [item.get_text() for item in xticklabels])
                #xticklabels[1::2] = ['']*len(xticklabels[1::2])
                #print('after ', xticklabels.get_text())
                #ax_vec[0,ip].set_xticklabels(xticklabels)

        #GCCif (Prof_Panel == 0)
        if (ip == 0):
        #GCCendif /* (Prof_Panel == 0) */
        #GCCif (Prof_Panel == 1)
        ip_legend_select = labelMoments_dict.get('ipc',0)
        if (ip == ip_legend_select):
        #GCCendif /* (Prof_Panel == 1) */

            if (label[0] is not None):
                if (iExtraLegend == 4):
                    lc = 'upper right'
                    fs = labelMoments_dict.get('fs',6)
                    leg_ncol = labelMoments_dict.get('ncol',1)
                    ipr_leg = 0
                    ipc_leg = 3
                    dummy_lines = []
                    for b_idx, b in enumerate(['-', ':','--', '-.',(0, (3, 10, 1, 10))]):
                        dummy_lines.append(ax_vec[0,ipc_leg].plot([],[], c="black", ls = b)[0])
                    legend2 = ax_vec[0,ipc_leg].legend([dummy_lines[i] for i in range(5)], ["0", "10", "20", "30", "60"], title="time / min", loc='upper left',fontsize=fs,handlelength = 5)
                    print('legend2 created', ipr_leg, ipc_leg)
                    #ax_vec[ipr_leg,ipc_leg].add_artist(legend2) this call is not necessary as the additional legend and the original legend are added in different panels.

                print('legend iMom,iiMom', iMom, iiMom)
                lc=labelMoments_dict.get('loc','lower right')
                fs=labelMoments_dict.get('fs',7)
                leg_ncol = labelMoments_dict.get('ncol',1)
                if (iMom != 5):
                    ax_vec[0,ip].legend(labelwithBin[:nr_sims+iaddsim1DREF],fontsize=fs,loc=lc, ncol=leg_ncol)
                else:
                    ax_vec[0,ip].legend(label[:nr_sims],fontsize=fs,loc=lc, ncol=leg_ncol)
            if text is not None:
                #print('len(text),type(text)',len(text),type(text))
                if type(text) is str:
                    ax_vec[0,ip].text(0.96, 0.92, text, horizontalalignment='right',verticalalignment='top', transform=ax_vec[0,ip].transAxes,bbox=dict(facecolor='none',pad=2),fontsize=6)
                else:
                    ha='right'
                    if (len(text) == 4):
                        ha=text[3]
                    ax_vec[0,ip].text(text[1],text[2], text[0], horizontalalignment=ha,verticalalignment='top', transform=ax_vec[0,ip].transAxes,bbox=dict(facecolor='none',pad=2),fontsize=6)

        if (title is not None): plt.suptitle(title)
        #GCCif (Prof_Panel == 0)
        if iTimes_select is None:
            plt.savefig(fp_out+'ProfMom'+ str(iMom) + '.png', format='png', bbox_inches='tight',dpi=600)
            plt.savefig(fp_out+'ProfMom'+ str(iMom) + '.pdf', format='pdf', bbox_inches='tight',dpi=600)
            plt.close()
        else:
            list_of_str = ["{0:04}".format(i) for i in iTimes_select]
            #print('list_of_str',list_of_str)
            plt.savefig(fp_out+'ProfMom'+ str(iMom) + '_t' + '_'.join(list_of_str)+'.png', format='png', bbox_inches='tight')
            plt.savefig(fp_out+'ProfMom'+ str(iMom) + '_t' + '_'.join(list_of_str)+'.pdf', format='pdf', bbox_inches='tight')
            plt.close()
        #GCCendif /* (Prof_Panel == 0) */

    #GCCif (Prof_Panel == 1)
    list_of_strMOM = ["{0:01}".format(i) for i in iMom_select]
    if iTimes_select is None:
        plt.savefig(fp_out+'ProfMom'+ '_'.join(list_of_strMOM) + '.png', format='png', bbox_inches='tight',dpi=300)
        plt.savefig(fp_out+'ProfMom'+ '_'.join(list_of_strMOM) + '.pdf', format='pdf', bbox_inches='tight',dpi=300)
        plt.close()
    else:
        list_of_str = ["{0:04}".format(i) for i in iTimes_select]
        #print('list_of_str',list_of_str)
        plt.savefig(fp_out+'ProfMom'+ '_'.join(list_of_strMOM) + '_t' + '_'.join(list_of_str)+'.png', format='png', bbox_inches='tight')
        plt.savefig(fp_out+'ProfMom'+ '_'.join(list_of_strMOM) + '_t' + '_'.join(list_of_str)+'.pdf', format='pdf', bbox_inches='tight')
        plt.close()
    #GCCendif /* (Prof_Panel == 1) */

# >>>>>>>>>>>>>>>>>>>> end function PlotMomentsProf >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#GCCif (RZ_scatter == 1 )
def Scatter_r_z(zEK_sip_plot,mEK_sip_plot,nr_SIPs_plot,t_vec_GVplot,fp_out=''):
    import Sedimentation as SD

    #fp_in="data/VBeard_fein/VBeard.dat"
    #VBeard=np.loadtxt(fp_in)
    #print(VBeard.shape)
    #print(VBeard[:,0].min(),VBeard[:,0].max())
    #print(VBeard[:,1].min(),VBeard[:,1].max())
    #print(VBeard[:,2].min(),VBeard[:,2].max())
    [nr_inst,nr_GVplot,nr_SIPs_max]= zEK_sip_plot.shape
    #iGV_select = range(0,nr_GVplot)
    #iGV_select = range(0,nr_GVplot,3)
    iGV_select = range(1,nr_GVplot)
    for i_plot in iGV_select:
        z_sum  = 0
        zm_sum = 0
        m_sum  = 0
        for i_inst in range(0,nr_inst):
            fig0 = plt.figure(figsize=(4,6), dpi=600)
            plt.set_cmap('Reds')
            plt.xlabel('Radius / '+ r'$\mu$' +'m')
            plt.ylabel('dz / m')
            plt.xlim([1e-1,300])
            plt.ylim([3500,4000])
            #plt.yscale('log')
            plt.xscale('log')
            nr_SIPs=nr_SIPs_plot[i_inst,i_plot]
            r=FK.m2r(mEK_sip_plot[i_inst,i_plot,:nr_SIPs],const_mass2rad)*1e6
            #w=SD.Fallg(r*1e-6)*1e2
            #print('min/max r',r.min(),r.max())
            #plt.xlim([r.min(),r.max()])
            plt.scatter(r,zEK_sip_plot[i_inst,i_plot,:nr_SIPs],s=2)  #,s=4
            z_sum  += np.sum(zEK_sip_plot[i_inst,i_plot,:nr_SIPs])
            zm_sum += np.sum(zEK_sip_plot[i_inst,i_plot,:nr_SIPs]*mEK_sip_plot[i_inst,i_plot,:nr_SIPs])
            m_sum  += np.sum(mEK_sip_plot[i_inst,i_plot,:nr_SIPs])
            print("z_mean: ", np.mean(zEK_sip_plot[i_inst,i_plot,:nr_SIPs]),
                              np.mean(zEK_sip_plot[i_inst,i_plot,:nr_SIPs]*mEK_sip_plot[i_inst,i_plot,:nr_SIPs])/np.mean(mEK_sip_plot[i_inst,i_plot,:nr_SIPs]))
            #plt.scatter(r,w)
            #plt.plot(VBeard[:,0],VBeard[:,2])
            #plt.plot( np.array([2e-3,5.e-2])*1e4 ,(4.709172 , 398.3176),'r')
            plt.savefig(fp_out+'Scatter_r_z_Inst'+ str(i_inst)+ '_t' + str(i_plot) + '.png', format='png', bbox_inches='tight',dpi=300)
        print('z_mean all inst: ', z_sum/np.sum(nr_SIPs_plot[:,i_plot]),zm_sum/m_sum)

        #plt.savefig(fp_out+'Scatter_r_z_t' + str(i_plot) + '.png', format='png', bbox_inches='tight',dpi=300)
#>>>>>>>>>>>>>>>>> end function Scatter_r_z >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#GCCendif  /* (RZ_scatter == 1)  */


#GCCif (FLUX_time == 1)
def PlotFluxesTime(FluxIn,FluxOut,FluxInAcc,FluxOutAcc,t_vec,fp_out=''
                   ,iplot_onlyMean=0,iMean=0,iMultipleSims=0,
                   io_sep=0,nt_vec=None,
                   text=None, label=[None], title=None,
                   iplotStdDev=0, FluxIn_STD = None, FluxOut_STD = None, FluxInAcc_STD = None, FluxOutAcc_STD = None
                   ):


    #io_sep:
    #       = 0 plot inflow and outflow in a single panel
    #       = 1 plot separate panels for outflow and inflow
    iaddsim1DREF=0
#GCCif (IREF == 9)
    ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF,rgridREF,mgridREF=REF.get_RefMetaData()
    print('REF meta data', ntREF,nzREF,dzREF,dtREF,TsimREF,nbinREF,scalREF,dlnrREF)
    t_vecREF=TsimREF/(ntREF-1)*np.arange(ntREF)
    if (itime_legend == 2): t_vecREF = t_vecREF/60
    print(t_vecREF)
    FluxesREF=REF.get_RefFluxData(ntREF,nbinREF)  #function output : np.zeros([4,nt,nz])
    iaddsim1DREF=1
#GCCendif /* (IREF == 9) */

    print('iMom_select_const:',iMom_select_const,'end')
    if iMom_select_const is None:
        iMom_select = [0,1,2,3]
        nrMom_select = 4
    else:
        iMom_select = iMom_select_const.copy()
        nrMom_select= len(iMom_select)
    print('iMom_select:',iMom_select,'end')


    skip_inflow=0
#GCCif (ICOMPPLOT == 1 && INFLUX_TOP == 0)
    if (io_sep==1): skip_inflow=1
#GCCendif /* (ICOMPPLOT == 1 && INFLUX_TOP = 0) */

    if (iMultipleSims == 0):
        FluxIn     = np.expand_dims(FluxIn    , axis=0) # add a dummy first dimension with length 1
        FluxOut    = np.expand_dims(FluxOut   , axis=0) # add a dummy first dimension with length 1
        FluxInAcc  = np.expand_dims(FluxInAcc , axis=0) # add a dummy first dimension with length 1
        FluxOutAcc = np.expand_dims(FluxOutAcc, axis=0) # add a dummy first dimension with length 1
        nt_vec = np.array([t_vec.size],dtype='int')
        t_vec = np.expand_dims(t_vec, axis=0) # add a dummy first dimension with length 1

    if (iMultipleSims == 1):
        if (nt_vec is None): print('provide nt_vec, when iMultipleSims = 1')
    print('FluxIn.shape', FluxIn.shape)

    if (iMean == 0):
        #FluxIn    =np.zeros([nr_sims,nr_inst,nr_GVplot-1,nr_MOMs+1])
        [nr_sims,nr_inst,nr_GVplot,nr_MOMs]=FluxIn.shape+np.array([0,0,1,-1])
        print('Shape Flux data: ', nr_sims,nr_inst,nr_GVplot,nr_MOMs)
        FluxInMean     = np.mean(FluxIn,     axis=1)
        FluxOutMean    = np.mean(FluxOut,    axis=1)
        FluxInAccMean  = np.mean(FluxInAcc,  axis=1)
        FluxOutAccMean = np.mean(FluxOutAcc, axis=1)
        if (iplotStdDev == 1 and FluxIn_STD is None):
            FluxIn_STD     = np.expand_dims(np.std(FluxIn,     axis=1), axis=0)
            FluxOut_STD    = np.expand_dims(np.std(FluxOut,    axis=1), axis=0)
            FluxInAcc_STD  = np.expand_dims(np.std(FluxInAcc,  axis=1), axis=0)
            FluxOutAcc_STD = np.expand_dims(np.std(FluxOutAcc, axis=1), axis=0)
        if (iplotStdDev == 2 and FluxIn_STD is None):
            FluxIn_STD     = np.abs(np.percentile(FluxIn    , [10,90], axis=1)-FluxInMean)
            FluxOut_STD    = np.abs(np.percentile(FluxOut   , [10,90], axis=1)-FluxOutMean)
            FluxInAcc_STD  = np.abs(np.percentile(FluxInAcc , [10,90], axis=1)-FluxInAccMean)
            FluxOutAcc_STD = np.abs(np.percentile(FluxOutAcc, [10,90], axis=1)-FluxOutAccMean)
    else:
        [nr_sims,nr_GVplot,nr_MOMs]=FluxIn.shape+np.array([0,1,-1])
        print('Shape Flux data: ', nr_sims,nr_GVplot,nr_MOMs)
        FluxInMean     = FluxIn
        FluxOutMean    = FluxOut
        FluxInAccMean  = FluxInAcc
        FluxOutAccMean = FluxOutAcc
        if (iplotStdDev > 0 and FluxIn_STD is None):
            print("if only Mean Values are provided, then StdDev/Percentiles can not be computed and must be provided as well")
            iplotStdDev = 0

    if (itime_legend == 2): t_vec = t_vec/60

    Simcycler_dict  = cycler_dict(nr_sims+iaddsim1DREF)
    #choose plot style for specific plot
    if ((nr_sims+iaddsim1DREF) == 1):
        pass
    else:
        plt.rc('axes', prop_cycle=Simcycler_dict[Sim_cycler])

    print('skip_inflow, io_sep', skip_inflow, io_sep)
    #plt.set_cmap('Reds')
    ylabel_vec=['$\lambda_0$','$\lambda_1$','$\lambda_2$','$\lambda_3$','$N_\mathrm{SIP}$']
    if (skip_inflow == 0):
        print('Here nr_MOMs nr_sims', nr_MOMs, nr_sims)
        fig, ax = plt.subplots(nrMom_select, 2, sharex=True,figsize=(8,8*(nrMom_select/5.)**0.8), dpi=600) #,sharey=True
        plt.xlabel(xlabel_Mom[itime_legend])

        #for i in range(nr_MOMs+1):
        for iiMom,iMom in enumerate(iMom_select):
            #ax[i,0].yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
            #ax[i,1].yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
            #GCCif (ICOMPPLOT == 2 || INFLUX_TOP != 0)
            for i_sim in range(nr_sims):
                if (iplotStdDev == 0):
                    ax[iiMom,0].plot(t_vec[i_sim,:nt_vec[i_sim]],FluxInMean[i_sim,:nt_vec[i_sim],iMom])
                    ax[iiMom,1].plot(t_vec[i_sim,:nt_vec[i_sim]],FluxInAccMean[i_sim,:nt_vec[i_sim],iMom])
                else:
                    print('FluxInMean[i_sim,:nt_vec[i_sim],iMom].shape', FluxInMean[i_sim,:nt_vec[i_sim],iMom].shape)
                    print('FluxIn_STD[:,i_sim,:nt_vec[i_sim],iMom].shape', FluxIn_STD[:,i_sim,:nt_vec[i_sim],iMom].shape)
                    ax[iiMom,0].errorbar(t_vec[i_sim,:nt_vec[i_sim]],FluxInMean[i_sim,:nt_vec[i_sim],iMom],yerr=np.squeeze(FluxIn_STD[:,i_sim,:nt_vec[i_sim],iMom]))
                    ax[iiMom,1].errorbar(t_vec[i_sim,:nt_vec[i_sim]],FluxInAccMean[i_sim,:nt_vec[i_sim],iMom],yerr=np.squeeze(FluxInAcc_STD[:,i_sim,:nt_vec[i_sim],iMom]))
                    print('IN INST', iMom, FluxIn_STD[:,i_sim,:nt_vec[i_sim],iMom]/FluxInMean[i_sim,:nt_vec[i_sim],iMom])
                    print('IN ACCU', iMom, FluxInAcc_STD[:,i_sim,:nt_vec[i_sim],iMom]/FluxInAccMean[i_sim,:nt_vec[i_sim],iMom])

                #print('FluxIn', FluxInMean[i_sim,:nt_vec[i_sim],iMom])
                #print('FluxInAcc', FluxInAccMean[i_sim,:nt_vec[i_sim],iMom])
                #GCCif (IREF == 9)
            if (iMom < 4): #skip iMom=4: no SIP number flux in bin modell
                ax[iiMom,0].plot(t_vecREF,FluxesREF[0,iMom,:],'k')
                ax[iiMom,1].plot(t_vecREF,FluxesREF[1,iMom,:],'k')
                #GCCendif /* (IREF == 9) */
            #GCCendif  /*(ICOMPPLOT == 2 || INFLUX_TOP != 0)*/

            if (io_sep == 0):
                #combine inflow and outflow in a single panel
                for i_sim in range(nr_sims):
                    if (iplotStdDev == 0):
                        ax[iiMom,0].plot(t_vec[i_sim,:nt_vec[i_sim]],FluxOutMean[i_sim,:nt_vec[i_sim],iMom])
                        ax[iiMom,1].plot(t_vec[i_sim,:nt_vec[i_sim]],FluxOutAccMean[i_sim,:nt_vec[i_sim],iMom])
                    else:
                        ax[iiMom,0].errorbar(t_vec[i_sim,:nt_vec[i_sim]],FluxOutMean[i_sim,:nt_vec[i_sim],iMom],yerr=np.squeeze(FluxOut_STD[:,i_sim,:nt_vec[i_sim],iMom]))
                        ax[iiMom,1].errorbar(t_vec[i_sim,:nt_vec[i_sim]],FluxOutAccMean[i_sim,:nt_vec[i_sim],iMom],yerr=np.squeeze(FluxOut_STD[:,i_sim,:nt_vec[i_sim],iMom]))
                        print('OUT INST', iMom, FluxOut_STD[:,i_sim,:nt_vec[i_sim],iMom]/FluxOutMean[i_sim,:nt_vec[i_sim],iMom])
                        print('OUT ACCU', iMom, FluxOutAcc_STD[:,i_sim,:nt_vec[i_sim],iMom]/FluxOutAccMean[i_sim,:nt_vec[i_sim],iMom])

                #GCCif (IREF == 9)
                if (iMom < 4): #skip iMom=4: no SIP number flux in bin modell
                    ax[iiMom,0].plot(t_vecREF,FluxesREF[2,iMom,:],'k')
                    ax[iiMom,1].plot(t_vecREF,FluxesREF[3,iMom,:],'k')
                #GCCendif /* (IREF == 9) */
            ax[iiMom,0].set_ylabel(ylabel_vec[iMom])

            if (iMom < 4):
                ax[iiMom,0].ticklabel_format(axis='y', style='sci', scilimits=(0,1))
                ax[iiMom,1].ticklabel_format(axis='y', style='sci', scilimits=(0,1))
            if (iMom == 0):
                ax[iiMom,0].set_title('flux kg$^l$/(m$^2$ s)')
                ax[iiMom,1].set_title('acc. flux kg$^l$/m$^2$')
            if (iMom == 4):
                ax[iiMom,0].set_xlabel(xlabel_Mom[itime_legend])
                ax[iiMom,1].set_xlabel(xlabel_Mom[itime_legend])
            #ax[iMom,0].xaxis.set_tick_position('bottom')
            #ax[i,1].xaxis.set_tick_position('bottom')
            #if (i < 4):
                #ax[i,0].get_xaxis().set_ticklabels([])
                #ax[i,1].get_xaxis().set_ticklabels([])
            #else:
                ##ax[i,0].get_xaxis().set_ticklabels()
                ##ax[i,1].get_xaxis().set_ticklabels()

        iMomLeg = 0
        if (label[0] is not None): ax[iMomLeg,0].legend(label[:nr_sims+iaddsim1DREF],fontsize=6,loc='lower left',ncol=leg_ncol)
        if text is not None:
            print('len(text),type(text)',len(text),type(text))
            if type(text) is str:
                ax[iMomLeg,0].text(0.96, 0.92, text, horizontalalignment='right',verticalalignment='top', transform=ax[iMomLeg,0].transAxes,bbox=dict(facecolor='none',pad=2),fontsize=6)
            else:
                ha='right'
                if (len(text) == 4):
                    ha=text[3]
                ax[iMomLeg,0].text(text[1],text[2], text[0], horizontalalignment=ha,verticalalignment='top', transform=ax[iMomLeg,0].transAxes,bbox=dict(facecolor='none',pad=2),fontsize=6)


        plt.tight_layout()
        strInOut=['','IN']
        plt.savefig(fp_out+'Fluxes'+strInOut[io_sep]+'.png', format='png', bbox_inches='tight',dpi=600)
        plt.savefig(fp_out+'Fluxes'+strInOut[io_sep]+'.pdf', format='pdf', bbox_inches='tight',dpi=600)

    if (io_sep == 1):
        fig, ax = plt.subplots(nrMom_select, 2, sharex=True,figsize=(8,8*(nrMom_select/5.)**0.8), dpi=600) #,sharey=True
        plt.xlabel(xlabel_Mom[itime_legend])

        #for i in range(nr_MOMs+1):
        for iiMom,iMom in enumerate(iMom_select):
            #inflow and outflow in separate panels
            for i_sim in range(nr_sims):
                if (iplotStdDev == 0):
                    ax[iiMom,0].plot(t_vec[i_sim,:nt_vec[i_sim]],FluxOutMean[i_sim,:nt_vec[i_sim],iMom])
                    ax[iiMom,1].plot(t_vec[i_sim,:nt_vec[i_sim]],FluxOutAccMean[i_sim,:nt_vec[i_sim],iMom])
                else:
                    ax[iiMom,0].errorbar(t_vec[i_sim,:nt_vec[i_sim]],FluxOutMean[i_sim,:nt_vec[i_sim],iMom],yerr=np.squeeze(FluxOut_STD[:,i_sim,:nt_vec[i_sim],iMom]))
                    ax[iiMom,1].errorbar(t_vec[i_sim,:nt_vec[i_sim]],FluxOutAccMean[i_sim,:nt_vec[i_sim],iMom],yerr=np.squeeze(FluxOut_STD[:,i_sim,:nt_vec[i_sim],iMom]))
                    print('OUT INST', iMom, FluxOut_STD[:,i_sim,:nt_vec[i_sim],iMom]/FluxOutMean[i_sim,:nt_vec[i_sim],iMom])
                    print('OUT ACCU', iMom, FluxOutAcc_STD[:,i_sim,:nt_vec[i_sim],iMom]/FluxOutAccMean[i_sim,:nt_vec[i_sim],iMom])
            #GCCif (IREF == 9)
            if (iMom < 4): #skip iMom=4: no SIP number flux in bin modell
                ax[iiMom,0].plot(t_vecREF,FluxesREF[2,iMom,:],'k')
                ax[iiMom,1].plot(t_vecREF,FluxesREF[3,iMom,:],'k')
            #GCCendif /* (IREF == 9) */
            ax[iiMom,0].set_ylabel(ylabel_vec[iMom])

            if (iMom < 4):
                ax[iiMom,0].ticklabel_format(axis='y', style='sci', scilimits=(0,1))
                ax[iiMom,1].ticklabel_format(axis='y', style='sci', scilimits=(0,1))
            if (iMom == 0):
                ax[iiMom,0].set_title('flux kg$^l$/(m$^2$ s)')
                ax[iiMom,1].set_title('acc. flux kg$^l$/m$^2$')
            if (iMom == 4):
                ax[iiMom,0].set_xlabel(xlabel_Mom[itime_legend])
                ax[iiMom,1].set_xlabel(xlabel_Mom[itime_legend])

        iMomLeg = 0
        if (label[0] is not None): ax[iMomLeg,0].legend(label[:nr_sims+iaddsim1DREF],fontsize=6,loc='upper left',ncol=leg_ncol)
        if text is not None:
            print('len(text),type(text)',len(text),type(text))
            if type(text) is str:
                ax[iMomLeg,0].text(0.96, 0.92, text, horizontalalignment='right',verticalalignment='top', transform=ax[iMomLeg,0].transAxes,bbox=dict(facecolor='none',pad=2),fontsize=6)
            else:
                ha='right'
                if (len(text) == 4):
                    ha=text[3]
                ax[iMomLeg,0].text(text[1],text[2], text[0], horizontalalignment=ha,verticalalignment='top', transform=ax[iMomLeg,0].transAxes,bbox=dict(facecolor='none',pad=2),fontsize=6)

        plt.tight_layout()
        plt.savefig(fp_out+'FluxesOUT.png', format='png', bbox_inches='tight',dpi=600)
        plt.savefig(fp_out+'FluxesOUT.pdf', format='pdf', bbox_inches='tight',dpi=600)
#>>>>>>>>>>>>>>>>> end function PlotFluxesTime >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#GCCendif /* (FLUX_time == 1) */

#GCCif (TRACKCENTER == 1)
def PlotCenter(zCenter, nr_SIPs_GB_save, t_vec_MOMsave, name,fp_out="",av=0):
    if (itime_legend == 2): t_vec_MOMsave = t_vec_MOMsave/60

    if (av == 0): [nr_inst,nr_MOMsave,nz,nr_MOMs,nr_eval]= zCenter.shape
    if (av == 1):
        print(zCenter.shape)
        zCenter=np.expand_dims(zCenter, axis=1)
        zCenter=np.expand_dims(zCenter, axis=0)
        nr_SIPs_GB_save=np.expand_dims(nr_SIPs_GB_save, axis=1)
        nr_SIPs_GB_save=np.expand_dims(nr_SIPs_GB_save, axis=0)
        [nr_inst,nr_MOMsave,nz,nr_MOMs,nr_eval]= zCenter.shape

    print('nr_inst,nr_MOMsave,nz,nr_MOMs,nr_eval',nr_inst,nr_MOMsave,nz,nr_MOMs,nr_eval)
    i_inst=0
    for iz in range(nz):
        fig, ax = plt.subplots(5, 1, sharex=True,figsize=(8,8), dpi=600)
        plt.xlabel(xlabel_Mom[itime_legend])
        for iMom in range(nr_MOMs):
            ax[iMom].plot(t_vec_MOMsave,zCenter[i_inst,:,iz,iMom,0],'k')
            ax[iMom].plot(t_vec_MOMsave,zCenter[i_inst,:,iz,iMom,1],'r')
            ax[iMom].plot(t_vec_MOMsave,zCenter[i_inst,:,iz,iMom,2],'g')
            ax[iMom].set_ylim([-1,1])
        ax[4].plot(t_vec_MOMsave,nr_SIPs_GB_save[i_inst,:,iz])
        plt.tight_layout()
        if (av == 0): plt.savefig(fp_out+name+"{0:02}.png".format(iz), format='png', bbox_inches='tight',dpi=600)
        if (av == 1): plt.savefig(fp_out+name+"_av.png", format='png', bbox_inches='tight',dpi=600)
#GCCendif  /* (TRACKCENTER == 1) */

def cycler_dict(n, nfull = 0):

    #from collections import OrderedDict
    #linestyles = OrderedDict(
    #[('solid',               (0, ())),
     #('loosely dotted',      (0, (1, 10))),
     #('dotted',              (0, (1, 5))),
     #('densely dotted',      (0, (1, 1))),

     #('loosely dashed',      (0, (5, 10))),
     #('dashed',              (0, (5, 5))),
     #('densely dashed',      (0, (5, 1))),

     #('loosely dashdotted',  (0, (3, 10, 1, 10))),
     #('dashdotted',          (0, (3, 5, 1, 5))),
     #('densely dashdotted',  (0, (3, 1, 1, 1))),

     #('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     #('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     #('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])
    ls_vec = ['-', ':','--', '-.',(0, (3, 10, 1, 10))]
    ls3 = ['-','--', ':'] # the last style is dotted (for bin model)
    ls4 = ['-','--', '-.', ':'] # the last style is dotted (for bin model)
    ls5 = ['-','--', '-.', (0, (3, 10, 1, 10)), ':'] # the last style is dotted (for bin model)
    #take colors from https://matplotlib.org/examples/color/named_colors.html
    col_vec_reg=['b','g','r','c','m','y','grey', 'orange','olive','skyblue']
    col_vec_reg3=['b','r','g','c','m','y','grey', 'orange','olive','skyblue']
    col_vec_reg2=['olive','b','g','r','c','m','y','grey', 'orange','skyblue']
    col_vec_reg4 = ['orange', 'lime', 'skyblue', 'grey']
    col_vec_reg5 = ['orange', 'skyblue', 'grey']
    col_vec_reg6 = ['lime', 'grey']
    col_vec_reg7 = ['orange', 'lime', 'skyblue', 'brown']
    col_vec_reg8 = ['grey','lime','purple','orange','k' ]
    col_vec_reg9 = ['grey','lime','orange','k' ]
    col_vec_reg1 = ['grey','lime','k' ]

    col_vec_v1=[plt.cm.cool(i) for i in np.linspace(0, 1, 7)]
    marker_vec=list('.sp*+xov^<>1234,')
    if (nfull == 0):
        nfull = n

    col_vec_cbrew = ['#000000']
    if (nfull == 8):
        col_vec_cbrew = ['#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4']
    if (nfull == 7):
        col_vec_cbrew = ['#d73027','#fc8d59','#fee090','#ffffbf','#e0f3f8','#91bfdb','#4575b4']
    if (nfull == 6):
        col_vec_cbrew = ['#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4']
    if (nfull == 5):
        col_vec_cbrew = ['#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6']
    if (nfull == 4):
        col_vec_cbrew = ['#d7191c','#fdae61','#abd9e9','#2c7bb6']
    if (nfull == 3):
        col_vec_cbrew = ['#fc8d59','#ffffbf','#91bfdb']
    if (nfull == 2):
        col_vec_cbrew = ['#d7191c', '#2c7bb6']
    if (nfull == 1):
        col_vec_cbrew = ['#d7191c']

    "YlOrRd"
    col_vec_cbrew2 = ['#000000']
    if (nfull == 7):
        col_vec_cbrew2 = ["#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"]
    if (nfull == 6):
        col_vec_cbrew2 = ["#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#b10026"]
    if (nfull == 5):
        col_vec_cbrew2 = ["#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#b10026"]
    if (nfull == 4):
        col_vec_cbrew2 = ["#feb24c","#fd8d3c","#f03b20","#bd0026"]
    if (nfull == 3):
        col_vec_cbrew2 = ["#fd8d3c","#f03b20","#bd0026"]
    if (nfull == 2):
        col_vec_cbrew2 = ["#fd8d3c","#e31a1c"]
    if (nfull == 1):
        col_vec_cbrew2 = ["#f03b20"]

    # list of possible cyclers (linestyle, markers, color_reg, color_v1)
    cycler_dict  = {'Line'     : cycler(linestyle=ls_vec[:n]),
                    'Marker'   : cycler(marker=marker_vec[:n]),
                    'ColorReg' : cycler(color=col_vec_reg[:n]),
                    'ColorReg4' : cycler(color=col_vec_reg4[:n]),
                    'ColorReg7' : cycler(color=col_vec_reg7[:n]),
                    'ColorReg9' : cycler(color=col_vec_reg9[:n]),
                    'Color_brg' : cycler(color=col_vec_reg3[:n]),
                    'ColorBrew': cycler(color=col_vec_cbrew[:n]),
                    'ColRegP1' : cycler(color=col_vec_reg2[:n]),
                    'ColorV1'  : cycler(color=col_vec_v1[:n]),
                    'C3M3'     : cycler(color=col_vec_reg[:3])*cycler(marker=marker_vec[:3]),
                    'C2M2'     : cycler(color=col_vec_reg[:2])*cycler(marker=marker_vec[:2]),
                    'C2L2'     : cycler(linestyle=ls_vec[:2])*cycler(color=col_vec_reg4[:2]),
                    'C3L2'     : cycler(linestyle=ls_vec[:2])*cycler(color=col_vec_reg4[:3]),
                    'C3L2n'     : cycler(linestyle=ls_vec[:2])*cycler(color=col_vec_reg9[:3]),
                    'L2C3'     : cycler(color=col_vec_reg5[:3])*cycler(linestyle=ls_vec[:2]),
                    'L2C3n'     : cycler(color=col_vec_reg9[:3])*cycler(linestyle=ls_vec[:2]),
                    'L2C2'     : cycler(color=col_vec_reg6)*cycler(linestyle=ls_vec[:2]),
                    'C4L2'     : cycler(linestyle=ls_vec[:2])*cycler(color=col_vec_reg[:4]),
                    'C4L2new'  : cycler(linestyle=ls_vec[:2])*cycler(color=col_vec_reg7[:4]),
                    'C4L3'     : cycler(linestyle=ls3)*cycler(color=col_vec_reg[:4]),
                    'C4L4'     : cycler(linestyle=ls4)*cycler(color=col_vec_reg[:4]),
                    'C4L5'     : cycler(linestyle=ls5)*cycler(color=col_vec_reg[:4]),
                    'C5L3'     : cycler(linestyle=ls3)*cycler(color=col_vec_reg8[:5]),
                    'C5L5'     : cycler(linestyle=ls5)*cycler(color=col_vec_reg[:5]),
                    'C6L2'     : cycler(linestyle=ls_vec[:2])*cycler(color=col_vec_reg[:6]),
                    'C7L2'     : cycler(linestyle=ls_vec[:2])*cycler(color=col_vec_reg[:7]),
                    'L4C3'     : cycler(color=col_vec_reg1[:3])*cycler(linestyle=ls_vec[:4]),
                    'custom1'     : cycler(color=['b','b','b','g','k'])+cycler(linestyle=['-', ':','--','-','-']),
                    'custom1_v2'  : cycler(color=['b','b','b','r','k'])+cycler(linestyle=['-', ':','--','-','-']),
                    'custom1_v3'  : cycler(color=['orange','orange','orange','lime','k'])+cycler(linestyle=['-', ':','--','-','-']),
                    'custom2'     : cycler(color=['k','k', 'orange','orange','orange', 'lime','lime','lime'])+cycler(linestyle=['-', ':','-', ':','--','-', ':','--']),
                    'custom3'  : cycler(color=['grey','orange'])
                    ,
                    'custom3rev'  : cycler(color=['lime','grey']),
                    'G_Gd_O_L'    : cycler(color=['grey','grey','orange','lime'])+cycler(linestyle=['-', ':','-','-']),
                    'G_L_Ld_O_Od': cycler(color=['grey','lime','lime','orange','orange'])+cycler(linestyle=['-','-',':','-',':']),
                    'blacksoliddash': cycler(color=['k','k'])+cycler(linestyle=['-', '--'])
                    }

    return cycler_dict

def flip(items, ncol):
    import itertools
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])
