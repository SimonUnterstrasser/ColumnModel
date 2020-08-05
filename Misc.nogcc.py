import numpy as np

"""
Mapping SIP to Bin
"""

#Mapping of SIP ensemble onto a bin grid
#SIP ensemble (of one realisation) with nr_SIPs SIP defined by nEK_sip and mEK_sip
#Properties of the bin grid are specified by n10_plot, r10_plot, min10_plot

def MapSIPBin(nEK_sip,mEK_sip,nr_SIPs,n10,r10,min10):
    n= n10 * r10
    mfix=np.zeros(n)
    mdelta=np.zeros(n)
    mfix[0]=10**(min10)
    for i in range(1,n):
        mfix[i]=10**((i-1)/n10+min10)
        mdelta[i-1]=mfix[i]-mfix[i-1]
    if (mfix[-1] < max(mEK_sip[0:nr_SIPs])):
        print("Note that the defined bin grid does not cover the full SIP ensemble ", mfix[-1] , max(mEK_sip[0:nr_SIPs]))

    #sort SIPs by size and check to which bin each SIP contributes
    z=1
    mEK_bin_tot=np.zeros(n)
    nEK_bin_ord=np.zeros(n)
    nEK_sip_tmp=nEK_sip[0:nr_SIPs]
    mEK_sip_tmp=mEK_sip[0:nr_SIPs]

    per=np.argsort(mEK_sip_tmp)

    for i in range(0,nr_SIPs):
        while(mfix[z]<mEK_sip_tmp[per[i]]):
            z=z+1

        nEK_bin_ord[z-1]=nEK_bin_ord[z-1]+nEK_sip_tmp[per[i]]/mdelta[z-1]
        #Hier entweder Binmittelpunkte (bei log muss noch mfix[i]=10^((i-1+0.5)/n10+min10)) oder so aehnlich abgeaendert werden
        if(nEK_bin_ord[z-1]!=0):
            mEK_bin_tot[z-1]=mEK_bin_tot[z-1]+nEK_sip_tmp[per[i]]*mEK_sip_tmp[per[i]]
        else:
            mEK_bin_tot[z-1]=(mfix[z]-mfix[z-1]/2)+mfix[z-1]

    mEK_bin_ord=np.zeros(n)
    for i in range(0,n):
        if(nEK_bin_ord[i]!=0):
            mEK_bin_ord[i]=mEK_bin_tot[i]/(nEK_bin_ord[i]*mdelta[i])
    return nEK_bin_ord,mEK_bin_ord,mEK_bin_tot,mdelta,z

#end subroutine MapSIPBin(nEK_sip,mEK_sip,nr_SIPs,n10,r10,min10):

def CountSIPs(nEK_sip,mEK_sip,nr_SIPs,nplot):
    # only for discrete applications with integer weights
    nEK_bin_plot=np.zeros(nplot)
    #print(max(mEK_sip),nplot)
    if (max(mEK_sip) > nplot):
        print('---------- Achtung: ', max(mEK_sip), nplot)
    for i in range(nr_SIPs):
        nEK_bin_plot[int(mEK_sip[i])-1] += nEK_sip[i]

    return nEK_bin_plot
#end subroutine CountSIPs(nEK_sip,mEK_sip,nr_SIPs,nplot)

def CIO_MOMmean(data_in=None,fp_in=None,fp_out=None,skal_m=1.0,dV=1.0,ikeep1D=None,igetSTD=0,iexcludetop=0):
    #CIO = "Compute or Input/Read mean Moments" and Output the data

    #Compute mean moments
    #Input: moment data of all realisations (given either by array data_in or read from files in folder  fp_in)
    #Output: return an array with averaged moments. Additionally write the data to a file, if  folder name fp_out is specified.

    #in box model simulations data_in.ndim is 3 and the average is taken over all instances => data_out.ndim=2
    #in colum model simulations data_in.ndim is 4
    #   if ikeep1D=None, take average over instances and column => data_out.ndim=2
    #   if ikeep1D=1,    take average over instances, produces mean vertical profiles of the moments => data_out.ndim=3

    #igetSTD = 1 return mean and standard deviation
    #igetSTD = 2 return and 10 and 90 percentile
    print('Compute mean moments')

    if (fp_in is not None):
        ndim_data=3
        fu = open(fp_in + 'Moments_meta.dat','r')
        dV = int(fu.readline())
        skal_m = int(fu.readline())
        nr_inst = int(fu.readline())
        nr_MOMsave = int(fu.readline())
        # in 1D read additional line with nz-value, then set ndim_data=4
        #t_vec_MOMsave= np.array(fu.readline().split(), dtype='float' )
        fu.close()
        fu = open(fp_in + 'Moments.dat','rb')
        data_in=np.loadtxt(fu).reshape((nr_inst,nr_MOMsave,4))
        fu.close()
        #print(data_in.shape)
    else:
        if (data_in is not None):
            ndim_data=data_in.ndim
            if (ndim_data == 3): [nr_inst,nt_MOMsave,nr_Mom]=data_in.shape
            if (ndim_data == 4): [nr_inst,nt_MOMsave,nz,nr_Mom]=data_in.shape
        else:
            print('supply either data or fp_in')


    if (ndim_data == 4 and ikeep1D is None):
        #case: column model data and compute totals over full column
        #two dimension is averaged over: over nr_inst and nz
        if (iexcludetop == 0):
            data_out_ColIntegr=np.mean(data_in,axis=2)  # additionally average over column
        else:
            nz_smallDom=int(nz/40*iexcludetop)
            print('nz_smallDom, nz: ', nz_smallDom, nz)
            data_out_ColIntegr=np.mean(data_in[:,:,:nz_smallDom,:],axis=2)  # additionally average over column
        data_out = np.mean(data_out_ColIntegr,axis=0)
        if (igetSTD == 1):
            MOM_StdDev = np.expand_dims(np.std(data_out_ColIntegr, axis=0), axis=0)
        if (igetSTD == 2):
            MOM_StdDev = np.abs(np.percentile(data_out_ColIntegr, [10,90], axis=0)-data_out) # given as positive deviations from the mean value
    else:
        #covers cases:
            # box model data
            # column model data, output are profiles
        data_out=np.mean(data_in,axis=0)
        if (igetSTD == 1):
            MOM_StdDev = np.expand_dims(np.std(data_in, axis=0), axis=0)
        if (igetSTD == 2):
            MOM_StdDev = np.abs(np.percentile(data_in, [10,90], axis=0)-data_out) # given as positive deviations from the mean value

    #print(data_out)
    if (fp_out is not None):
        fu = open(fp_out + 'MomentsMean.dat','wb')
        #data_in=np.savetxt(fu,data_out, fmt='%1.6e')
        np.savetxt(fu,data_out, fmt='%1.6e')
        fu.close()

    if (ndim_data == 4 and ikeep1D is not None):
        for i in range(4): data_out[:,:,i] = data_out[:,:,i]*skal_m**i/dV
    else:
        for i in range(4): data_out[:,i] = data_out[:,i]*skal_m**i/dV

    if (igetSTD > 0):
        return data_out, MOM_StdDev
    else:
        return data_out

def Moments(nr_inst,nEK_sip_ins,mEK_sip_ins,nr_moment):
    moment_inst=np.array(np.zeros(nr_inst))
    for k in range(0,nr_inst):
        moment_inst[k]=sum(nEK_sip_ins[k,:]*(mEK_sip_ins[k,:]**nr_moment))
    moment_all=sum(moment_inst[:])/nr_inst

    return moment_all, moment_inst

##############################################################################

def Moments_k0_3(nEK_sip_ins,mEK_sip_ins):
    moments_k=np.zeros(4)
    for k in range(0,4):
        moments_k[k]=sum(nEK_sip_ins*mEK_sip_ins**k)

    return moments_k
##############################################################################

def get_first_second_highest_SIPmass(mEK_sip_tmp):
    [nr_inst,nr_GVplot,nr_SIPs]= mEK_sip_tmp.shape
    maxM  = np.zeros([nr_inst,nr_GVplot])
    max2M = np.zeros([nr_inst,nr_GVplot])
    for i_time in range(nr_GVplot):
        for i_inst in range(nr_inst):
            max2M[i_inst,i_time], maxM[i_inst,i_time] = np.partition(mEK_sip_tmp[i_inst,i_time,:], -2)[-2:]

    maxM_mean  = np.mean(maxM,axis=0)
    max2M_mean = np.mean(max2M,axis=0)

    return maxM_mean, max2M_mean

def get_isep(zEK_sip_ins,zGBsep,nz):
    # works only, if a final SIP with z=zNan>zGBsep[nz] is present
    iSIP_GBsep = np.zeros(nz+1,dtype='int')
    nr_SIPs_GB = np.zeros(nz,dtype='int')
    iSIP=0
    iSIP_GBsep[0]=0
    #iz=0
    #print(iz,iSIP_GBsep[iz],zGBsep[iz])
    for iz in range(1,nz+1):
        while (zEK_sip_ins[iSIP] < zGBsep[iz]):
            iSIP+=1
        iSIP_GBsep[iz]=iSIP
        #print(iz,iSIP_GBsep[iz],zGBsep[iz])
    #print(iSIP_GBsep)
    nr_SIPs_GB = iSIP_GBsep[1:nz+1]-iSIP_GBsep[0:nz]
    #print(nr_SIPs_GB)
    return iSIP_GBsep,nr_SIPs_GB

def get_isep2(zEK_sip_ins,zGBsep,nz,nrSIPs):
    iSIP_GBsep = np.zeros(nz+1,dtype='int')
    nr_SIPs_GB = np.zeros(nz,dtype='int')
    iSIP=0
    iSIP_GBsep[0]=0
    iz=0
    #print(iz,iSIP_GBsep[iz],zGBsep[iz])
    for iz in range(1,nz+1):
        while (zEK_sip_ins[iSIP] < zGBsep[iz]):
            if (iSIP < nrSIPs-1):
                iSIP += 1
            else:
                break
        #print(iz,iSIP_GBsep[iz],zGBsep[iz])
        iSIP_GBsep[iz]=iSIP

    nr_SIPs_GB = iSIP_GBsep[1:nz+1]-iSIP_GBsep[0:nz]
    return iSIP_GBsep,nr_SIPs_GB

def m2r(mass_vec,const_mass2rad):
    #mass_vec  IN: mass in kg
    #         OUT: radius in m
    return (mass_vec*const_mass2rad)**(1./3.)

def r2m(r_vec,const_mass2rad):
    #mass_vec  IN: radius in m
    #         OUT: mass in kg
    return (r_vec**3.0)/const_mass2rad

def trackcenters(n,m,z,iSIP_GBsep,nr_SIPs_GB,nz,dz):
    # computes centroid of SIPs position inside grid boxes
    tmp   = np.zeros([2,nz,4])
    nom   = np.zeros([2,4])
    denom = np.zeros([2,4])
    for iz in range(0,nz):
        if (nr_SIPs_GB[iz] > 1):
            ia=iSIP_GBsep[iz]
            ie=iSIP_GBsep[iz+1]
            znorm= (z[ia:ie]/dz)-iz
            for iMom in range(4):
                weightIC=n[ia:ie]*(m[ia:ie]**iMom)
                weightSIP=m[ia:ie]**iMom
                a=(znorm*weightIC ).sum()
                b= weightIC.sum()
                c=(znorm*weightSIP).sum()
                d=weightSIP.sum()
                tmp[0,iz,iMom] = a/b
                tmp[1,iz,iMom] = c/d
                nom[0,iMom]+=a
                nom[1,iMom]+=c
                denom[0,iMom]+=b
                denom[1,iMom]+=d
    return tmp,nom,denom

def as_list(x):
    # converts, e.g., an integer value into an iterable list with a single integer value
    if type(x) is list:
        return x
    else:
        return [x]


def openfile(fn,options='r'):
    import os
    import gzip
    if os.path.isfile(fn+'.gz'):
        dat=gzip.open(fn+'.gz',options)
    elif os.path.isfile(fn):
        dat=open(fn,options)
    else:
        print('Datei nicht gefunden: ', fn+'(.gz)')
        return None
    return dat

def filename(fn,fn_base=''):
    import os
    igzip   = 0
    inormal = 0
    ifile = 0
    if os.path.isfile(fn_base+fn+'.gz'):
        igzip = 1
        ifile = 1
    if os.path.isfile(fn_base+fn):
        inormal = 1
        ifile = 2

    if (igzip + inormal == 2):
        'Both file types exist. Use Zip file (=1), use normal file (=2)'
        print(fn_base+fn)
        ifile = int(input("Both file types exist. Use Zip file (=1), use normal file (=2):  "))

    if (ifile > 0):
        suffix = [".gz",""]
        fn_return = fn+suffix[ifile-1]
        return fn_return

    print("File does not exist")
    fn_return = None
    return fn_return

def smooth(y, box_pts):
    #example call: plot(x, smooth(y,3), 'r-', lw=2)
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
