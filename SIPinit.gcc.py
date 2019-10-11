import numpy as np
import math

#GCCif (INITV == 1)
import random
import sys

def InitSIP_ExpVert_singleSIP_WS(imlow,n10,r10,min10,eta_nu,xf0,N0,dV,nr_sip_max,dV_skal=1,silent=0):
#imlow      IN : Integer Code
#n10        IN : Anzahl an Bins pro Massendekade
#r10        IN : Anzahl an Massendekaden
#min10      IN : log von unterster Bingrenze
#xf0        IN : Mean mass in kg
#N0         IN : Anfangskonzentration in 1/m**3
#dV         IN : volume of grid box in m**3 or fractional grid box if dV_skal > 1
#nr_sip_max IN: maximale Anzahl an SIPs


    nr_SIPs = 0
    nEK_sip=np.zeros(nr_sip_max)
    mEK_sip=np.zeros(nr_sip_max)

    #Eigenschaften Bingitter
    nr_bins=n10*r10
    mfix=np.zeros(nr_bins)
    mdelta=np.zeros(nr_bins)
    m=np.zeros(nr_bins-1)

    m_low = 0
    if(imlow == 1):
        m_low=(4/3)*math.pi*1e3*(1.5625e-6)**3.0
    if(imlow == 2):
        m_low=(4/3)*math.pi*1e3*(0.6e-6)**3.0

    mfix[0]=10**min10

    iSIP=0

    for i in range(1,nr_bins):
        mfix[i]=10**((i-1)/n10+min10)
        mdelta[i-1]=mfix[i]-mfix[i-1]

    for irep in range(dV_skal):
        nEK_sip_krit_low=0
        nEK_sip_tmp=np.zeros(nr_bins)

        for i in range(1,nr_bins):
            m[i-1]=mfix[i-1]+random.random()*mdelta[i-1]

        s=m/xf0
        for i in range(0,nr_bins-1-1):
            nEK_sip_tmp[i]=N0/xf0*math.exp(-s[i])*mdelta[i]    # entweder hier mit dV multiplizieren or a few lines below

        nEK_sip_krit_low=max(nEK_sip_tmp)*eta_nu

        for i in range(0,nr_bins):
            if(nEK_sip_tmp[i]>nEK_sip_krit_low and m[i]>m_low):
                nEK_sip[iSIP]=nEK_sip_tmp[i]*dV/dV_skal    #mit dV multiplizieren um von Konzentration auf Einheit 1 zu kommen
                mEK_sip[iSIP]=m[i]
                iSIP=iSIP+1
                if(iSIP>nr_sip_max-1):
                    sys.exit("nr_sip_ins(k) ist groesser als nr_sip_max  "+str(iSIP)+"   "+str(i)+"  "+str(irep))
    nr_SIPs=iSIP
    if (silent == 0):
        print('Initialized Exponential Distribution, kappa = ', n10, ' NSIP = ', nr_SIPs, ' dV_skal = ', dV_skal)

    return  nr_SIPs,nEK_sip[0:nr_SIPs],mEK_sip[0:nr_SIPs],nEK_sip_krit_low,mdelta,mfix
# end function InitSIP_ExpVert_singleSIP_WS -----------------------------------------------------------------------------------


def InitSIP_uniform(n10,r10,min10,const_mass2rad, gew=1):
    nr_SIPs=220
    min10=-14
    exponents=min10+np.arange(nr_SIPs)/30
    mEK_sip=10**exponents
    nEK_sip=(np.zeros(nr_SIPs)+1e6) * gew

    #import Misc as FK
    #print(FK.m2r(mEK_sip,const_mass2rad)*1e6)

    #import Misc as FK
    #rr=np.array([2e-3,5.e-2])*1e-2 # in m 
    #nr_SIPs=2
    #mEK_sip=FK.r2m(rr,const_mass2rad)
    #nEK_sip=mEK_sip*0 +1e6
    return  nr_SIPs,nEK_sip,mEK_sip


# end function InitSIP_uniform -----------------------------------------------------------------------------------
#GCCendif /* if (INITV == 1) */


#GCCif (INITV == 2)
#the following statement includes a parameter file via a preprocessor directive
#GCCinclude "params.txt"

def InitSIP_Alfonso(dV_skal=1):
    if (iSIPinit_discrete == '1a'):
        nr_SIPs = nr_sip_max
        nEK_sip_ins=np.zeros(nr_SIPs,dtype='int')
        mEK_sip_ins=np.zeros(nr_SIPs,dtype='int')
        nEK_sip_ins[:]=1
        mEK_sip_ins[:nr_sips_17um]=1
        mEK_sip_ins[nr_sips_17um:]=2

    if (iSIPinit_discrete == '1b'):
        nr_SIPs = nr_sip_max * dV_skal
        nr_SIPs_half=nr_SIPs//2
        nEK_sip_ins=np.zeros(nr_SIPs,dtype='int')
        mEK_sip_ins=np.zeros(nr_SIPs,dtype='int')
        nEK_sip_ins[:nr_SIPs_half*dV_skal]=2
        mEK_sip_ins[:nr_SIPs_half*dV_skal]=1
        nEK_sip_ins[nr_SIPs_half*dV_skal:]=1
        mEK_sip_ins[nr_SIPs_half*dV_skal:]=2
    if (iSIPinit_discrete == '2'):
        nr_SIPs = 100
        nEK_sip_ins=np.zeros(nr_SIPs,dtype='int')
        mEK_sip_ins=np.zeros(nr_SIPs,dtype='int')
        nEK_sip_ins[:]=1
        mEK_sip_ins[:]=1
    
    return nr_SIPs, nEK_sip_ins, mEK_sip_ins

#def InitSIP_Alfonso():
    #nr_SIPs = 40
    #nEK_sip_ins=np.zeros(nr_SIPs,dtype='int')
    #mEK_sip_ins=np.zeros(nr_SIPs,dtype='int')
    #nEK_sip_ins[:]=1
    #mEK_sip_ins[:]=np.arange(40)+1
    #return nr_SIPs, nEK_sip_ins, mEK_sip_ins

#GCCendif /* if (INITV == 2) */