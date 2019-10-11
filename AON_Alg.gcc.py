
import math
#GCCif (SPEED_VECTOR == 0)
import random
#GCCendif  /* (SPEED_VECTOR == 0) */
import numpy as np


#the following statement includes a parameter file via a preprocessor directive
#GCCinclude "params.txt"

#GCCif (WELLMIXED == 0)
def Aggregation(nr_SIPs,nEK_sip_tmp,mEK_sip_tmp,     count_colls,cck,m_low,eta_indize,m_kernel,
                fLog_p=None, fLog_currColls=None, fLog_accColls=None, fLog_currColls_WM2D=None,
                fLog_Combs=None):
#GCCendif /* (WELLMIXED == 0)*/
#GCCif (WELLMIXED > 0)
def Aggregation(nr_SIPs,nEK_sip_tmp,mEK_sip_tmp,z,zn,count_colls,cck,m_low,eta_indize,m_kernel,
                fLog_p=None, fLog_currColls=None, fLog_accColls=None, fLog_currColls_WM2D=None,
                fLog_Combs=None):
    #z old position
    #zn new position
    #it is assumed that the SIP list is sorted by z. smaller z first!
#GCCendif /* (WELLMIXED > 0)*/

    #cck            : Kernel-Matrixwerte
    #nr_SIPs        :
    #nEK_sip_tmp    : Vektor mit SIP-Gewichtsfaktor (will be multiplied by skal_n in KCE)
    #mEK_sip_tmp    : Vektor mit SIP-Troepfchenmassen (will be multiplied by skal_m in KCE)
    #m_low          : kleinste Masse fuer die Kernelwert in Matrix gegeben ist
    #eta_indize     : Anzahl der Kernelwerte in Matrix pro Massendekade
    #m_kernel       : mass grid on which kernel values are given
    #count_colls    : 10-array: counts the number of collections, only active if PPD COUNT_COLLS = 1

#GCCif (COUNT_COLLS == 1)
    ncoll=count_colls.size
    count_colls_before=count_colls.copy()
    p_col_sum=0.
#GCCendif  /* (COUNT_COLLS == 1) */

    ibreak = 0

#GCCif (KERNEL == 1 || KERNEL == 2)
    #GCCif (KERNEL_INTPOL <= 1) /* logarithmic mass bin*/
    indize=np.zeros(mEK_sip_tmp.shape)
    for i in range(0,max(mEK_sip_tmp.shape)):
        if(math.log(mEK_sip_tmp[i]*skal_m)>math.log(m_low)):
            indize[i]=(math.log(mEK_sip_tmp[i]*skal_m)-math.log(m_low))/math.log(eta_indize)
    #print(min(indize),max(indize),nr_SIPs)
    np.clip(indize,0,398,out=indize)

        #GCCif (KERNEL_INTPOL == 1)
    indize_floored=np.floor(indize).astype(int)  # floor liefert float values, in int umwandeln, damit die Werte als Indizes verwendet werden koennen
        #GCCendif /* (KERNEL_INTPOL == 1) */
        #GCCif (KERNEL_INTPOL == 0)
    indize_floored=np.floor(indize+0.5).astype(int)  # floor liefert float values, in int umwandeln, damit die Werte als Indizes verwendet werden koennen
    #nun zentriert
        #GCCendif /* (KERNEL_INTPOL == 0) */
    #GCCendif  /* (KERNEL_INTPOL <= 1) */
    #GCCif (KERNEL_INTPOL == 2) /* linear mass bin*/
    indize_floored=mEK_sip_tmp-1
    #GCCendif  /* (KERNEL_INTPOL ==2) */

    #GCCif (WELLMIXED == 0)
    const=dt*dVi
    #GCCendif /* (WELLMIXED == 0)*/
    #GCCif (WELLMIXED > 0)
    const=dAi
    nOverTakes=0
    nrCombsTested=0
    #GCCendif /* (WELLMIXED > 0)*/
#GCCendif /* (KERNEL == 1 || KERNEL == 2) */

#GCCif (LINEAR == 0)
    nr_Combs = int(0.5*nr_SIPs*(nr_SIPs-1))
    #GCCif (SPEED_VECTOR == 1)
    #nr_Combs = int(0.5*nr_SIPs*(nr_SIPs-1))
    p_vec = np.random.random(nr_Combs)
    iii = 0
    #GCCendif /* (SPEED_VECTOR == 1) */

    #GCCif (REVERSE == 0)
    for i in range(0,nr_SIPs):
        #GCCif (DISCRETE == 2)
        if (nEK_sip_tmp[i] == 0): continue
        #GCCendif /* (DISCRETE == 2) */
        for j in range(i+1,nr_SIPs):

        #GCCif (WELLMIXED > 0)
            #implicitly z[i]<z[j]
            if (zn[i]<zn[j]): continue
            nOverTakes +=1
            #if you have not exited the loop, then SIP i overtook SIP j
        #GCCendif /* (WELLMIXED > 0)*/
    #GCCendif /* (REVERSE == 0) */

    #GCCif (REVERSE == 1)
    for i in range(nr_SIPs-1,0,-1):
        #GCCif (DISCRETE == 2)
        if (nEK_sip_tmp[i] == 0): continue
        #GCCendif /* (DISCRETE == 2) */
        for j in range(i-1,0-1,-1):
            #print(i,j,nr_SIPs)
        #GCCif (WELLMIXED > 0)
            nrCombsTested += 1
            #implicitly z[i]>z[j]
            if (zn[i]>z[j]): break
            if (zn[i]>zn[j]): continue
            nOverTakes +=1
            #if you have not exited the loop, then SIP i overtook SIP j
        #GCCendif /* (WELLMIXED > 0)*/

    #GCCendif /* (REVERSE == 1) */

#GCCendif /* (LINEAR == 0) */

#GCCif (LINEAR == 1)
    index_list=np.argsort(np.random.random(nr_SIPs))
    nr_Combs = math.floor(nr_SIPs/2)
    factor_upscale = 0.5*nr_SIPs*(nr_SIPs-1)/nr_Combs

    #GCCif (SPEED_VECTOR == 1)
    p_vec = np.random.random(nr_Combs)
    iii = 0
    #GCCendif /* (SPEED_VECTOR == 1) */
    for i_comb in range(0,nr_Combs):
            i = index_list[i_comb * 2]
            j = index_list[i_comb * 2 + 1]
            #GCCif (LINEAR_LIMIT == 1)
            i_limit_active=0
            #GCCendif /* (LINEAR_LIMIT == 1) */
#GCCendif /* (LINEAR == 1) */

        #GCCif (KERNEL == 1 || KERNEL == 2)
            #GCCif (KERNEL_INTPOL == 1)
                #Hall/Long Kernel Bilineare Interpolation
            iact=indize[i]
            jact=indize[j]
            iiu=indize_floored[i]
            iju=indize_floored[j]
            cckiuju=cck[iiu,iju]
            cckiujo=cck[iiu,iju+1]
            cckioju=cck[iiu+1,iju]
            cckiojo=cck[iiu+1,iju+1]
            gewi=iact-iiu
            gewj=jact-iju
            cck_value=(cckiuju*(1-gewi)*(1-gewj)+cckiujo*(1-gewi)*gewj+cckioju*gewi*(1-gewj)+cckiojo*gewi*gewj)*const
            #GCCendif  /*(KERNEL_INTPOL == 1) */
            #GCCif (KERNEL_INTPOL == 0 || KERNEL_INTPOL == 2)
            iiu=indize_floored[i]
            iju=indize_floored[j]
            cck_value=cck[iiu,iju]*const
            #GCCendif  /*(KERNEL_INTPOL == 0 || KERNEL_INTPOL == 2) */
        #GCCendif /* (KERNEL == 1 || KERNEL == 2) */
        #GCCif (KERNEL == 0)
            cck_value = b_golovin*(mEK_sip_tmp[i] + mEK_sip_tmp[j])*skal_m
        #GCCendif /* (KERNEL == 0) */
        #GCCif (KERNEL == 3)
            cck_value = C_prod * mEK_sip_tmp[i] * mEK_sip_tmp[j]*skal_m2
        #GCCendif /* (KERNEL == 3) */

            nEK_agg_tmp = cck_value*nEK_sip_tmp[i]*nEK_sip_tmp[j]
            mEK_agg_tmp = mEK_sip_tmp[i]+mEK_sip_tmp[j]
        #GCCif (LINEAR == 1)
            nEK_agg_tmp = nEK_agg_tmp * factor_upscale
            #GCCif (LINEAR_LIMIT == 1)
            tmp_max = max(nEK_sip_tmp[i],nEK_sip_tmp[j])
            if (nEK_agg_tmp > tmp_max):
                #print('oha',nEK_agg_tmp , factor_upscale,nEK_sip_tmp[i],nEK_sip_tmp[j])
                i_limit_active=1
                #GCCif (AGG_MC == 1)
                print('not yet tested!')
                tmp_min=min(nEK_sip_tmp[i],nEK_sip_tmp[j])
                nEK_agg_tmp = math.ceil(tmp_max/tmp_min)*tmp_min
                #GCCendif  /*(AGG_MC == 1) */
                #GCCif (AGG_MC == 2)
                print('A n_agg larger than n_i and n_j', i, j, nEK_agg_tmp, nEK_sip_tmp[i], nEK_sip_tmp[j])
                nEK_agg_tmp = .99 * tmp_max
                    #GCCif (COUNT_COLLS == 1)
                count_colls[ncoll-1] += 1
                    #GCCendif  /* (COUNT_COLLS == 1) */
                #GCCendif  /*(AGG_MC == 2) */
            #GCCendif /* (LINEAR_LIMIT == 1) */
        #GCCendif /* (LINEAR == 1) */

        #GCCif (DISCRETE <= 1)
            tmp_max = max(nEK_sip_tmp[i], nEK_sip_tmp[j])
            if (nEK_agg_tmp > tmp_max):
                print('B n_agg larger than n_i and n_j', i, j, nEK_agg_tmp, nEK_sip_tmp[i], nEK_sip_tmp[j])
                nEK_agg_tmp = 0.99 * tmp_max
                    #GCCif (COUNT_COLLS == 1)
                count_colls[ncoll-1] += 1
                    #GCCendif  /* (COUNT_COLLS == 1) */

            #GCCif (WARN >= 1)
            nEK_save=(nEK_sip_tmp[i],nEK_sip_tmp[j])
            #GCCendif /* (WARN >= 1) */
            tmp_min=min(nEK_sip_tmp[i],nEK_sip_tmp[j])
            #GCCif (SPEED_ZERO == 0)
            if(tmp_min>=0):
                c_mc = nEK_agg_tmp/tmp_min
                #GCCif (COUNT_COLLS == 1)
                p_col_sum += c_mc
                #GCCendif  /* (COUNT_COLLS == 1) */
                c_mcCEIL = math.ceil(c_mc)
                p_col=c_mc-math.floor(c_mc)
            else:
                #GCCif (WARN >= 1)
                print('oha',c_mc,nEK_sip_tmp[i],nEK_sip_tmp[j])
                #GCCendif /* (WARN >= 1) */
                c_mcCEIL=0
                p_col=0
            #GCCendif /* (SPEED_ZERO == 0) */
            #GCCif (SPEED_ZERO == 1)
            c_mc = nEK_agg_tmp/tmp_min
            #GCCif (WARN >= 1)
            if (c_mc > 1e5 or c_mc < 0):
                print('oha',c_mc,nEK_sip_tmp[i],nEK_sip_tmp[j])
            #GCCendif /* (WARN >= 1) */
            #GCCif (COUNT_COLLS == 1)
            p_col_sum += c_mc
            #GCCendif  /* (COUNT_COLLS == 1) */
            c_mcCEIL = math.ceil(c_mc)
            p_col=c_mc-math.floor(c_mc)
            #GCCendif /* (SPEED_ZERO == 1) */
            #Erzeugung Zufallszahl (auf Grundlage des Mersenne Twisters)
            #GCCif (SPEED_VECTOR == 0)
            pkrit=random.random()
            #GCCendif /* (SPEED_VECTOR == 0) */
            #GCCif (SPEED_VECTOR == 1)
            pkrit=p_vec[iii]
            iii+=1
            #GCCendif /* (SPEED_VECTOR == 1) */

            if (c_mcCEIL == 1):
                if(p_col>pkrit):
                    #GCCif (COUNT_COLLS == 1)
                    count_colls[1] += 1
                    #GCCendif  /* (COUNT_COLLS == 1) */
                    if(nEK_sip_tmp[i]<nEK_sip_tmp[j]):
                        i1=i; i2=j
                    else:
                        i1=j; i2=i
                    nEK_sip_tmp[i2]=nEK_sip_tmp[i2]-nEK_sip_tmp[i1]
                    mEK_sip_tmp[i1]=mEK_agg_tmp
                #GCCif (DISCRETE == 1)
                    if (nEK_sip_tmp[i2] == 0): mEK_sip_tmp[i2]=0
                #GCCendif /* (DISCRETE == 1) */
                #GCCif (COUNT_COLLS == 1)
                else:
                    count_colls[0] += 1
                #GCCendif  /* (COUNT_COLLS == 1) */
            elif (c_mcCEIL > 1):
                # Multiple Collection case
                if(nEK_sip_tmp[i]<nEK_sip_tmp[j]):
                    i1=i; i2=j
                else:
                    i1=j; i2=i
                #i1 ist das SIP mit kleinerem Gewichtsfaktor
            #GCCif (AGG_MC == 0)
            #keine Multiple Collection, einfache Kollektion also n_agg wird auf nEK_sip_tmp[i1] gedeckelt
                nEK_sip_tmp[i2]=nEK_sip_tmp[i2]-nEK_sip_tmp[i1]
                mEK_sip_tmp[i1]=mEK_agg_tmp
                #GCCif (COUNT_COLLS == 1)
                index = 2
                #GCCendif  /* (COUNT_COLLS == 1) */
            #GCCendif  /*(AGG_MC == 0) */
            #GCCif (AGG_MC == 1)
            #Integer Multiple Collection, n_agg ist ganzzahliges Vielfaches von nEK_sip_tmp[i1], c_mc wird entweder auf- oder abgerundet
                if(pkrit < p_col):
                    c_mc_pick = c_mcCEIL
                    if (c_mc_pick * nEK_sip_tmp[i1] > tmp_max):
                        print('AGG_MC1: n_agg larger than n_i and n_j, choose c_mcCEIL-1', i, j, nEK_agg_tmp, nEK_sip_tmp[i], nEK_sip_tmp[j])
                        c_mc_pick = c_mcCEIL - 1
                else:
                    c_mc_pick = c_mcCEIL - 1

                mEK_sip_tmp[i1] = (mEK_sip_tmp[i1]+ c_mc_pick*mEK_sip_tmp[i2])
                nEK_sip_tmp[i2] = nEK_sip_tmp[i2] - nEK_sip_tmp[i1]*c_mc_pick
                #GCCif (COUNT_COLLS == 1)
                index=min([c_mc_pick, ncoll-3])
                #GCCendif  /* (COUNT_COLLS == 1) */
            #GCCendif  /*(AGG_MC == 1) */
            #GCCif (AGG_MC == 2)
            #Floating Point Multiple Collection, das berechnete c_mc wird verwendet, diese Operation generiert nicht-ganzzahlige Vielfache der Ausgangsmassen
                mEK_sip_tmp[i1]=(mEK_sip_tmp[i1]+c_mc*mEK_sip_tmp[i2])
                #nEK_sip_tmp[i2]=nEK_sip_tmp[i2]-nEK_sip_tmp[i1]*c_mc
                nEK_sip_tmp[i2]=nEK_sip_tmp[i2]-nEK_agg_tmp
                #GCCif (LINEAR_LIMIT*LINEAR == 1)
                #if (i_limit_active==1): nEK_sip_tmp[i2]=0   # avoids that due to numerical rounding errors a tiny non-zero value is established.
                #GCCendif /* (LINEAR_LIMIT*LINEAR == 1) */
                #GCCif (COUNT_COLLS == 1)
                index=min([math.ceil(c_mc -0.5),ncoll-3])  # 1 < c_mc < 1.5 -> Index 1; 1.5 < c_mc < 2.5 -< Index 2
                #GCCendif  /* (COUNT_COLLS == 1) */
            #GCCendif  /*(AGG_MC == 2) */
                #GCCif (COUNT_COLLS == 1)
                count_colls[index] += 1
            else:
                #this considers case c_mcCEIL = 0
                #happens if nEK_sip_tmp[i] or nEK_sip_tmp[j] is zero OR
                #in case of no interpolation of hydrodynamic kernel, SIPs with similar radii pick diagonal element iiu=iju (which is zero!).
                count_colls[ncoll-2] += 1
                #GCCendif  /* (COUNT_COLLS == 1) */
            #GCCif (WARN >= 1)
            if (nEK_sip_tmp[i] < 0.0) or (nEK_sip_tmp[j] < 0.0):
                print('nEK_sip_tmp regular:', i,j,i1,i2, nEK_sip_tmp[i],nEK_sip_tmp[j],c_mc,c_mcCEIL,nEK_save,nEK_agg_tmp)
            #GCCendif /* (WARN >= 1) */
        #GCCendif /* (DISCRETE <= 1) */

        #GCCif (DISCRETE == 2)
        #only SIP weighting factor of 1 (and 0) occur
        #-> multiple collections are senseless
        #if collection occurs, then set SIP j to zero and exit iteration over inner loop over j

            if (nEK_sip_tmp[j] > 0):  # nEK_sip_tmp[i] is non-zero, only SIP j can be zero
                p_col = nEK_agg_tmp

                #Erzeugung Zufallszahl (auf Grundlage des Mersenne Twisters)
                pkrit=random.random()

                if(p_col>pkrit):
                    #print('i, j', i,j)
                    #print('singleColl  vor: n_i =', nEK_sip_tmp[i],' n_j =', nEK_sip_tmp[j], ' m_i =',mEK_sip_tmp[i], ' m_j =',mEK_sip_tmp[j])
                    nEK_sip_tmp[j]=0
                    mEK_sip_tmp[j]=0
                    mEK_sip_tmp[i]=mEK_agg_tmp
                    #print('singleColl nach: n_i =', nEK_sip_tmp[i],' n_j =', nEK_sip_tmp[j], ' m_i =',mEK_sip_tmp[i], ' m_j =',mEK_sip_tmp[j])
                    #GCCif (LINEAR == 0)
                    continue
                    #GCCendif /* (LINEAR == 0) */
        #GCCendif /* (DISCRETE == 2) */
    #print('mEK_sip_tmp', mEK_sip_tmp)
    #print('nEK_sip_tmp', nEK_sip_tmp)

    ##GCCif (WELLMIXED > 0)
    #nrCombs=0.5*nr_SIPs*(nr_SIPs-1)
    #print(nOverTakes, nrCombs,nOverTakes/nrCombs)
    ##GCCendif /* (WELLMIXED > 0)*/
    nOverTakes_pBC=0
    ##GCCif (WELLMIXED == 3 && INFLUX_TOP == 2)
    # version with full column overtakes and periodic boundary conditions
    #check collisions across lower boundary in this second separate pass of the SIP list
    #print(type(zn))
    #print(zn.shape)
    #print(zn.size)
    #nOverTakes_pBC=0
    #find all SIPs that cross the lower boundary:
    index_list_crossingSIPs = np.where(zn<0)[0]
    #print('index_list_crossingSIPs')
    #print(index_list_crossingSIPs)
    #print('zn',zn[index_list_crossingSIPs])
    nr_SIPs_crossing = len(index_list_crossingSIPs)
    #print('nr_SIPs_crossing',nr_SIPs_crossing)
    #all other SIPs are potential collision partners.
    index_list_candidateSIPs = np.where(zn>0)[0]
    nr_SIPs_candidates=nr_SIPs - nr_SIPs_crossing
    #Note that SIP collisions between two SIPs that both cross the lower boundary was already tested in the original pass
    #print('nr_SIPs_WM3_LBC',nr_SIPs_crossing,nr_SIPs_candidates)
    #print('ABC',type(index_list_crossingSIPs))
    #print('DEF',type(index_list_candidateSIPs))

    for i in index_list_crossingSIPs:
        for j in index_list_candidateSIPs:
            #GCCif (WARN >= 1)
            nEK_save=(nEK_sip_tmp[i],nEK_sip_tmp[j])
            #GCCendif /* (WARN >= 1) */
            #print('i,j', i,j)
            # upward-shift of all crossing SIPs by zCol
            # it is sure that z[i] > zCol > z[j]
            # so the only necesary test is, if zn[i]+zCol< zn[j], skip iteration if this is not the case!
            nrCombsTested += 1
            if ((zn[i]+zCol) > zn[j]): continue
            nOverTakes_pBC += 1
            #print('i,j', i,j,zn[i]+zCol, zn[j])
            #GCCif (KERNEL_INTPOL == 1)
                #Hall/Long Kernel Bilineare Interpolation
            iact=indize[i]
            jact=indize[j]
            iiu=indize_floored[i]
            iju=indize_floored[j]
            cckiuju=cck[iiu,iju]
            cckiujo=cck[iiu,iju+1]
            cckioju=cck[iiu+1,iju]
            cckiojo=cck[iiu+1,iju+1]
            gewi=iact-iiu
            gewj=jact-iju
            cck_value=(cckiuju*(1-gewi)*(1-gewj)+cckiujo*(1-gewi)*gewj+cckioju*gewi*(1-gewj)+cckiojo*gewi*gewj)*const
            #GCCendif  /*(KERNEL_INTPOL == 1) */
            #GCCif (KERNEL_INTPOL == 0 || KERNEL_INTPOL == 2)
            iiu=indize_floored[i]
            iju=indize_floored[j]
            cck_value=cck[iiu,iju]*const
            #GCCendif  /*(KERNEL_INTPOL == 0 || KERNEL_INTPOL == 2) */

            nEK_agg_tmp = cck_value*nEK_sip_tmp[i]*nEK_sip_tmp[j]
            mEK_agg_tmp = mEK_sip_tmp[i]+mEK_sip_tmp[j]

        #GCCif (DISCRETE <= 1)
            tmp_max = max(nEK_sip_tmp[i], nEK_sip_tmp[j])
            if (nEK_agg_tmp > tmp_max):
                print('pBCs, n_agg larger than n_i and n_j', i, j, nEK_agg_tmp, nEK_sip_tmp[i], nEK_sip_tmp[j])
                nEK_agg_tmp = 0.99 * tmp_max
                count_colls[ncoll-1] += 1

            tmp_min=min(nEK_sip_tmp[i],nEK_sip_tmp[j])
            #GCCif (SPEED_ZERO == 0)
            if(tmp_min>0):
                c_mc = nEK_agg_tmp/tmp_min
                #GCCif (COUNT_COLLS == 1)
                p_col_sum += c_mc
                #GCCendif  /* (COUNT_COLLS == 1) */
                c_mcCEIL = math.ceil(c_mc)
                p_col=c_mc-math.floor(c_mc)
            else:
                #GCCif (WARN >= 1)
                print('oha',c_mc,nEK_sip_tmp[i],nEK_sip_tmp[j])
                #GCCendif /* (WARN >= 1) */
                c_mcCEIL=0
                p_col=0
            #GCCendif /* (SPEED_ZERO == 0) */
            #GCCif (SPEED_ZERO == 1)
            c_mc = nEK_agg_tmp/tmp_min
            #GCCif (WARN >= 1)
            if (c_mc > 1e5 or c_mc < 0):
                print('oha',c_mc,nEK_sip_tmp[i],nEK_sip_tmp[j])
            #GCCendif /* (WARN >= 1) */

            c_mcCEIL = math.ceil(c_mc)
            p_col=c_mc-math.floor(c_mc)
            #GCCif (COUNT_COLLS == 1)
            p_col_sum += c_mc
            #GCCendif  /* (COUNT_COLLS == 1) */
            #GCCendif /* (SPEED_ZERO == 1) */

            #Erzeugung Zufallszahl (auf Grundlage des Mersenne Twisters)
            pkrit=np.random.random()
            if (c_mcCEIL == 1):
                if(p_col>pkrit):
                    #GCCif (COUNT_COLLS == 1)
                    count_colls[1] += 1
                    #GCCendif  /* (COUNT_COLLS == 1) */
                    if(nEK_sip_tmp[i]<nEK_sip_tmp[j]):
                        i1=i; i2=j
                    else:
                        i1=j; i2=i
                    nEK_sip_tmp[i2]=nEK_sip_tmp[i2]-nEK_sip_tmp[i1]
                    mEK_sip_tmp[i1]=mEK_agg_tmp
                    #GCCif (DISCRETE == 1)
                        if (nEK_sip_tmp[i2] == 0): mEK_sip_tmp[i2]=0
                    #GCCendif /* (DISCRETE == 1) */
                #GCCif (COUNT_COLLS == 1)
                else:
                    count_colls[0] += 1
                #GCCendif  /* (COUNT_COLLS == 1) */
            elif (c_mcCEIL > 1):
                # Multiple Collection case
                if(nEK_sip_tmp[i]<nEK_sip_tmp[j]):
                    i1=i; i2=j
                else:
                    i1=j; i2=i
                #i1 ist das SIP mit kleinerem Gewichtsfaktor
            #GCCif (AGG_MC == 0)
            #keine Multiple Collection, einfache Kollektion also n_agg wird auf nEK_sip_tmp[i1] gedeckelt
                #print('Oha', nEK_sip_tmp[i],nEK_sip_tmp[j],nEK_agg_tmp,c_mc)

                nEK_sip_tmp[i2]=nEK_sip_tmp[i2]-nEK_sip_tmp[i1]
                mEK_sip_tmp[i1]=mEK_agg_tmp
                #GCCif (COUNT_COLLS == 1)
                index = 2
                #GCCendif  /* (COUNT_COLLS == 1) */
            #GCCendif  /*(AGG_MC == 0) */
            #GCCif (AGG_MC == 1)
            #Integer Multiple Collection, n_agg ist ganzzahliges Vielfaches von nEK_sip_tmp[i1], c_mc wird entweder auf- oder abgerundet
                if(pkrit<p_col):
                    c_mc_pick=c_mcCEIL
                    if (c_mc_pick * nEK_sip_tmp[i1] > tmp_max):
                        print('AGG_MC1: n_agg larger than n_i and n_j, choose c_mcCEIL-1', i, j, nEK_agg_tmp, nEK_sip_tmp[i], nEK_sip_tmp[j])
                        c_mc_pick = c_mcCEIL - 1
                else:
                    c_mc_pick=c_mcCEIL - 1

                mEK_sip_tmp[i1]=(mEK_sip_tmp[i1]+c_mc_pick*mEK_sip_tmp[i2])
                nEK_sip_tmp[i2]=nEK_sip_tmp[i2]-nEK_sip_tmp[i1]*c_mc_pick
                #GCCif (COUNT_COLLS == 1)
                index=min([c_mc_pick,ncoll-3])
                #GCCendif  /* (COUNT_COLLS == 1) */
            #GCCendif  /*(AGG_MC == 1) */
            #GCCif (AGG_MC == 2)
            #Floating Point Multiple Collection, das berechnete c_mc wird verwendet, diese Operation generiert nicht-ganzzahlige Vielfache der Ausgangsmassen
                mEK_sip_tmp[i1]=(mEK_sip_tmp[i1]+c_mc*mEK_sip_tmp[i2])
                #nEK_sip_tmp[i2]=nEK_sip_tmp[i2]-nEK_sip_tmp[i1]*c_mc
                nEK_sip_tmp[i2]=nEK_sip_tmp[i2]-nEK_agg_tmp
                #GCCif (LINEAR_LIMIT*LINEAR == 1)
                if (i_limit_active==1): nEK_sip_tmp[i2]=0
                #GCCendif /* (LINEAR_LIMIT*LINEAR == 1) */
                #GCCif (COUNT_COLLS == 1)
                index=min([math.ceil(c_mc -0.5),ncoll-3])  # 1 < c_mc < 1.5 -> Index 2; 1.5 < c_mc < 2.5 -< Index 3
                #GCCendif  /* (COUNT_COLLS == 1) */
            #GCCendif  /*(AGG_MC == 2) */
                #GCCif (COUNT_COLLS == 1)
                count_colls[index] += 1
            else:
                #this considers case c_mcCEIL = 0
                #happens if nEK_sip_tmp[i] or nEK_sip_tmp[j] is zero OR
                #in case of no interpolation of hydrodynamic kernel, SIPs with similar radii pick diagonal element iiu=iju (which is zero!).
                count_colls[ncoll-2] += 1
                #GCCendif  /* (COUNT_COLLS == 1) */
            #GCCif (WARN >= 1)
            if (nEK_sip_tmp[i] < 0.0) or (nEK_sip_tmp[j] < 0.0):
                print('nEK_sip_tmp regular:', i,j, nEK_sip_tmp[i],nEK_sip_tmp[j],c_mc,nEK_save,nEK_agg_tmp)
            #GCCendif /* (WARN >= 1) */

        #GCCendif /* (DISCRETE <= 1) */
    #print('nr_SIPs_WM3_LBC',nr_SIPs_crossing, nr_SIPs_candidates,
          #nOverTakes_pBC, nOverTakes_pBC/(nr_SIPs_crossing*nr_SIPs_candidates))

    ##GCCendif /* (WELLMIXED == 3 && INFLUX_TOP == 2)*/

    ##GCCif (WELLMIXED > 0)
    nrCombs=nr_Combs  #   0.5*nr_SIPs*(nr_SIPs-1)
    ##GCCif (WELLMIXED < 3)
    #print(nOverTakes, nrCombs,nOverTakes/nrCombs)
    ##GCCelse
    #print(nOverTakes,nOverTakes_pBC, nrCombs,nOverTakes/nrCombs,nOverTakes_pBC/nrCombs)
    ##GCCendif /*(WELLMIXED < 3) */
    ##GCCendif /* (WELLMIXED > 0)*/

    #GCCif (COUNT_COLLS == 1)
    i_fLog_currColls = 0
    if (fLog_currColls is None):
        fLog_currColls= open('log_currColls.dat','a')
        i_fLog_currColls = 1

    i_fLog_accColls = 0
    if (fLog_accColls is None):
        fLog_accColls= open('log_accColls.dat','a')
        i_fLog_accColls = 1

    nrCombs = nr_Combs  #0.5*nr_SIPs*(nr_SIPs-1)
    #print('mean p',p_col_sum/nrCombs,p_col_sum,nrCombs)
    #count_colls[0]=nrCombs-(count_colls[:ncoll-1].sum()-count_colls_before[:ncoll-1].sum())
    fLog_Combs.write("{}\n".format(nr_Combs))
    cc=(count_colls-count_colls_before,count_colls)
    fLog_cc=(fLog_currColls,fLog_accColls)
    for cc,fLog_cc in zip(cc,fLog_cc):
        fLog_cc.write(" ".join("{}".format(x) for x in cc)+'\n')
        #print('CHECK',cc.sum(),nrCombs,nr_SIPs)
        if (nrCombs > 0):
            cc_frac = cc/nrCombs
        else:
            cc_frac = cc * 0.0
        fLog_cc.write(" ".join("{:.4}".format(x) for x in cc_frac)+'\n')

    if (i_fLog_currColls == 1):
        fLog_currColls.close()
    if (i_fLog_accColls == 1):
        fLog_accColls.close()

    p_OT = 1.0
    if (nrCombs > 0):
        p_coll = p_col_sum / nrCombs
    else:
        p_coll = 0.0

    ##GCCif (WELLMIXED > 0)
    i_fLog_currColls_WM2D = 0
    if (fLog_currColls_WM2D is None):
        fLog_currColls_WM2D= open('log_currColls_WM2D.dat','a')
        i_fLog_currColls_WM2D = 1
    cc=count_colls-count_colls_before
    cc_sum=cc.sum()
    nOT_tot=nOverTakes+nOverTakes_pBC

    if (nrCombs > 0):
        p_OT = nOT_tot / nrCombs
    else:
        p_OT = 0.0

    if (nOT_tot > 0): p_coll = p_col_sum / nOT_tot
    #print('nOT: ',cc_sum,nOT_tot)
    #print('mean p WM2D:', p_col_sum/nOT_tot)
    vec = np.array((nOverTakes,nrCombsTested,nOT_tot, nrCombs))
    fLog_currColls_WM2D.write(" ".join("{}".format(x) for x in vec)+'\n')
    if (nrCombs > 0):
        vec = vec/nrCombs
    fLog_currColls_WM2D.write(" ".join("{}".format(x) for x in vec)+'\n')
    if (i_fLog_currColls_WM2D == 1):
        fLog_currColls_WM2D.close()
    ##GCCendif /* (WELLMIXED > 0)*/

    i_fLog_p = 0
    if (fLog_p is None):
        fLog_p= open('log_p.dat','a')
        i_fLog_p = 1

    vec=(p_OT, p_coll, p_coll * p_OT)
    outp_str = " ".join("{:.3e}".format(x) for x in vec) + "  "
    outp_str += " ".join("{}".format(x) for x in [int(nrCombs),nr_SIPs])
    fLog_p.write(outp_str+'\n')
    if (i_fLog_p == 1):
        fLog_p.close()
    #GCCendif  /* (COUNT_COLLS == 1) */

    #GCCif (WARN >= 1)
    min_nEK = min(nEK_sip_tmp)
    if (min_nEK < 0.0):
        print('min_nEK: ', min_nEK)
    #GCCendif /* (WARN >= 1) */


    return ibreak
#end subroutine Aggregation(dt,nr_SIPs,nEK_sip_tmp,mEK_sip_tmp,dV,count_colls,cck,m_low,eta_indize)



