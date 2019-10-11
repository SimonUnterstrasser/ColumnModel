#!/bin/csh
## script that controls the plotting procedures
## first create the files params_gcc.txt and params_plot.txt, which contain the preprocessor definitions and parameter deifinitions
## then create pre-processed interpretable source files in the given subdir
## then execute the Python program

setenv IMODE 2
# = 1 plotte Abbildungen einer einzelnen Simulation in plot-Unterordner
     # setze FILEPATH auf SimOrdner
# = 2 plotte Abbildungen mehrerer Simulationen
    # setze FILEPATH auf Ausgabepfad und die Liste an Simordnern in filepath_input
#the selection of the plots to be produced is given the preprocessor variables inside params_plot_gcc.txt

set CurrDir =  `pwd`  #"/athome/unte_si/_data/Themen/Aggregation/PythonAggregation/"

setenv FILEPATH 'SimsLong/LongK_ExpV_dt1s_FPMC/'
setenv FILEPATH 'SimsLong/LongK_ExpV_dt10s_FPMC/'
setenv FILEPATH 'SimsLong/LongK_ExpV_dt10s/'
setenv FILEPATH 'Sims1D/Test_08/'
setenv FILEPATH 'Sims1D/fulldom_zeroinflow_1km_dz20m/'
#setenv FILEPATH 'Plots/grafik008/onlyAgg/kappa_var_nz50_inst020/'
setenv FILEPATH 'Plots/grafik008/SediAgg/dt_var/'

echo $CurrDir

cat > params_plot_gcc.txt << '\eof'
#define PLOT  /* don't touch */
#/* the following PPDs control which plots are generated.*/
#define MOM_meanTD 1   /* TD = TotalDomain, average over all instances, optionally plot all ensemble members
#                       = 1 plot mean moments, average over full column:
#                            in units kg^k m^3
#                       = 2 plot integrated moments, integrate over z-direction, (multiply by Lz relative to option 1)
#                            in units kg^k m^2
#                       = 3 plot total moments, simply sum up SIP-data
#                            in units kg^k  */

#define GV_meanTD  0   /* plot mean size distribution, average over all instances and full column, TD = TotalDomain  */
#define GV_out     0   /* plot size distribution over all outfalling SIPs, average over all instances */
#define MOM_prof   0   /* plot vertical profiles of mean moments, average over instances */
#define RZ_scatter 0   /* (r,z)-scatter plot */
#define FLUX_time  0   /* plot fluxes over time */
#define rSIGMA_time 0 /* plot relative STD of largest SIP mass + mass of largest and second largest SIP  */

#define addSIPplot 0 /* create extra SIP profile plot, add extra panel with nr_SIP(t) in Mom(t)-plot
#                       = 1 profile plot shows N_SIP per 100m segment
#                       = 2 profile plot shows N_SIP,GB          */
#define Prof_Panel 0 /* plot all profile plots of selected quantities in a row */
#define MOM_Panel2Dsims 0 /*  = 0 create a matrix of various moment quantities (rows) and various subsets of simulations (column)
#                             = 1 create a 2D matrix of various subsets of simulations, display only one quantity, define nr_col, nr_rows */
#define MOM_NORMmass 0 /* normalise first moment (=mass) by initial value */
#define MOM_skal_m 0  /* = 1 show moments in terms of elemental mass */
#define STD 0 /* = 1 add standard deviation to plots of MOM_meanTD, MOM_prof and FLUX_time
#                = 2 add 10, 90 percentile to plots of MOM_meanTD, MOM_prof and FLUX_time*/
#define KERNEL 1  /* =0 Golovin kernel, =1 Long kernel, =2 Hall kernel, =3 Produkt Kernel*/
#define IMODE %IMODE%
#define DISCRETE 0
#define IREF 9
# /*  = 0 keine Referenzloesung,
#     = 1 Golovin analytische bzw Wang Long/Hall Referenzloesung,
#     = 2 Alfonso Hall Referenzloesung
#     = 3 Alfonso ProduktK
#     = 9 read Bott/Wang(?) simulation data, set path fp_ref */
#define COLUMN 1  /* = 0 klassisches BoxModel, = 1 Column model mit klassischer Aggregation in jeder GB */
'\eof'

sed --in-place=tmp "s:%IMODE%:${IMODE}:g" params_plot_gcc.txt

cat > params_plot.txt  << '\eof'
import math
#-----------Parameterdefinitionen---------------------------
 # plotte GV und Momente einer einzelnen Simulation in plot1 Unterordner

#GCCif (IMODE == 1)
filepath = '../'
label=None
text=None
title=None
Time_cycler='ColorReg'
Sim_cycler='Line'
#GCCendif

#GCCif (IMODE > 1)

#grafik008/SediAgg
fp_sf='%CURRDIR%'+ '/Sims1D/'
# strs = ["" for x in range(size)]nz10_inst10_k40_dz10
# filepath_input =  [fp_sf + 'Verification/NoSedi_INTPOL1/nz50_inst20_k'+ i + '_dz10/' for i in ['05','10','20','40', '60']  ]
# label=['05','10','20','40', '60']
filepath_input =  [fp_sf + 'Verification/SediAgg/nz50_inst20_k40_dz10'+ i + '/' for i in ['_dt2', '_dt5', '', '_dt20', '_dt50', '_dt100']  ]
label=['2', '5', '10', '20', '50', '100']

title=r'$dt$'
ilastGVonly= 1
# choose cyclers from 'Line', 'Marker', 'ColorReg', 'ColorV1'
# Time_cycler:
Time_cycler='Line'
# Sim_cycler:
Sim_cycler='ColorReg'
#define a list of lists with times to be displayed in a single plot
Times_list=[[0,1800,3600]]
#end of block
#GCCendif

iMom_select_const = [0,2]
itime_legend = 2  # = 1 in Sekunden, =2 in Minuten

    #Density of water in kg/m**3
rho_w=1e3
    #Konversionskonstante
const_mass2rad=1./((4./3.)*math.pi*rho_w)


#GCCif (DISCRETE == 0)
TintervalGV = 200
    #Volumen einer Gitterbox
dV=1. #in m^-3
dVi=1./dV
    #Definition Bin-Gitter fuer GV-Plot
n10_plot=5
r10_plot=17
nplot=n10_plot*r10_plot
min10_plot=-18
#ilastGVonly = 0

nr_sip_max    = 40000    #maximale SIP-Anzahl
nr_inst_max   = 100
nr_GVplot_max = 37

iaddsimREF=0
iaddsim1DREF=0
#GCCif (IREF > 0)
iaddsimREF = 1
#GCCif (IREF == 1)
if label is not None: label.append('Bin')
#GCCendif /* (IREF == 1) */
#GCCendif /* (IREF > 0) */
#GCCif (IREF == 1)
iPanelREF=[0,1,2,3,4,5,6,7]
#GCCendif /* (IREF == 1) */

#GCCif (IREF == 9)
iaddsim1DREF = 1
#### CHANGE fp_Bott 
fp_Bott = "/athome/unte_si/_data/Themen/Aggregation/Bott_1D/"
# fp_ref  = fp_Bott + "testcases_onlysedi/quarterdomSIN_zeroinflow_1km_dz10m_FD3/"
#fp_ref  = fp_Bott + "testcases_onlysedi/quarterdomSIN_zeroinflow_1km_dz5m_FD3_iscal16_R050/"
fp_ref  = fp_Bott + "testcases/quarterdomSIN_zeroinflow_1km_dz10m_FD3/"
fp_ref  = fp_Bott + "testcases/fulldom_zeroinflow_1km_dz20m_FD3/"
dVREF = 1
fp_ref  = [fp_Bott + "VerificationLong/nz10_scal4_dz100_dt10s" + i + "/" for i in ['', '_LWC1.5_rfix','_LWC2.0_rfix','_LWC2.5_rfix'] ]
#### CHANGE fp_Wang
fp_Wang = "/athome/unte_si/_data/Themen/Aggregation/Wang_1D/"
fp_ref  = [fp_Wang + "testcases_onlyagg/test04_iscal16/"]

nzREFmax=10
ntREFmax=200
nrsimREF=len(fp_ref)
#labelREF=['LWC1.5', 'LWC2.0','LWC2.5' ]
#GCCendif /* (IREF == 9) */
#GCCendif /* (DISCRETE == 0) */

#GCCif (DISCRETE >= 1)
    #Volumen einer Gitterbox
dV=1e-6 #in m^-3
dVi=1./dV

#GCCif (IREF == 2)
TintervalGV = 250
iaddsimREF = 1
if label is not None: label.append('Alfonso')

    #normalisation constant for mass
skal_rad=17.0e-6 # represents radius of unit particle in m
skal_m=(skal_rad**3.0)/const_mass2rad  # 1 represents an elemental mass, i.e mass of a 17um particle in SI units
nplot = 40
nr_sip_max    = 30    #maximale SIP-Anzahl
nr_inst_max   = 10000
nr_GVplot_max = 11
#GCCendif /* (IREF == 2) */
#GCCif (IREF == 3)
iaddsimREF = 1
if label is not None: label.append('Alfonso')
    #normalisation constant for mass
skal_rad    = 14.0e-6 # represents radius of unit particle in m
skal_m      = (skal_rad**3.0)/const_mass2rad  # 1 represents an elemental mass, i.e mass of a 17um particle in SI units
nplot       = 100
nr_sip_max    = 100    #maximale SIP-Anzahl
nr_inst_max   = 10000
nr_GVplot_max = 41
#GCCendif /* (IREF == 3) */

ilastGVonly = 1


#GCCendif /* (DISCRETE >= 1) */


#GCCif (KERNEL == 0)
b_original=1500   # in (cm^3)*(g^(-1))*(s^(-1))
b_golovin=b_original*(1e-6)*(1e3)*dt*dVi #  SI units (m^3*kg^(-1)*(s^(-1))  * s/m^3
#GCCendif /* (KERNEL == 0) */
#GCCif (KERNEL == 3)
C_original=5.49e10   # in (cm^3)*(g^(-2))*(s^(-1)) bzw auch (m^3*kg^(-2)*(s^(-1))

#GCCendif /* (KERNEL == 3) */
'\eof'

sed --in-place=tmp "s:%CURRDIR%:${CurrDir}:g" params_plot.txt
#--------------------------------generation of parameter files finished----------------------------------------


# ----------------------------------preprocess and run Python code -------------------------------------------------
# the full source files will be saved in subfolder full_source
# the preprocessed files will be generated and executed in filepath

#set echo
#lege fest, welche Programmdateien zum Ablauf des Programms benoetigt werden (separate Listen fuer Dateien mit/ohne Praeprozessierung)
set list_module_gcc_files = "Kernel.gcc.py Referenzloesung.gcc.py SIPinit.gcc.py PlotResults.gcc.py PlotSim.gcc.py InputOutput.gcc.py"
set list_module_nogcc_files = "Misc.nogcc.py"

setenv FILEPATH_OUT $FILEPATH'/'
if ($IMODE == 1) setenv FILEPATH_OUT $FILEPATH'/plot/'

echo 'current working directory'
echo $FILEPATH_OUT
    #copy files
mkdir -p $FILEPATH_OUT
mv params_plot.txt params_plot_gcc.txt $FILEPATH_OUT
cp $list_module_gcc_files $FILEPATH_OUT
cp $list_module_nogcc_files $FILEPATH_OUT
cp $0 $FILEPATH_OUT
pwd
cd $FILEPATH_OUT
pwd
    #pre-process all Python code that is listed in list_module_gcc_files
    #workaround for application of gcc necessary, as both GCC and Python use # to indicate directive and comments respectively
    #temporarily use a different character sequence instead of # for Python comments
    #additionally, GCC directives in the .gcc.py files use #GCC as marker

sed -re 's/^.+#GCC/#GCC/' params_plot.txt > params_plot.tmp.tmptxt
sed 's/#/\!comment#/g;s/\!comment#GCC/#/g;' params_plot.tmp.tmptxt > params_plot.tmptxt
cat params_plot_gcc.txt params_plot.tmptxt > params_plot.fpp
gcc -E -P -C params_plot.fpp > params_plot.txt

foreach file_gcc_py  ($list_module_gcc_files)
    echo process file $file_gcc_py
    set base = `basename $file_gcc_py .gcc.py`
    sed -re 's/^.+#GCC/#GCC/' $file_gcc_py > ${base}.tmp.tmptxt
    sed 's/#/\!comment#/g;s/\!comment#GCC/#/g;' ${base}.tmp.tmptxt > ${base}.tmptxt
    cat params_plot_gcc.txt ${base}.tmptxt > ${base}.fpp
    gcc -E -P -C ${base}.fpp > ${base}.tmptxt
    sed 's/\!comment#/#/g;' ${base}.tmptxt > ${base}.py
    rm $file_gcc_py
end
sed --in-place 's/\!comment#/#/g;' params_plot.txt

    #simply copy all Python code that is listed in list_module_gcc_files
foreach file_nogcc_py  ($list_module_nogcc_files)
     echo copy file $file_nogcc_py
     set base = `basename $file_nogcc_py .nogcc.py`
     mv ${base}.nogcc.py ${base}.py
end

# start program
echo "start execution"
# set startdate = `date`
pwd

#use local python anaconda installation
#### CHANGE: call python 3 with argument PlotResults.py
/net/lx116/export/home/unte_si/anaconda3/bin/python3 PlotResults.py
#oder vorinstalliertes Python, kann sich von Rechner zu Rechner unterscheiden oder auch gar nicht existieren
#python3 PlotResults.py

# set enddate = `date`

#echo "execution start and end time"
# echo $startdate
# echo $enddate

#clean up
rm *.fpp *.s *.tmptxt
