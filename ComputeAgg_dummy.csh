#!/bin/csh
## first create the file params_gcc.txt andparams.txt, which contain the preprocessor definitions and relevant parameters included in the Python script
## then create pre-processed interpretable source files in the given subdir
## then execute the Python program

#set working directory in which Python program is executed
setenv filepath 'Sims1D/Verification/SediAgg/nz50_inst20_k40_dz10/Instanz_00_09/'
cat > params_gcc.txt << '\eof'
#define COMP  /* don't touch */

#/*  ----------Kernel-------------   */
#define KERNEL 1  /* =0 Golovin kernel, =1 Long kernel, =2 Hall kernel, =3 product kernel */
#define LongK_Options 1 /* =1 tabulated data, =2 computed data */
#define KERNEL_INTPOL 1 /* bei Hall/Long-Kernel: =1 Kernelwert bilinear interpoliert; =0 ohne Interpolation (ein Wert aus LookUp-Table wird eingelesen) , =2 ohne Interpolation, aber lineares statt logarithmisches Gitter */*/

#/* ---------- Initialization -------*/
#define INITV 1  /* = 1 ExpVerteilung, =2 diskrete Verteilung
# /* bei diskreten Verteilung beim Plotten kein Einordnen der SIPs in gegebenes Bin-Gitter, sondern Vielfache der Elementarmasse verwenden. */


#/*----------- AON options ----------*/
#define AGG_MC 2 /* 0 keine Multiple Collections, = 1 Integer Multiple Collection (der IntegerFaktor wird gewrfelt) = 2 Float Multiple Collection */
#define DISCRETE 0 /* 0 Kontinuierliche Fall, = 1 und 2 diskreter Fall (verwende INITV =2 und AGG_MC = 1), Vorgeben einer Elementarmasse, bei 1 Gewichte Integerwerte, bei 2 nur Gewichte 1 */
#define LINEAR 0 /* =0  check all SIP combinations, = 1 linear sampling of SIP combinations */
#define LINEAR_LIMIT 0 /* n_agg can be larger than max(nu_i,nu_j): limit n_agg by max_i,j (=1), no special treatment, can lead to negative weights (=0) */
#define WELLMIXED 0 /* 0: classical 3D well mixed assumption, 1: 2D well mixed assumption, consider overtakes in each gridbox, 2, 3: 2D wellmixed assumption, consider overtakes in full column, 3: additionally accounts for overtakes across lower boundary, only relevant for period BCs */
#define REVERSE 1 /* =1: Reverse SIP processing start from highest SIP not lowest */
#define COUNT_COLLS 0 /* track number of collisions and output to log file*/
#define SPEED_ZERO 1 /* = 1 do no test for zero weight SIPs inside AON */
#define SPEED_VECTOR 0 /* = 1 more vector evaluations */
#define WARN 0 /* =1 warnings in AON-algorithmus */

#/* ---------- 1D options ------------*/
#define COLUMN 1  /* = 0 klassisches BoxModel, = 1 Column model mit klassischer Aggregation in jeder GB */
#define INFLUX_TOP 2  /* influx across top boundary: = 0 no influx , = 1 influx with given SD, =2 periodic BC, outfalling SIPs re-enter domain */
#define PROCESS 0 /* Sedi+Koaleszenz aktiv (=0), Sedimentation (=2) bzw Koaleszenz (=1) ausschalten */
#define RANDOMZ 0 /* random placement of SIPs inside column after each time step */
#define TRACKCENTER 0 /* tracks SIPs centers in each grid box */
#define TRACKOUT 0 /* tracks all SIPs that cross lower boundary */

#/* ----------Plot options -----------*/
#define MOM_meanTD 1   /* plot mean moments, average over all instances and full column, TD = TotalDomain  */
#define GV_meanTD  0   /* plot mean size distribution, average over all instances and full column, TD = TotalDomain  */
#define MOM_prof   0   /* plot vertical profiles of mean moments, average over instances */
#define RZ_scatter 0   /* (r,z)-scatter plot */
#define FLUX_time  0   /* plot fluxes over time */

#/* -------- Reference Solution -----*/
#define IREF 1
# /*  = 0 keine Referenzloesung,
#     = 1 Golovin analytische bzw Wang Long/Hall Referenzloesung,
#     = 2 Alfonso Hall Referenzloesung
#     = 3 Alfonso ProduktK
#     = 9 read Bott/Wang(?) simulation data, set path fp_ref */
'\eof'

cat > params.txt  << '\eof'
import math
#-----------Parameterdefinitionen---------------------------


    #steuert Art und Anzahl der Ausgaben waehrend des Programmdurchlaufs
iPM = 0

    #Zeitschritt
dt=  10.

    #Zeitdauer
Tsim= 3600.  #3600

    #Volumen einer Gitterbox
dV=1. #in m^-3


dV_skal=1  #increase volume by this factor and call init routine as often
dV=dV*dV_skal
dVi=1./dV

    #Anzahl der Instanzen
nr_inst= 10  #00
nr_ins_start=0   # bei verteilter Berechnung der Instanzen kann Seed-Parameter angepasst werden

    #Ausgabe und Speicherung der Simulationsergebnisse
        #Zeitpunkte, an denen Groessenverteilung geplottet werden soll
t_start_GVplot = 0   # in s
t_intervall_GVplot = 200   # in s

        #Zeitpunkte, an denen Momente gespeichert und geplottet werden soll
t_start_MOMsave = 0  # in s
t_intervall_MOMsave = 50 # in s

        #units of time axis in Moment and Flux plots
itime_legend = 2  # = 1 in Sekunden, =2 in Minuten

#GCCif (KERNEL == 0)
b_original=1500   # in (cm^3)*(g^(-1))*(s^(-1))
b_golovin=b_original*(1e-6)*(1e3)*dt*dVi #  SI units (m^3*kg^(-1)*(s^(-1))  * s/m^3
#GCCendif /* (KERNEL == 0) */
#GCCif (KERNEL == 3)
C_original=5.49e10   # in (cm^3)*(g^(-2))*(s^(-1)) bzw auch (m^3*kg^(-2)*(s^(-1))
C_prod=C_original*dt *dVi #  SI units kg^(-2)

#GCCendif /* (KERNEL == 3) */


#GCCif (INITV == 1)
    #Density of water in kg/m**3
rho_w=1e3
    #Konversionskonstante
const_mass2rad=1./((4./3.)*math.pi*rho_w)

    #Festlegung Anfangsverteilung (SingleSIP of ExpVerteilung)
        #physikalische Parameter
            #Mean Droplet radius in m
r0=9.3e-6   # 9.3e-6  #50.e-6  #9.3e-6
            #Massenkonzentration in kg/m^3
LWC=1e-3  #1e-3
        #numerische Parameter
            #Anzahl der Bins pro Massen-Groessenordnung
n10=40
            #Anzahl an Massen-Groessenordnungen
r10=18
            #log10 von unterer Massen-Grenze
min10=-18
            #legt minimale SIP-Troepfchenmasse fest
            # = 0 verwende 0 bzw min10 als untere Grenze, =1 verwende 1.5625um als untere Grenze, =2 verwende 0.6um als untere Grenze
imlow=2
            #minimal zulaessiges SIP-Gewicht relativ zu vorhandenem Maximalgewicht
eta_nu=1e-9

            #maximale SIP-Anzahl
nr_sip_max=15000
    #normalisation constant for SIP droplet masses, usually relevant for the discrete case
skal_m= 1.  # 1 represents an elemental mass, i.e mass of a 17um particle in SI units

#GCCendif /* (INITV == 1) */
#GCCif (INITV == 2)
    #Density of water in kg/m**3
rho_w=1e3
    #Konversionskonstante
const_mass2rad=1./((4./3.)*math.pi*rho_w)

iSIPinit_discrete = '1b'
if (iSIPinit_discrete == '1a'):
    #Maximum number of SIPs
    nr_sip_max=30
    #normalisation constants for concentration and mass
    skal_rad=17.0e-6 # represents radius of unit particle in m
    skal_m=(skal_rad**3.0)/const_mass2rad  # 1 represents an elemental mass, i.e mass of a 17um particle in SI units
    nplot = 40

if (iSIPinit_discrete == '1b'):
    #Maximum number of SIPs
    nr_sip_max=20*dV_skal

    #normalisation constants for concentration and mass
    skal_rad=17.0e-6 # represents radius of unit particle in m
    skal_m=(skal_rad**3.0)/const_mass2rad  # 1 represents an elemental mass, i.e mass of a 17um particle in SI units
    nplot = 40*dV_skal

if (iSIPinit_discrete == '2'):
    #Maximum number of SIPs
    nr_sip_max  = 100
    #normalisation constants for concentration and mass
    skal_rad    = 14.0e-6 # represents radius of unit particle in m
    skal_m      = (skal_rad**3.0)/const_mass2rad  # 1 represents an elemental mass, i.e mass of a 17um particle in SI units
    nplot       = 100
    skal_m2     = skal_m*skal_m
#GCCendif /* (INITV == 2) */

#GCCif (DISCRETE == 0)
    #Definition Bin-Gitter fuer GV-Plot
n10_plot=4
r10_plot=17
nplot=n10_plot*r10_plot
min10_plot=-18
#GCCendif /* (DISCRETE == 0) */

iaddsimREF=0
#GCCif (IREF > 0)
iaddsimREF = 1
#GCCendif /* (IREF > 0) */


#GCCif (COLUMN == 1)
#Anzahl der Gitterboxen
nz = 50  #25

#Art der raeumlichen Initialization
i_init_1D = 6
    #  !    1 = empty domain
    #  !    2 = top GB only
    #  !    3 = top half domain
    #  !    4 = linearly decaying over total column from 2*g_init t
    #  !    5 = linearly decaying over top half from 2*g_init to 0
    #  !    6 = total domain
    #  !    7 = linearly decaying over top quarter from 2*g_init to
    #  !    8 = sin()-hill over top half with 2*g_init max peak
    #  !    9 = sin()-hill over top quarter with 2*g_init max peak

    #Maschenweite in z-Richtung
dz=10. # in m

    #Grundflaeche einer Gitterbox
dA=dV/dz
dAi=1./dA

    #dummy z-position for unused SIPs, do not change
zNan=100000

zCol=nz*dz

#### CHANGE fp_Bott
fp_Bott = "/athome/unte_si/_data/Themen/Aggregation/Bott_1D/"
# fp_ref  = fp_Bott + "testcases_onlysedi/quarterdomSIN_zeroinflow_1km_dz10m_FD3/"
fp_ref  = fp_Bott + "testcases_onlysedi/quarterdomSIN_zeroinflow_1km_dz5m_FD3_iscal16_R050/"

#GCCendif /* (COLUMN == 1) */

'\eof'
#--------------------------------generation of parameter files finished----------------------------------------


# ----------------------------------preprocess and run Python code -------------------------------------------------
# the full source files will be saved in subfolder full_source
# the preprocessed files will be generated and executed in filepath

#set echo

#lege fest, welche Programmdateien zum Ablauf des Programms benoetigt werden (separate Listen fuer Dateien mit/ohne Praeprozessierung)
set list_module_gcc_files = "CompSim.gcc.py AON_Alg.gcc.py Kernel.gcc.py PlotSim.gcc.py Referenzloesung.gcc.py SIPinit.gcc.py Sedimentation.gcc.py"
set list_module_nogcc_files = "Misc.nogcc.py"

    #copy files
mkdir -p $filepath'/full_source/'
mv params.txt params_gcc.txt $filepath'/full_source/'
cp $list_module_gcc_files $list_module_nogcc_files $filepath'/full_source/'
cp $0 $filepath

cd $filepath
echo "Current working directory"
pwd

    #pre-process all Python code that is listed in list_module_gcc_files
    #workaround for application of gcc necessary, as both GCC and Python use # to indicate directive and comments respectively
    #temporarily use a different character sequence instead of # for Python comments
    #additionally GCC directives in the .gcc.py files use #GCC as marker

sed 's/#/\!comment#/g;s/\!comment#GCC/#/g;' full_source/params.txt > params.tmptxt
cat full_source/params_gcc.txt params.tmptxt > params.fpp
gcc -E -P -C params.fpp > params.txt

foreach file_gcc_py  ($list_module_gcc_files)
    echo process file $file_gcc_py
    set base = `basename $file_gcc_py .gcc.py`
    sed -re 's/^.+#GCC/#GCC/' 'full_source/'$file_gcc_py > ${base}.tmp.tmptxt
    sed 's/#/\!comment#/g;s/\!comment#GCC/#/g;' ${base}.tmp.tmptxt > ${base}.tmptxt
    cat full_source/params_gcc.txt ${base}.tmptxt > ${base}.fpp
    gcc -E -P -C ${base}.fpp > ${base}.tmptxt
    sed 's/\!comment#/#/g;' ${base}.tmptxt > ${base}.py
end
sed --in-place 's/\!comment#/#/g;' params.txt

    #simply copy all Python code that is listed in list_module_nogcc_files
foreach file_nogcc_py  ($list_module_nogcc_files)
    echo copy file $file_nogcc_py
    set base = `basename $file_nogcc_py .nogcc.py`
    cp 'full_source/'/${base}.nogcc.py ${base}.py
end

echo $0 > log.txt
hostname >> log.txt
# start program
echo start execution
set startdate = `date`

#### CHANGE: call python 3 with argument CompSim.py
#use local python anaconda installation
/net/lx116/export/home/unte_si/anaconda3/bin/python3 CompSim.py

#mit Profiling
#/net/lx116/export/home/unte_si/anaconda3/bin/python3 -m cProfile -o Profiling.dat -s time CompSim.py

#oder vorinstalliertes Python, kann sich von Rechner zu Rechner unterscheiden oder auch gar nicht existieren
    #python3 CompSim.py
set enddate = `date`

echo "execution start and end time"
echo $startdate
echo $enddate

#clean up
echo "Clean Up"
rm *.fpp *.s *.tmptxt
mkdir -p CodeCompute/full_source
mv *.py *.csh *txt CodeCompute
mv full_source/* CodeCompute/full_source
rmdir full_source

gzip SIP.dat
gzip Moments.dat