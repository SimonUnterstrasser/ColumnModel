# script that controls the computation procedure
# first create the file params_gcc.txt and params.txt, which contain the preprocessor definitions and relevant parameters included in the Python script
# then create pre-processed interpretable source files in the given subdir
# then execute the Python program

import os
import shutil
import sys
import re
import glob
import platform
import time

#set working directory in which Python program is executed
FILEPATH = 'Sims1D/Verification/SediAgg/nz50_inst20_k05_dz10_py/'
iverboseoutput = 1
ilog_linux = 1

def intermediate_format(tmp):
    tmp = re.sub('^.+#GCC', '#GCC', tmp,flags=re.MULTILINE)
        #include flag, then each line is tested. otherwise only the first line of the multiline string matches the pattern;
        #"By default in python, the "^" and "$" special characters (these characters match the start and end of a line, respectively) only apply to the start and end of the entire string."
    tmp = tmp.replace('#','!comment#')
    tmp = tmp.replace('!comment#GCC','#')
    return(tmp)

def py_format(tmp):
    tmp = tmp.replace('!comment#','#')
    return(tmp)

def py_format_file_inplace(fn):
    f = open(fn,'r')
    c = f.read()
    f.close()
    c_pyf = py_format(c)
    f = open(fn,'w')
    f.write(c_pyf)
    f.close()

CurrDir =  os.getcwd()
print(CurrDir)

output_to_file_params_gcc = """
#define COMP  /* don't touch */
#/* the following PPDs (PreProcessor Directives) control how computations are performed.*/
#/*  ----------Kernel-------------   */
#define KERNEL 1  /* =0 Golovin kernel, =1 Long kernel, =2 Hall kernel, =3 product kernel */
#define LongK_Options 1 /* in case of Long kernel: =1 tabulated data, =2 computed data */
#define KERNEL_INTPOL 0 /* relevant for Hall/Long kernel: =1 kernal values are bilinearly interpolated ; =0 no nterpolation (a particular value is taken from the LookUp-Table) , =2 no interpolation, but linear instead of logarithmic mass grid is used */

#/* ---------- Initialization -------*/
#define INITV 1  /* = 1 Exponential distribution, =2 discrete distribution
# /* if computations are based on (small-size) discrete distributions, then no continuous size distributions are deduced. Instead, a discrete mass grid with multiples of the elemental mass is used. */


#/*----------- AON options ----------*/
#define AGG_MC 2 /* 0 no Multiple Collections, = 1 Integer Multiple Collection (conserves integer weights) = 2 Float Multiple Collection */
#define DISCRETE 0 /* = 0 continous use case, = 1 and 2 discrete use case (specify INITV =2 und AGG_MC = 1), Only multiples of a specified elemental mass are possible; = 1 weights are integer valued, = 2 only weights 1 */
#define LINEAR 0 /* =0  check all SIP combinations (quadratic sampling), = 1 linear sampling of SIP combinations */
#define LINEAR_LIMIT 0 /* options when nu_coll>nu_i,nu_j (the options are only active for LINEAR = 1):
#        =0: no special treatment, can lead to negative weights, nu_coll can be larger than max(nu_i,nu_j):
#        =1: limit nu_coll by 0.99*max(nu_i,j),
#        =2: equal partition SIPs with nu_i,j=0.5 min(nu_i,j) (recommended by Shima during GMD,2020 review)
#            in the present setup with nu-values being floats, equality of n_i and nu_j is not tested,
#            as this is an extremely rare event for the current implementation of SIP ensemble generation.
#            Currently, a collection between two equal weights SIPs generates a zero weight SIP which is removed.
#        =3: use an uneven splitting with 0.4 and 0.6 min(nu_i,j)
# NOTE:  for LINEAR = 0, always the limiter 0.99*max(nu_i,j) is used */
#define WELLMIXED 0 /* 0: classical 3D well mixed assumption, 1: 2D well mixed assumption, consider overtakes in each gridbox, 2, 3: 2D wellmixed assumption, consider overtakes in full column, 3: additionally accounts for overtakes across lower boundary, only relevant for period BCs */
#define REVERSE 1 /* =1: Reverse SIP processing start from highest SIP not lowest, benefitial for WM2D-variant */
#define COUNT_COLLS 0 /* track number of collisions and output to log file*/
#define SPEED_ZERO 1 /* = 1 do no test for zero weight SIPs inside AON */
#define SPEED_VECTOR 0 /* = 1 more vector evaluations */
#define WARN 0 /* =1 warnings in AON-algorithmus */

#/* ---------- 1D options ------------*/
#define COLUMN 1  /* = 0 classical box model, = 2 column model with additional sedimentation */
#define INFLUX_TOP 2  /* influx across top boundary: = 0 no influx , = 1 influx with given SD, =2 periodic BC, outfalling SIPs re-enter domain */
#define PROCESS 0 /* sedimentation and collisional growth are both active (=0), switch off sedimentation (=2) or collisional growth (=1) */
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
# /*  = 0 no reference solution,
#     = 1 analytical Golovin solution or Long/Hall reference solution provided by Wang,
#     = 2 Alfonso Hall kernel reference solution
#     = 3 Alfonso product kernel reference solution
#     = 9 read Bott/Wang(?) simulation data, set path fp_ref */
"""

f = open('params_gcc.txt','w')
f.write(output_to_file_params_gcc)
f.close()

output_to_file_params = """
import math
#----------- parameter definitions--------------------------

    # time step in seconds
dt = 10.

    # simulated time in seconds
Tsim = 3600.  #3600

    #grid box volume
dV = 1. #in m^-3

dV_skal = 1  #increase volume by this factor and call init routine multiple times (increase by factor dV_skal)
dV = dV*dV_skal
dVi = 1./dV

    #number of realisations
nr_inst = 20
nr_ins_start = 0   # each realisation uses a different seed parameter. If a simulation is divided in several subsimulations, then provide a starting index of the current realisations range

    #storage interval of output
        #points in time, at which SIP data is saved and size disitributions can be produced
t_start_GVplot = 0   # in s
t_intervall_GVplot = 200   # in s

        #points in time, at which moment data is saved (usually a finer time grid is used here)
t_start_MOMsave = 0  # in s
t_intervall_MOMsave = 50 # in s

        #units of time axis in Moment and Flux plots
itime_legend = 2  # = 1 in seconds, =2 in minutes

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
    # conversion constant mass to radius
const_mass2rad=1./((4./3.)*math.pi*rho_w)

    #Properties of initial size distribution
        #the SingleSIP-method as described in Unterstrasser et al, 2017 is used to generate a SIP ensemble

        #physical parameters
            #Mean Droplet radius in m
r0=9.3e-6   # 9.3e-6  #50.e-6  #9.3e-6
            #mass concentration in kg/m^3
LWC=1e-3  #1e-3

        #numerical parameters
            #number of bins per mass decade (called kappa in the GMD-papers)
n10=5
            #number of mass decades
r10=18
            #starting mass of bin grid (min10=log10(mass in kg))
min10=-18
            #determines smallest mass/radius of SIP at which SIPs are eventually created
            # = 0: min10 is natural lower threshold, = 1: use 1.5625um, = 2: use 0.6um;
imlow=2
            #determines smallest SIP weight relative to the maximum SIP weight of the SIP ensemble Maximalgewicht
eta_nu=1e-9

            #upper bound for SIP number (determines the size of the SIP arrays)
nr_sip_max=15000

    #normalisation constant for SIP droplet masses, usually relevant for the discrete case
skal_m= 1.  # = 1 => no scaling of SIP droplet masses

#GCCendif /* (INITV == 1) */
#GCCif (INITV == 2)
    #Density of water in kg/m**3
rho_w=1e3
    # conversion constant mass to radius
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
    #definition bin grid for size disitrbution plots (can be chosen independtly of n10)
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
#number of grid boxes in the column
nz = 50  #25

#type of spatial initialisation
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

    #Mesh size
dz=10. # in m

    #base plane of grid box
dA=dV/dz
dAi=1./dA

    #probabilistic SIP influx, average number of SIP per GB and bin
nr_SIPs_skal = 1

    #dummy z-position for unused SIPs, do not change
zNan=100000

zCol=nz*dz

#### CHANGE fp_Bott
fp_Bott = "/athome/unte_si/_data/Themen/Aggregation/Bott_1D/"
# fp_ref  = fp_Bott + "testcases_onlysedi/quarterdomSIN_zeroinflow_1km_dz10m_FD3/"
fp_ref  = fp_Bott + "testcases_onlysedi/quarterdomSIN_zeroinflow_1km_dz5m_FD3_iscal16_R050/"

#GCCendif /* (COLUMN == 1) */

    #controls the extent of logging output onto screen
iPM = 0


"""
f = open('params.txt','w')
f.write(output_to_file_params)
f.close()
#--------------------------------generation of parameter files finished----------------------------------------


# ----------------------------------preprocess and run Python code -------------------------------------------------
# the full source files will be saved in subfolder full_source (FILEPATH_SOURCE)
# the preprocessed files will be generated and executed in folder FILEPATH

#list of source code files (separate lists for files with/without preprocessor directives)
list_module_gcc_files = ('CompSim.gcc.py', 'AON_Alg.gcc.py', 'Kernel.gcc.py', 'PlotSim.gcc.py', 'Referenzloesung.gcc.py', 'SIPinit.gcc.py', 'Sedimentation.gcc.py')
list_module_nogcc_files = ('Misc.nogcc.py',) # a 1-element needs a comma

FILEPATH_SOURCE = os.path.join(FILEPATH,'full_source')

print('current working directory:\n', FILEPATH)

    #copy files
if not os.path.isdir(FILEPATH_SOURCE):
    os.makedirs(FILEPATH_SOURCE)

for file2move in ('params.txt','params_gcc.txt'):
    # if dst is given as folder path plus filename, then move overwrites a potentially existing file with the same name in the dst folder
    shutil.move(file2move, os.path.join(FILEPATH_SOURCE, file2move))

name_current_script = sys.argv[0]
shutil.copy(name_current_script, FILEPATH)

filelist2copy = list_module_gcc_files + \
                list_module_nogcc_files
for file2copy in filelist2copy:
    shutil.copy(file2copy, FILEPATH_SOURCE)

os.chdir(FILEPATH)

#pre-process all Python code that is listed in list_module_gcc_files
#workaround for application of gcc necessary, as both GCC and Python use # to indicate directive and comments respectively
#temporarily use a different character sequence instead of # for Python comments
#additionally, GCC directives in the .gcc.py files use #GCC as marker

tmp = intermediate_format(output_to_file_params)

if (iverboseoutput == 1):
    f = open('params.tmptxt','w') #verbose output
    f.write(tmp) #verbose output
    f.close() #verbose output

f = open('params.fpp','w')
f.write(output_to_file_params_gcc+tmp)
f.close()
os.system('gcc -E -P -C params.fpp > params.txt')

print(os.listdir())
print(os.listdir('full_source'))

for file_gcc_py in list_module_gcc_files:
    print('process file {}'.format(file_gcc_py))
    filebasename = file_gcc_py.replace('.gcc.py','')
    f = open(os.path.join('full_source', file_gcc_py),'r')
    content = f.read()
    f.close()
    content_imf = intermediate_format(content)
        #suffix _imf = intermediate format

    if (iverboseoutput == 1):
        f = open(filebasename+'.tmp.tmptxt','w') #verbose output
        f.write(content_imf) #verbose output
        f.close() #verbose output

    f = open(filebasename+'.fpp','w')
    f.write(output_to_file_params_gcc+content_imf)
    f.close()

    #gcc -E -P -C ${base}.fpp > ${base}.tmptxt
    command_eval = 'gcc -E -P -C {0}.fpp -o {0}.tmptxt'.format(filebasename)
    #command_eval = 'gcc -E -P -C {0}.fpp -o {0}.tmptxt'.format(FILEPATH_OUT+filebasename)
    #print('call preprocessor with command:\n', command_eval)
    os.system(command_eval)

    f = open(filebasename+'.tmptxt','r')
    content = f.read()
    f.close()

    content_pyf = py_format(content)
    f = open(filebasename+'.py','w')
    f.write(content_pyf)
    f.close()

py_format_file_inplace('params.txt')

    #simply copy all Python code that is listed in list_module_gcc_files
for file_nogcc_py in list_module_nogcc_files:
     print('copy file {}'.format(file_nogcc_py))
     filebasename = file_nogcc_py.replace('.nogcc.py','')
     os.rename(os.path.join('full_source',filebasename+'.nogcc.py'), filebasename+'.py')

#preparation of source code finished----------------------------------------

f = open('log.txt','w')
f.write(name_current_script + '\n')
f.write(platform.uname()[1] + '\n')
f.close()

starttime = time.time()

# start program
print('start execution')
#### CHANGE: call python 3 with argument CompSim.py
#use local python anaconda installation
os.system('/net/lx116/export/home/unte_si/anaconda3/bin/python3 CompSim.py')
#with Profiling
#os.system('/net/lx116/export/home/unte_si/anaconda3/bin/python3 -m cProfile -o Profiling.dat -s time CompSim.py')

print("execution start and end time")
endtime = time.time()
print(time.asctime(time.localtime(starttime)),
      time.asctime(time.localtime(endtime  )))
f = open('log.txt','a')
f.write('total computing time in sec: '+ str(int(endtime-starttime)) + '\n')
f.close()

#clean up
print("Clean Up")
filelist_tobedeleted = glob.glob("*.fpp") + glob.glob("*.s") + glob.glob("*.tmptxt")
for file in filelist_tobedeleted:
    os.remove(file)

subdir = os.path.join('CodeCompute','full_source')
if not os.path.isdir(subdir):
    os.makedirs(subdir)

filelist_tobemoved = glob.glob("*.py") + glob.glob("*txt") + glob.glob("full_source/*")

for file2move in filelist_tobemoved:
    # if dst is given as folder path plus filename, then move overwrites a potentially existing file with the same name in the dst folder
    shutil.move(file2move, os.path.join('CodeCompute', file2move))

# filelist_tobemoved =
# for file2move in filelist_tobemoved:
    # # if dst is given as folder path plus filename, then move overwrites a potentially existing file with the same name in the dst folder
    # shutil.move(file2move, os.path.join('CodeCompute',  file2move))

os.rmdir('full_source')
#dest = shutil.move('full_source', 'CodeCompute/', copy_function = shutil.copytree)

os.system('gzip SIP.dat; gzip Moments.dat')