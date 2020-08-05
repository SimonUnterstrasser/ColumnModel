# script that controls the plotting procedure
# first create the files params_gcc.txt and params_plot.txt, which contain the preprocessor definitions and parameter definitions
# then create pre-processed interpretable source files in the given subdir
# then execute the Python program

import os
import shutil
import sys
import re
import glob

def intermediate_format(tmp):
    tmp = re.sub('^.+#GCC', '#GCC', tmp,flags=re.MULTILINE)
        #include flag re.MULTILINE, then each line is tested. otherwise only the first line of the multiline string matches the pattern;
        #"By default in python, the ‘^’ and ‘$’ special characters (these characters match the start and end of a line, respectively) only apply to the start and end of the entire string."
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

IMODE  = 1
# = 1 generate plots of a single simulation. The plots will be stored in the subfolder "plot" of the according simulation folder FILEPATH
# = 2 generate summary plots of multiple simulations
    # FILEPATH determines the output path of the figure file and the list of simulations is given in list filepath_input
#the selection of the plots to be produced is given by preprocessor variables inside params_plot_gcc.txt

iverboseoutput = 1

CurrDir =  os.getcwd()

FILEPATH = 'Sims1D/Verification/SediAgg/nz50_inst20_k05_dz10_py/'

print(CurrDir)

output_to_file_params_plot_gcc = """
#define PLOT  /* don't touch */
#/* the following PPDs control which plots are generated.*/
#define MOM_meanTD 1   /* TD = TotalDomain, average over all instances, optionally plot all ensemble members
#                       = 1 plot mean moments, average over full column:
#                            in units kg^k m^3
#                       = 2 plot integrated moments, integrate over z-direction, (multiply by Lz relative to option 1)
#                            in units kg^k m^2
#                       = 3 plot total moments, simply sum up SIP-data
#                            in units kg^k  */
#define FINAL 0 /* plot only moment value at the end of the simulation */
#define GV_meanTD  0   /* plot mean size distribution, average over all instances and full column, TD = TotalDomain  */
#define GV_out     0   /* plot size distribution over all outfalling SIPs, average over all instances */
#define MOM_prof   1   /* plot vertical profiles of mean moments, average over instances */
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
#define IMODE {}
#define DISCRETE 0
#define IREF 9
# /*  = 0 no reference solution,
#     = 1 analytical Golovin solution or Long/Hall reference solution provided by Wang,
#     = 2 Alfonso Hall kernel reference solution
#     = 3 Alfonso product kernel reference solution
#     = 9 read Bott/Wang(?) simulation data, set path fp_ref */
#define COLUMN 1  /* = 0 classical box model, = 2 column model with additional sedimentation */
""".format(IMODE)

#output_to_file_params_plot_gcc.replace("%IMODE%")

f = open('params_plot_gcc.txt','w') #verbose output
f.write(output_to_file_params_plot_gcc) #verbose output
f.close() #verbose output


output_to_file_params_plot = """
import math

#----------- parameter definitions---------------------------
 # plot moments of a single file in subfolder "plot/"

#GCCif (IMODE == 1)
filepath = '../'
label=None
text=None
title=None
Times_list=[[0,1800,3600]]
Time_cycler='ColorReg'
Sim_cycler='Line'
#GCCendif

#GCCif (IMODE > 1)

#grafik008/SediAgg
fp_sf='{}'+ '/Sims1D/'
filepath_input =    [ fp_sf + 'Verification/NoSedi_LinSamp_LL2/nz50_inst20_k'+ i + '_dz10/' for i in ['05','10','20','40','60', '100', '200' ]] + \
                    [ fp_sf + 'Verification/SediAgg_LinSamp_LL2/nz50_inst20_k'+ i + '_dz10/' for i in ['05','10','20','40','60', '100','200' ]] + \
                    [ fp_sf + 'Verification/NoSedi_INTPOL1/nz50_inst20_k'+ i + '_dz10/' for i in ['05','10','20','40','60', '100', '200' ]] + \
                    [ fp_sf + 'Verification/SediAgg_INTPOL1/nz50_inst20_k'+ i + '_dz10/' for i in ['05','10','20','40','60', '100','200' ]] + \
                    [ fp_sf + 'Verification/NoSedi_LinSamp_LL2/nz50_inst20_k'+ i + '_dz10_dt1/' for i in ['05','10','20','40','60', '100' ,'200']] + \
                    [ fp_sf + 'Verification/SediAgg_LinSamp_LL2/nz50_inst20_k'+ i + '_dz10_dt1/' for i in ['05','10','20','40','60', '100' ,'200']] + \
                    [ fp_sf + 'Verification/NoSedi_LinSamp_LL2/nz50_inst20_k40_dz10'+ i+ '/' for i in ['_dt1','_dt2','_dt5','','_dt20']] + \
                    [ fp_sf + 'Verification/SediAgg_LinSamp_LL2/nz50_inst20_k40_dz10'+ i+ '/' for i in ['_dt1','_dt2','_dt5','','_dt20']]+ \
                    [ fp_sf + 'Verification/NoSedi/nz50_inst20_k40_dz10'+ i+ '/' for i in ['_dt2','_dt5','','_dt20','_dt50','_dt100']] + \
                    [ fp_sf + 'Verification/SediAgg/nz50_inst20_k40_dz10'+ i+ '/' for i in ['_dt2','_dt5','','_dt20','_dt50','_dt100']]+ \
                    [ fp_sf + 'Verification/NoSedi_LinSamp_LL2/nz50_inst20_k100_dz10'+ i+ '/' for i in ['_dt1','_dt2','_dt5','','_dt20']] + \
                    [ fp_sf + 'Verification/SediAgg_LinSamp_LL2/nz50_inst20_k100_dz10'+ i+ '/' for i in ['_dt1','_dt2','_dt5','','_dt20']]


pos_x = 0.23
pos_y = 0.08
pos_xNo=0.04
pos_x2 = 0.65
# text_vec = [[r'$nz = 10$'+ '\nnoSedi\nLinSamp' , pos_x,pos_y], [r'$nz = 10$'+ '\nLinSamp', pos_xNo,pos_y] , \
# [r'$nz = 10$'+'\n' + '$dt = 1$s'+ '\nnoSedi\nLinSamp', pos_xNo,pos_y], [r'$nz = 10$'+'\n' + '$dt = 1$s'+ '\nLinSamp', pos_xNo,pos_y] , \
# [r'$nz = 10$'+ '\nLinSamp' , pos_x,pos_y], [r'$nz = 10$'+ '\n noSedi \nLinSamp', pos_xNo,pos_y] , \
# [r'$nz = 10$'+ '\n' + r'$\kappa=100$' + '\n noSedi \nLinSamp', pos_xNo,pos_y] , [r'$nz = 10$'+ '\n' + r'$\kappa=100$' + '\nLinSamp', pos_xNo,pos_y] ]
text_vec = [ 'noSedi, $dt = 10$s',  '$dt = 10$s', \
             'nS, $dt = 10$s, QuadSamp',  '$dt = 10$s, QuadSamp', \
             'noSedi, $dt = 1$s'  ,  '$dt = 1$s',  \
             'noSedi, $\kappa=40$' ,  '$\kappa=40$', \
             'nS, $\kappa=40$, QuadSamp' ,  '$\kappa=40$, QuadSamp', \
             'noSedi, $\kappa=100$' ,  '$\kappa=100$']

titleMAIN ='AON WM3D, LinSamp'

iPanel=[0,0,0,0,0,0,0,
        1,1,1,1,1,1,1,
        2,2,2,2,2,2,2,
        3,3,3,3,3,3,3,
        4,4,4,4,4,4,4,
        5,5,5,5,5,5,5,
        6,6,6,6,6,
        7,7,7,7,7,
        8,8,8,8,8,
        9,9,9,9,9,
        10,10,10,10,10,10,
        11,11,11,11,11,11]

iseries_indices=[0,7,14,21,28,35,42,47,52,58,64,69,74]
nr_series = 12
iseriestype=  [0,0,0,0,0,0,1,1,1,1,1,1]
iseries_shift=[0,0,0,0,0,0,0,0,1,1,0,0]
nr_seriestypes=2
series_vecs=[['5','10','20','40','60','100','200'],['1s','2s','5s','10s','20s','50s','100s']]
iperm = [0,1,4,5,2,3, 6,7,10,11,8,9]

# iPanelAnnotate=[0.91,0.11]
# iExtraLegend = 2
ilastGVonly= 1
# choose cyclers from 'Line', 'Marker', 'ColorReg', 'ColorV1'
# Time_cycler:
Time_cycler='Line'
# Sim_cycler:
Sim_cycler='L2C3'
#define a list of lists with times to be displayed in a single plot
Times_list=[[0,1800,3600]]
#end of block
#GCCendif

iMom_select_const = [0,2]
itime_legend = 2  # = 1 in seconds, =2 in minutes

    # density of water in kg/m**3
rho_w=1e3
    # conversion constant mass to radius
const_mass2rad=1./((4./3.)*math.pi*rho_w)


#GCCif (DISCRETE == 0)
TintervalGV = 200
    #dV: grid box volume
dV=1. #in m^-3
dVi=1./dV
    #definition of bin grid for SD plots
n10_plot=5
r10_plot=17
nplot=n10_plot*r10_plot
min10_plot=-18
#ilastGVonly = 0

nr_sip_max    = 40000    #maximum SIP number
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
    #dV grid box volume
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
nr_sip_max    = 30    #maximum SIP number
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
nr_sip_max    = 100    #maximum SIP number
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
""".format(CurrDir)

#sed --in-place=tmp "s:%CURRDIR%:${CurrDir}:g" params_plot.txt

f = open('params_plot.txt','w') #verbose output
f.write(output_to_file_params_plot) #verbose output
f.close() #verbose output


#--------------------------------generation of parameter files finished----------------------------------------


# ----------------------------------preprocess and run Python code -------------------------------------------------
# the full source files will be saved in subfolder full_source
# the preprocessed files will be generated and executed in filepath

#list of source code files (separate lists for files with/without preprocessor directives)
list_module_gcc_files = ('Kernel.gcc.py',  'Referenzloesung.gcc.py', 'SIPinit.gcc.py', 'PlotResults.gcc.py', 'PlotSim.gcc.py', 'InputOutput.gcc.py')
list_module_nogcc_files = ('Misc.nogcc.py',) # a 1-element needs a comma

FILEPATH_OUT = FILEPATH + '/'
if (IMODE == 1):
    FILEPATH_OUT = FILEPATH + '/plot/'

print('current working directory:\n', FILEPATH_OUT)

    #copy files
if not os.path.isdir(FILEPATH_OUT):
    os.makedirs(FILEPATH_OUT)

for file2move in ('params_plot.txt','params_plot_gcc.txt'):
    # if dst is given as folder path plus filename, then move overwrites a potentially existing file with the same name in the dst folder
    shutil.move(file2move, os.path.join(FILEPATH_OUT, file2move))

name_current_script =sys.argv[0]
filelist2copy = list_module_gcc_files + \
                list_module_nogcc_files + \
                (name_current_script, )
for file2copy in filelist2copy:
    print('copy file: ', file2copy, FILEPATH_OUT)
    shutil.copy(file2copy, FILEPATH_OUT)

os.chdir(FILEPATH_OUT)

#pre-process all Python code that is listed in list_module_gcc_files
#workaround for application of gcc necessary, as both GCC and Python use # to indicate directive and comments respectively
#temporarily use a different character sequence instead of # for Python comments
#additionally, GCC directives in the .gcc.py files use #GCC as marker
tmp = intermediate_format(output_to_file_params_plot)

if (iverboseoutput == 1):
    f = open('params_plot.tmptxt','w') #verbose output
    f.write(tmp) #verbose output
    f.close() #verbose output


f = open('params_plot.fpp','w')
f.write(output_to_file_params_plot_gcc+tmp)
f.close()
os.system('gcc -E -P -C params_plot.fpp > params_plot.txt')

for file_gcc_py in list_module_gcc_files:
    print('process file {}'.format(file_gcc_py))
    filebasename = file_gcc_py.replace('.gcc.py','')
    f = open(file_gcc_py,'r')
    content = f.read()
    f.close()
    content_imf = intermediate_format(content)
        #suffix _imf = intermediate format

    if (iverboseoutput == 1):
        f = open(filebasename+'.tmptxt','w') #verbose output
        f.write(content_imf) #verbose output
        f.close() #verbose output

    f = open(filebasename+'.fpp','w')
    f.write(output_to_file_params_plot_gcc+content_imf)
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
    os.remove(file_gcc_py)

py_format_file_inplace('params_plot.txt')

    #simply copy all Python source files that are listed in list_module_gcc_files
for file_nogcc_py in list_module_nogcc_files:
     print('copy file {}'.format(file_nogcc_py))
     filebasename = file_nogcc_py.replace('.nogcc.py','')
     os.rename(filebasename+'.nogcc.py', filebasename+'.py')

#preparation of source code finished----------------------------------------

# start program
print('start execution')

#### CHANGE: call python 3 with argument PlotResults.py
os.system('/net/lx116/export/home/unte_si/anaconda3/bin/python3 PlotResults.py')
# os.system('/net/lx116/export/home/unte_si/anaconda3/bin/python3 caller_function.py')

#clean up
print("Clean Up")
filelist_tobedeleted = glob.glob("*.fpp") + glob.glob("*.s") + glob.glob("*.tmptxt")
for file in filelist_tobedeleted:
    os.remove(file)
