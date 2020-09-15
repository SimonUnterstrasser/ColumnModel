This repository contains Python 3 code of a

coalescence/aggregation column model

The code was used to produce the data and generate figures of the study

Collisional growth in a particle-based cloud microphysical model: Insights from column model simulations using LCM1D (v1.0)

The code was run on Linux system with a C-shell, sed, gcc and Python 3

Copy the repository to any local directory.
Execute "csh ComputeAgg_dummy.csh" for running a box/column model simulation.
Execute "csh PlotAgg_dummy.csh" for a-posteriori plot gereration.
All parameters are set within the csh-scripts.
The folder Documentation contains a documentation.

The revised code contains Python scripts ComputeAgg_dummy.py and PlotAgg_dummy.py. Those are equivalent to the csh counterparts.
The advantage is that no csh script has to be executed any longer. Moreover, calls to sed (which caused trouble in some cases) are replaced by calls to python intrinsic module re.
Then only Python and gcc are necessary for a successful execution of the code.

p.s.:
Before running the program on your machine,
you may need to adapt a few calls/folders in the csh scripts
(the locations can be found by searching for "#### CHANGE"). In particular, fp_Wang and fp_Bott should point to the directories in which the reference bin model data is stored (these data are available from the Zenodo data set with DOI 10.5281/zenodo.4030878)
See file "ProgramVersions.txt" for a list of program versions for which the code was successfully exectued.
