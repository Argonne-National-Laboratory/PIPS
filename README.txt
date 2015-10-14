PIPS - suite of parallel optimization solvers mainly for stochastic optimization problems
consisting of the following solvers:
 i.   PIPS-IPM - parallel distributed memory interior-point for stochastic LPs and convex QPs
 ii.  PIPS-S   - distributed memory implementation of the revised simplex method
 iii. PIPS-NLP - parallel interior-point equipped with filter line-search for structured NLPs

###############################################################################################
# LICENSE
###############################################################################################
See LICENSE file.

###############################################################################################
# CONTRIBUTIONS
###############################################################################################
PIPS-IPM
  Cosmin G. Petra - Argonne National Laboratory
  Miles Lubin - Argonne National Laboratory
  Naiyuan Chiang - Argonne National Laboratory

PIPS-S
  Miles Lubin - Argonne National Laboratory
  Cosmin G. Petra - Argonne National Laboratory
  Geoffrey Oxberry - Lawrence Livermore National Laboratory

PIPS-NLP 
 Naiyuan Chiang - Argonne National Laboratory
 Victor Zavala - Argonne National Laboratory and Univ. of Wisconsin-Madison
 Cosmin G. Petra - Argonne National Laboratory	 

###############################################################################################
# INSTALATION Instructions
###############################################################################################

1. Install package wget, cmake, mpich2, and boost.
You can get them via the following command (xxx stands for the name of the package):
In Linux(Ubuntu): apt-get install xxxx

2. Go to the following folders and run the script wgetXXX.sh
ThirdPartyLibs/ASL  
ThirdPartyLibs/CBC 
ThirdPartyLibs/ConicBundle   
ThirdPartyLibs/METIS
ThirdPartyLibs/UMFPACK
For an example, use command "sh wgetASL.sh" in the folder ThirdPartyLibs/ASL  

3. Download MA27 and MA57 from HSL and put the source code in the correct folder. 
(See ThirdPartyLibs/MA27/README.txt and ThirdPartyLibs/MA57/README.txt for more details.)

4. Assuming we are trying to install PIPS in the folder PIPSMAINPATH/build_pips, where 
PIPSMAINPATH is the root instalation folder, use the following commands in the PIPSMAINPATH
folder to configure and install PIPS:
mkdir build_pips
cd build_pips
cmake ..
make

5. The build system will install executables from three resources: PIPS-IPM, PIPS-S and PIPS-NLP. 
For the usages of these executables, please follow the README.txt files in the corresponding sub-folder.

###############################################################################################
# ACKNOWLEDGMENTS
###############################################################################################

PIPS has been developed under the financial support of: 
- Department of Energy, Office of Advanced Scientific Computing Research
- Department of Energy, Early Career Program 
- Department of Energy, Office of Electricity Delivery and Energy Reliability

