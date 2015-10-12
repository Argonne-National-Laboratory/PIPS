1. Install package wget, cmake, mpich2, and boost.
You can get them via the following command (xxx stands for the name of the package):
In Linux(Ubuntu): apt-get install xxxx

2. Go to the following folders and run the script wgetXXX.sh
ThirdPartyLib/ASL  
ThirdPartyLib/CBC 
ThirdPartyLib/ConicBundle   
ThirdPartyLib/METIS
ThirdPartyLib/UMFPACK
For an example, use command "sh wgetASL.sh" in the folder ThirdPartyLib/ASL  

3. Download MA27 and MA57 from HSL and put the source code in the correct folder. (See ThirdPartyLob/MA27/README.txt and ThirdPartyLob/MA57/README.txt for more details.)

4. Assuming we are trying to install PIPS in the folder PIPSMAINPATH/build_pips, where PIPSMAINPATH is the folder you decompress the .gz file. 
In the folder PIPSMAINPATH, we can use the following commands to configure and install PIPS:
mkdir build_pips
cd build_pips
cmake ..
make

5. The build system will install executables from three resouces: PIPS-IPM, PIPS-S and PIPS-NLP. For the usages of these executables, please follow the README.txt files in the corresponding sub-folder.

