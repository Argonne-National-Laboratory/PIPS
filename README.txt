1. Go to the following folders and run the script wgetXXX.sh
ThirdPartyLib/ASL  
ThirdPartyLib/CBC 
ThirdPartyLib/ConicBundle   
ThirdPartyLib/METIS
ThirdPartyLib/UMFPACK
For an example, use command sh wgetASL.sh in the folder ThirdPartyLib/ASL  

2. Download MA27 and MA57 from HSL and put the source code ma27.f and ma57.f in the folder ThirdPartyLib/MA27/src/ and ThirdPartyLib/MA57/src, respectively. Then you can use the Makefile in the folder to install it. (See ThirdPartyLob/MA27/README.txt for more details.)

3. Use following command to install pips
cmake PATH_TO_PIPS
where PATH_TO_PIPS is the main folder where you decompressed pipsnlp.tar.gz
