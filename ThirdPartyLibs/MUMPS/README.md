
1) Download MUMPS from http://mumps.enseeiht.fr/index.php?page=dwnld

2) Build parallel libraries with double ('d') arithmetic

3) Copy static libraries libdmumps.a, libmumps_common.a, and libpord.a into 'ThirdPartyLibs/MUMPS/lib' directory; copy include directory into 'ThirdPartyLibs/MUMPS/include/' 

4) Build ParMetis static library and copy to 'ThirdPartyLibs/METIS/src/lib/libparmetis.a'

5) Link parmetis and mumps against a metis library that exports the symbol metis_setdefaultoptions_ (not the case for metis recieved via wgetMETIS.sh) rather use METIS shipped with parmetis
