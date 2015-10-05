#include "stdio.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>

#include "ps.h"
#include "graphPart.h"


int main (int argc, char *argv[])
{

 int i;

 
 FILE *file=NULL;
 int   arg;
 char *gname;
 graphPart graphP;


 DCPS dcpowersys;
 if (argc <= 2) {
   fprintf(stderr, "%s:  usage: -f <datafilename>\n", argv[0]);
   exit(1);
 }

 for (arg=1; arg<argc && argv[arg] != NULL; arg++) {
   if (strcmp(argv[arg], "-f") == 0) {
	 if (arg + 1 >= argc) {
       fprintf(stderr, "%s:  usage: -f <datafilename>\n", argv[0]);
	   exit(1);
	 }
 
	 file = fopen(argv[arg+1], "r");
	 if (file== NULL) {
       fprintf(stderr,"%s: could not open file %s\n", argv[0], argv[arg+1]);
	   exit(1);
	 }
		   /* store the problem name */
	 gname = strdup(argv[arg + 1]);
	 break;
   }
 }


  dcpowersys.PSReadDCData(gname);
  dcpowersys.writeGraphFile(gname,2);

  graphP.computeGPfromGraphFile(gname);
  dcpowersys.setPartitioning(graphP.id_part,graphP.actparts,graphP.ncut_line);


  free(gname);
  
}

