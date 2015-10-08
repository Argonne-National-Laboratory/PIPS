/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "ps.h"
#include "stdio.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>


#include "math.h"

//#include "metis.h"

#include "graphPart.h"
#include <assert.h>


using namespace std;


extern "C" void METIS_PartGraphKway ( int *, int *, int *, int *, int *, int *,
											int *, int *, int *, int *, int * );

extern "C" void METIS_PartGraphRecursive ( int *, int *, int *, int *, int *, int *,
											int *, int *, int *, int *, int * );


static int maxnpart = 100;


graphPart::graphPart()
	:num_node(0),
	 num_balCon(1),
	 num_line(0),
	 npartsReq(0),
	 actparts(0),
	 ncut_line(0),
	 DefineBlkFlag(0),
	 xadj(NULL),
	 adjncy(NULL),
	 adjwgt(NULL),
	 vwgt(NULL),
	 id_part(NULL),
	 line_status(NULL),
  	 line_start(NULL),
  	 line_end(NULL),	
	 gname(NULL),
	 line_name(NULL),	
	 bus_name(NULL),
	 dataname(NULL),
	 vtxdist(NULL),
	 gen_status(NULL),
	 gen_bus(NULL),
	 num_gen(0),
	 local_adjncy(NULL),
	 local_xadj(NULL),
	 local_gvnTxs(0),
	 local_numTxs(-1),
	 vsize(NULL),
	 ubvec(NULL),
	 tpwgts(NULL),
	 local_vwgt(NULL),
	 local_vsize(NULL),
	 local_adjwgt(NULL),
	 local_id_part(NULL)
{}

graphPart::~graphPart()
{
  if(xadj) 		
  	delete[] xadj;
  if(adjncy) 	
  	delete[] (adjncy);
  if(id_part) 	
  	delete[] (id_part);  
  if(vwgt)		
  	delete[] (vwgt);  
  if(adjwgt)	
  	delete[] (adjwgt);  
}



/* ---------------------------------------------------------------------------
_ReadGraph: read graph file and generate 2 vectors.
---------------------------------------------------------------------------- */
void
graphPart::_ReadGraph(const char gfilename[]) 
{
  char myString[10000]; 
  char *p;
  char *save_ptr=NULL;  

  int i,j;

  FILE *gfile;
  gfile = fopen(gfilename, "r");

  /* scan First Line to get number of nodes and lines */
  fgets(myString, 10000, gfile);

  p = strtok_r(myString, " ",&save_ptr);
  num_node = atoi(p);

  printf("token:% s", p);

  
  p = strtok_r(NULL, " \n",&save_ptr);
  num_line = atoi(p);
  num_balCon = num_line;
  printf(" %s", p);

  /* get number of parts required */
  p = strtok_r(NULL, " \n",&save_ptr);
  npartsReq = atoi(p);
  printf(" %s\n", p);

  xadj = new int[num_node+1];
  adjncy = new int[2*num_line];
  id_part= new int[num_node];


  /* ------------ scan graph file defined by METIS  -------------- */
  j=0;
  xadj[0]=0;

  for(i=0;i<num_node;i++){

  	fgets(myString, 10000, gfile);
	p = strtok_r(myString, " \n",&save_ptr);

	printf("token: %s", p);
	
    while(p) {
	  adjncy[j] = atoi(p)-1;
	  j++;	  
	  p = strtok_r(NULL, " \n",&save_ptr);
	  if(p) printf(" %s", p);	  
    }
	printf("\n");
	
	xadj[i+1]=j;
  }

  assert(j==num_line*2);
  
}


// assign partID to line and generator
void 
graphPart::distributeLocalAndGlobalVar()
{

  int i,j;
  int k=0;
  int num_cutLine=0;
  int *num_localLineInEachPart;
		


  num_localLineInEachPart = (int*) malloc((actparts)*sizeof(int));
  line_status = (int*) malloc((num_line)*sizeof(int));
  if(!line_start)
  	line_start = (int*) malloc((num_line)*sizeof(int));
  if(!line_end)
	line_end = (int*) malloc((num_line)*sizeof(int));

  for(i=0;i<actparts;i++){
	num_localLineInEachPart[i]=0;
  }


  for(i=0;i<num_node;i++){
  	
	for(j=xadj[i]; j<xadj[i+1];j++){

	  if( i<adjncy[j] && id_part[i]==id_part[adjncy[j]]){
		num_localLineInEachPart[id_part[i]]++;
		num_localLineInEachPart[id_part[i]]++;
		if(line_ID){
		  line_status[line_ID[j]]=id_part[i];
		  line_start[line_ID[j]]=i;
		 line_end[line_ID[j]]=adjncy[j];
		}else{
			 line_status[k]=id_part[i];
			 line_start[k]=i;
			line_end[k]=adjncy[j];
		}
		k++;
	  }
	  else if(i<adjncy[j] && id_part[i]!=id_part[adjncy[j]]){
		num_cutLine++;
		if(line_ID){
		  line_status[line_ID[j]]=-1;
		  line_start[line_ID[j]]=i;
		  line_end[line_ID[j]]=adjncy[j];
		}else{
			 line_status[k]=-1;
			 line_start[k]=i;
			line_end[k]=adjncy[j];
		}
		k++;
	  }
	}
  }

  if(gen_status){
	for(i=0; i<num_gen;i++){
	  for(j=0; j<num_node;j++){
		if(gen_bus[i] == bus_name[j])
		  gen_status[i] = id_part[j];
	  }
	}
  }

  //assert(ncut_line == num_cutLine);
  assert(k == num_line);
  
}



void
graphPart::computeGPfromGraphFile(const char gfilename[])
{
  int retval;

  int i,j,k=0;
  int flagT;
  int *partIDFlag;
  	
  //int ncon = 1; //

//  int *vwgt = NULL, *vsize = NULL;
//  double *ubvec = NULL, *tpwgts = NULL;

  std::string graphFile (gfilename);

  graphFile.append(".graph");

  _ReadGraph(graphFile.c_str());


  num_balCon=1;

#if 0
// this is for metis version > 5.0.0
  int options[40];

  METIS_SetDefaultOptions (options);
  
/*
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;

  options[METIS_OPTION_CTYPE]   = 1;
  options[METIS_OPTION_IPTYPE]  = 4;
  options[METIS_OPTION_RTYPE]   = 1;
  options[METIS_OPTION_MINCONN] = 0;
  options[METIS_OPTION_CONTIG]  = 0;
  options[METIS_OPTION_SEED]    = -1;
  options[METIS_OPTION_NITER]   = 10;
  options[METIS_OPTION_NCUTS]   = 1;
  options[METIS_OPTION_UFACTOR] = 30;
  options[METIS_OPTION_DBGLVL]  = 0;
*/

options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // = 0
options[METIS_OPTION_CTYPE]   = METIS_CTYPE_SHEM; // = 1
options[METIS_OPTION_IPTYPE]  = METIS_IPTYPE_METISRB; // = 4
options[METIS_OPTION_RTYPE]   =  METIS_RTYPE_GREEDY; //= 1
options[METIS_OPTION_NO2HOP]  = 0;
options[METIS_OPTION_MINCONN] = 0;
options[METIS_OPTION_CONTIG]  = 0;
options[METIS_OPTION_SEED]	  = -1;
options[METIS_OPTION_NITER]   = 10;
options[METIS_OPTION_NCUTS]   = 1;
options[METIS_OPTION_UFACTOR] = 30;
options[METIS_OPTION_DBGLVL]  = 0;

tpwgts = (double *)malloc((npartsReq*num_balCon)*sizeof(double));

for (i=0; i<npartsReq; i++) {
  for (j=0; j<num_balCon; j++)
	tpwgts[i*num_balCon+j] = 1.0/(double)npartsReq;
}

vwgt = (int *)malloc((num_node*num_balCon)*sizeof(int));
for (i=0; i<num_node*num_balCon; i++) {
	vwgt[i] = 1;
}

adjwgt = (int *)malloc((2*num_line)*sizeof(int));
for (i=0; i<2*num_line; i++) {
	adjwgt[i] = 1;
}

vsize = (int *)malloc(num_node*sizeof(int));
for (i=0; i<num_node; i++) {
	vsize[i] = 1;
}


  
  printf("size of int %d, size of int * %d, size of double %d, size of double * %d\n", 
	  sizeof(int),sizeof(int *),sizeof(double),sizeof(double *));



  printf("\n Metis: Require %d parts \n", npartsReq);


  retval = METIS_PartGraphKway(&num_node, &num_balCon, xadj, adjncy, 
	  vwgt, vsize, adjwgt, &npartsReq, tpwgts, ubvec, options, &ncut_line, id_part);


  // check how many parts we can have after METIS
  partIDFlag = (int*) malloc(npartsReq*sizeof(int));
  for(i=0; i<npartsReq;i++){
	partIDFlag[i] = -1;
  }  

  for(i=0; i<num_node;i++){
  	flagT=0;
  	for(j=0; j<k; j++){
	  if(id_part[i]==partIDFlag[j]){
		flagT=1;
		break;
	  }
	}
	if(flagT==0){
	  partIDFlag[k++] = id_part[i];
	  actparts++;
	}
  }

  switch (retval)
  {
	case 1://METIS_OK:
	  printf ("\nMETIS_OK\n");
	  printf("\n Metis: Actual %d parts \n", actparts);

	  distributeLocalAndGlobalVar();

	  break;
	case -2://METIS_ERROR_INPUT:
	  printf ("\nMETIS_ERROR_INPUT\n");
	  exit(1);
	case -3://METIS_ERROR_MEMORY:
	  printf ("\nMETIS_ERROR_MEMORY\n");
	  exit(1);
	case -4://METIS_ERROR:
	  printf ("\nMETIS_ERROR\n");
	  exit(1);
  }
#else 
  int options[5];
  options[0]=0;

  int useWgt   = 0;
  int idxFrom1 = 0;

//  vwgt = (int *)malloc((num_node*num_balCon)*sizeof(int));
//  for (i=0; i<num_node*num_balCon; i++) {
//	vwgt[i] = 1;
//  }
  
//  adjwgt = (int *)malloc((2*num_line)*sizeof(int));
//  for (i=0; i<2*num_line; i++) {
//	adjwgt[i] = 1;
//  }

  actparts = npartsReq;
  
  if(npartsReq<=8)
    METIS_PartGraphRecursive(&num_node, xadj, adjncy, NULL, NULL, &useWgt, &idxFrom1, 
	    				&npartsReq, options, &ncut_line, id_part);
  else
  	METIS_PartGraphKway(&num_node, xadj, adjncy, NULL, NULL, &useWgt, &idxFrom1, 
	    				&npartsReq, options, &ncut_line, id_part);
  
  if(npartsReq)
  	printf("--- Generate %d parts by metis. \n", npartsReq);

  if(actparts != npartsReq){
    printf("--- Generate %d parts by metis. (Require %d parts) \n", npartsReq, actparts);
  	actparts = npartsReq;
  }
    
#endif 


//  delete (graphFile);
}






