/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "ps.h"
#include "stdio.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>

#include <cassert>


using namespace std;

DCBUS::DCBUS()
	: bus_i(0), va(0), ngen(0),nconnlines(0),loads(0.),
	  startloc(0), startloc_Part(-1), partID(-1),haveCutLine(false)
{}

void DCBUS::DCBUSGetNGen(int *ngen)
{
  *ngen = this->ngen; 
}

void DCBUS::DCBUSGetGen(int gen_num,DCGEN **gen)
{
  if(gen_num < 0 || gen_num > this->ngen-1){
  	printf("Ngen at bus = %d, generator number requested %d",this->ngen,gen_num);
	exit(0);
  }
  *gen = this->gens[gen_num];
}

void DCBUS::DCBUSGetLoad(double *load)
{ 
  *load = this->loads;
}

void DCBUS::DCBUSGetVariableLocation(int *loc)
{
  *loc = this->startloc;
}

void DCBUS::DCBUSGetVariableLocation_Part(int *loc)
{
  *loc = this->startloc_Part;
}

void DCBUS::DCBUSGetSupportingLines(int *nsupplines, DCLINE ***supplines)
{
  *nsupplines = this->nconnlines;
  *supplines  = this->connlines;
}

DCLINE::DCLINE()
	: fbus(0), tbus(0), resistance(0.),reactance(0.),tapratio(0.),
	  partID(-1)
{}


void DCLINE::DCLINEGetConnectedBuses(DCBUS ***connbuses)
{
  *connbuses = this->connbuses;
}

DCGEN::DCGEN()
	: bus_i(0), pg(0.), pt(0.),pb(0.),
	  cost_alpha(0.), cost_beta(0.), cost_gamma(0.),partID(-1)
{}



DCPS::DCPS()
	: MVAbase(100.0), Nbus(0), Ngen(0),Nbranch(0),
	  refBusID(0),
	  bus(NULL), gen(NULL), line(NULL),
	  setupcalled(false), 
	  Nparts(1),Ncuts(0),num_LineInEachPart(NULL),num_CutInEachPart(NULL),num_BusInEachPart(NULL),BusInEachPart(NULL)
{}

DCPS::~DCPS()
{
  if(setupcalled){
	free(bus);
	free(line);
	free(gen);
  }

  if(num_LineInEachPart){
	free(num_LineInEachPart);
  }	
  if(num_BusInEachPart){
	free(num_BusInEachPart);
  }	  
  if(num_CutInEachPart){
	free(num_CutInEachPart);
  }	 
  

  if (BusInEachPart){
    for(int i=0;i<Nparts;i++){
	  delete [] BusInEachPart[i];
    }
    delete [] BusInEachPart;
  }
}

void DCPS::PSGetNumGenerators(int *Ngen_)
{
  if(Ngen_) *Ngen_ = Ngen;
}


void DCPS::PSReadDCData(const char netfile[])
{
  FILE           *fp;
  DCBUS          *Bus;
  DCGEN          *Gen;
  DCLINE         *Branch;
  int       line_counter=0;
  int       bus_start_line=-1,bus_end_line=-1; /* xx_end_line points to the next line after the record ends */
  int       gen_start_line=-1,gen_end_line=-1;
  int       br_start_line=-1,br_end_line=-1;
  int       gencost_start_line=-1,gencost_end_line=-1;
  char      fileline[PS_MAXLINE];
  int       loadi=0,geni=0,bri=0,busi=0,gencosti=0,i;
  int       extbusnum,bustype_i;
  double    Pd,Qd;
  int       intbusnum;

  fp = fopen(netfile,"r");
  /* Check for valid file */
  if (fp == NULL) {
    printf("Cannot open file %s",netfile);
	exit(0);
  }

  while(fgets(fileline,PS_MAXLINE,fp) != NULL) {
    if(strstr(fileline,"param: 	Bus") != NULL)    bus_start_line = line_counter+1; /* Bus data starts from next line */
    if(strstr(fileline,"param: 	Generator") != NULL && gen_start_line == -1)    gen_start_line = line_counter+1; /* Generator data starts from next line */
    if(strstr(fileline,"param: 	Line") != NULL) br_start_line = line_counter+1; /* Branch data starts from next line */
    if(strstr(fileline,";") != NULL) {
      if (bus_start_line != -1 && bus_end_line == -1) bus_end_line = line_counter+1;
      if (gen_start_line != -1 && gen_end_line == -1) gen_end_line = line_counter+1;
      if (br_start_line  != -1 && br_end_line == -1)   br_end_line = line_counter+1;
    }
    line_counter++;
  }
  fclose(fp);

  Nbus        = bus_end_line - bus_start_line ;
  Ngen        = gen_end_line - gen_start_line ;
  Nbranch     = br_end_line  - br_start_line  ;

  printf("nbus = %d, ngen = %d, nbranch = %d\n",Nbus,Ngen,Nbranch);

  bus = (DCBUS*)malloc(Nbus*sizeof(DCBUS));
  gen = (DCGEN*)malloc(Ngen*sizeof(DCGEN));
  line = (DCLINE*)malloc(Nbranch*sizeof(DCLINE));
  
  Bus = bus; Gen = gen; Branch = line;

  for(i=0; i < Nbus; i++) {
    bus[i].ngen = 0;
  }

  /* Setting sbase to 100 */
  MVAbase = 100.0;


  /* Start to read data*/  
  int tempID=0;
  int tempInt=0;
  int tempDouble=0;
  
  fp = fopen(netfile,"r");
  for(i=0;i<line_counter;i++) {
    fgets(fileline,PS_MAXLINE,fp);

	/* Read bus data */
    if((i >= bus_start_line) && (i < bus_end_line)) {
	  tempID=0;
      sscanf(fileline,"%*s %d %lf %*lf %*lf %*lf",		\
	     &tempInt,&Bus[busi].loads);

	  Bus[busi].bus_i      = busi;
	  if(tempInt==1) refBusID = busi;
	  
      busi++;
    }

    /* Read generator data */
    if(i >= gen_start_line && i < gen_end_line) {
      sscanf(fileline,"%*d %d %lf %lf %*lf %*lf %lf %lf %lf", \
	  	 &Gen[geni].bus_i, \
	     &Gen[geni].pt,&Gen[geni].pb,&Gen[geni].cost_gamma,&Gen[geni].cost_beta,&Gen[geni].cost_alpha);

      geni++;
    }
	
    /* Read line data */
    if(i >= br_start_line && i < br_end_line) {
      sscanf(fileline,"%*d %d %d %lf %lf %*lf %lf %*lf %*d", \
	  	 &Branch[bri].fbus, &Branch[bri].tbus, \
	     &Branch[bri].resistance, &Branch[bri].reactance,		 \
	     &Branch[bri].flowLimit);
	  
      Branch[bri].tapratio = 1.0;

      bri++;
    }
  }
  fclose(fp);

  SetUpPS();
  
}

void DCPS::SetUpPS(){
   /* Set up
     (a) connectivity information for lines and buses 
  */
  int i;
   
  /* Set bus  */
  for( i=0; i < Nbranch; i++ ){  
	bus[line[i].fbus].connlines[ bus[line[i].fbus].nconnlines ] = &line[i];
	bus[line[i].tbus].connlines[ bus[line[i].tbus].nconnlines ] = &line[i];
  
	bus[line[i].fbus].nconnlines++; 	
	bus[line[i].tbus].nconnlines++;

    if (bus[line[i].fbus].nconnlines > NLOAD_AT_BUS_MAX) {
	  printf("Exceeded maximum number of lines allowed at bus");
	  exit(0);
    }	
  }
  for( i=0; i < Ngen; i++ ){  
	bus[gen[i].bus_i].gens[ bus[gen[i].bus_i].ngen ] = &gen[i];
  
	bus[gen[i].bus_i].ngen++; 	

    if (bus[gen[i].bus_i].ngen > NGEN_AT_BUS_MAX) {
	  printf("Exceeded maximum number of generators allowed at bus");
	  exit(0);
    }
  }

  /* Set lines  */
  for( i=0; i < Nbranch; i++ ){  
	line[i].connbuses[0] = &bus[line[i].fbus];
	line[i].connbuses[1] = &bus[line[i].tbus];
  }

  int xloc=0;
  for( i=0; i < Nbus; i++ ){  
    bus[i].startloc = xloc;
	xloc += 1 + bus[i].ngen;
  }

  setupcalled = true;

}



void DCPS::writeGraphFile(const char gfilename[], const int numPart )
{
  int nLink;
  int fBus, tBus;

  int   mloc=0;

  DCBUS          *busf,*bust;
  DCBUS **connbuses;
  
  DCLINE   *linePt;  
  DCLINE   **connlines;

  std::string graphFile (gfilename);
  graphFile.append(".graph");

  ofstream myfile (graphFile.c_str());
  if (myfile.is_open())
  {
    myfile << Nbus << " " << Nbranch << " " << numPart << " \n";

	for(int i=0;i<Nbus;i++){
	  bus[i].DCBUSGetSupportingLines(&nLink,&connlines);
	  for(int k=0; k < nLink; k++) {
		linePt = connlines[k];
		linePt->DCLINEGetConnectedBuses(&connbuses);
		busf = connbuses[0];
		bust = connbuses[1];

		if(bus[i].bus_i == busf->bus_i) {
		  myfile << bust->bus_i+1 << " ";
		} else {
		  myfile << busf->bus_i+1 << " ";
		}
	  }
	  myfile << "\n";
	}

    myfile.close();
  }
  else cout << "Unable to open file";

  cout << "Generate " << graphFile << " file!\n";


}


// assign partID to buses and lines
void 
DCPS::setPartitioning(const int *partID_, const int actparts_, const int num_cuts)
{
  int nLink;

  int fBus, tBus;
  DCBUS  *busf,*bust;
  DCBUS **connbuses;
  
  DCLINE   *linePt;  
  DCLINE   **connlines;

  Nparts = actparts_;
	
  num_LineInEachPart 	= (int*) malloc((Nparts)*sizeof(int));
  num_BusInEachPart 	= (int*) malloc((Nparts)*sizeof(int));  
  num_CutInEachPart 	= (int*) malloc((Nparts)*sizeof(int));
  
  for(int i=0;i<Nparts;i++){
	num_LineInEachPart[i]=0;
	num_BusInEachPart[i]=0;
	num_CutInEachPart[i]=0;	
  }

  // assign partID to buses and lines  
  for(int i=0;i<Nbus;i++){
	bus[i].DCBUSGetSupportingLines(&nLink,&connlines);
	bus[i].partID = partID_[i];
	num_BusInEachPart[partID_[i]]++;
	
	for(int k=0; k < nLink; k++) {
	  linePt = connlines[k];
	  linePt->DCLINEGetConnectedBuses(&connbuses);
	  busf = connbuses[0]; fBus = busf->bus_i;
	  bust = connbuses[1]; tBus = bust->bus_i;

	  if( partID_[fBus] == partID_[tBus] ){
		linePt->partID = partID_[i];
		num_LineInEachPart[partID_[i]]++;
	  }else{
		Ncuts++;
		num_CutInEachPart[partID_[i]]++;	
		linePt->partID = -1;
		busf->haveCutLine = true;
		bust->haveCutLine = true;		
	  }
	}
  }

  Ncuts = Ncuts/2;
  assert(Ncuts==num_cuts);

  int *xloc = new int[Nparts];
  for( int i=0; i < Nparts; i++ ){	
	xloc[i] = 0;
  }

  for( int i=0; i < Nbus; i++ ){	
  	int partID = bus[i].partID;
    bus[i].startloc_Part = xloc[partID];
    xloc[partID] += 1 + bus[i].ngen;
  }

  for( int i=0; i < Ngen; i++ ){	
  	gen[i].partID = partID_[gen[i].bus_i];
  }

  BusInEachPart 		= new int*[Nparts];
  for(int i=0;i<Nparts;i++){
	BusInEachPart[i]	= new int[num_BusInEachPart[i]];
  }
  int *findbus_part = new int[Nparts];
  for(int i=0;i<Nparts;i++) findbus_part[i]=0;
  
  for(int i=0;i<Nbus;i++){
  	int partID_temp = bus[i].partID;
  	BusInEachPart[partID_temp][findbus_part[partID_temp]++] = i;
  }
  for(int i=0;i<Nparts;i++){
	assert( findbus_part[i] == num_BusInEachPart[i]); 
  } 

  delete [] findbus_part;

  delete [] xloc;
}

