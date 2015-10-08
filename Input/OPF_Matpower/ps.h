/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef __PS_H
#define __PS_H

#include <constants.h>
/* PS.H is the header file containing the forward declarations for the PS data structs.
   It also contains its API declarations 
*/
  
class DCBUS;
class DCGEN;
class DCLINE;


/* Bus data struct */
class DCBUS{
 public:
  DCBUS();
  ~DCBUS(){};
  
  int      bus_i; /* Integer bus number */
  char     i[20]; /* Bus Number */
  char     name[20]; /* Bus Name */
  double   va; /* Bus voltage phase angle */


  /* Info of connected edges to this bus */
  int      nconnlines;
  DCLINE   *connlines[MAXCONNLINES];

  int      ngen; /* Number of generators incident at this bus */
  DCGEN    *gens[NGEN_AT_BUS_MAX];
  double   loads;

  /* Starting location for the variables for this bus in the application residual "local" array.
     The local array also contains ghosted values. startloc is typically used to access values in
     the local vector, while startlocglob is used for setting the values in the matrix */
  int startloc; 

  int startloc_Part;

  int partID;

  bool haveCutLine;

  void DCBUSGetNGen(int*);
  void DCBUSGetGen(int,DCGEN**);
  void DCBUSGetLoad(double*);
  void DCBUSGetVariableLocation(int*);
  void DCBUSGetVariableLocation_Part(int*);
  void DCBUSGetSupportingLines(int*, DCLINE***);
  
};


/* Generator data struct */
class DCGEN{
 public:
 	
  DCGEN();
  ~DCGEN(){};

  int      bus_i;
  char     i[20]; /* Bus Number or extended bus name*/
  char     id[2]; /* Generator identifier, in case of multiple generators at same bus. 1 by default */
  
  double   pg; /* Generator active power output */
  
  double   pt; /* Gen max active power output: MW */
  double   pb; /* Gen min active power output: MW */

  /* Quadratic cost function is alpha*Pg^2 + beta*Pg + gamma */
  double   cost_gamma;
  double   cost_beta;
  double   cost_alpha;

  
  int partID;
  
};

/* Line data structure */
class DCLINE{
	public:

   DCLINE();
   ~DCLINE(){};
  
  int      fbus;
  int      tbus;
  char     i[20]; /* Bus Number or extended bus name*/
  char     j[20]; /* Bus Number or extended bus name*/

  double   resistance; /* Branch resistance: pu */
  double   reactance; /* Branch reactance: pu */

  double   flowLimit; /* Branch flow limit*/
  
  double   tapratio;

  /* From and to buses */
  DCBUS  *connbuses[2];

  int partID;


  void DCLINEGetConnectedBuses( DCBUS***);
  
};

/* Power system data structure */
class DCPS{
public:
  double MVAbase; /* System base MVA */
  int    Nbus,Ngen,Nbranch; /* global # of buses,gens,branches, and loads (includes elements which are
                                          out of service */
//  int    nbus,ngen,nbranch; /* local # of buses,gens,branches,and loads */
  DCBUS   *bus;
  DCGEN   *gen;
  DCLINE  *line;

  /* keys for components */
  int compkey[10];

  bool setupcalled; /* Is setup called on PS? */

  int refBusID;



  DCPS();
  ~DCPS();
  
//  void PSReadMatPowerData(const char[]){};
  
  void PSReadDCData(const char[]);

  void PSGetNumGenerators(int*);

  void writeGraphFile(const char[], const int);
  
  void setPartitioning(const int *partID_, const int actparts_,const int num_cuts);



  int Nparts;
  int Ncuts;
  int *num_LineInEachPart;

  int *num_BusInEachPart;

  int *num_CutInEachPart;



  int **BusInEachPart;

private:
  void SetUpPS();
  
};


#endif

