/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/
 

#include "dcopflow.hpp"
#include <stddef.h>
#include <assert.h>

DCOPFLOW::DCOPFLOW(const DCPS* ps_)
	: ps(ps_), 
	  obj(0.0),
	  gradobj(NULL), 
	  Xl(NULL), Xu(NULL), Gl(NULL), Gu(NULL),
	  Nvar(-1), Nconeq(-1), Nconineq(-1), Ncon(-1),
	  n(-1),m(-1),nnz_jac_g(-1),Nparts(0),
	  setupcalled(false), setupcalled_part(false),
	  lambda_g(NULL), lambda_xl(NULL), lambda_xu(0)
{
  if(ps->Nparts==1)
  	DCOPFLOWSetUp();
  else
  	DCOPFLOWSetUp_Partition();
}

DCOPFLOW::~DCOPFLOW()
{
  if(setupcalled){
	delete(Xl);
	delete(Xu);
	delete(Gl);
	delete(Gu);
	delete(gradobj);
	delete(lambda_g);
	delete(lambda_xl);
	delete(lambda_xu);	
  }	
  if(setupcalled_part){
    delete(busMap_AllTo1st);
    delete(Nvar_2nd);  
    delete(Ncon_2nd); 
    delete(Nconeq_2nd); 
    delete(Nconineq_2nd); 
	delete(numDummyVar);
	delete(numDummyCon);
  }
}


/*
  OPFLOWSetUp - Sets up a power flow application object

  Input Parameters:
. opflow - the OPFLOW object

  Notes:
  This routine sets up the OPFLOW object and the underlying PS object. It
  also distributes the PS object when used in parallel.
*/
void
DCOPFLOW::DCOPFLOWSetUp()
{

  Nvar		= ps->Nbus + ps->Ngen;
  Nconeq 	= ps->Nbus + 1;
  Nconineq 	= ps->Nbranch;
  Ncon 		= Nconeq + Nconineq;

  /* Create the solution vector and upper and lower bounds */
  Xl = new double[Nvar];
  Xu = new double[Nvar];
  
  /* Create the constraint function vector and its bounds */
  Gl = new double[Ncon];
  Gu = new double[Ncon];

  /* Create the gradient vector */
  gradobj = new double[Nvar];

  /* Create vectors for multipliers */
  lambda_g  = new double[Ncon];
  lambda_xl = new double[Nvar];
  lambda_xu = new double[Nvar];

  /* Constraint jacobian */
  GetJacNNZ();

  /* Set var/con bounds */
  SetVariableandConstraintBounds();

  setupcalled = true;
}



/*
  OPFLOWSetVariableandConstraintBounds - Sets the bounds on variables and constraints

  Input Parameters:
. opflow - the OPFLOW object

  Output Parameters:
+ Xl     - vector of lower bound on variables
. Xu     - vector of upper bound on variables
. Gl     - vector of lower bound on constraints
. Gu     - vector of lower bound on constraints
*/
void 
DCOPFLOW::SetVariableandConstraintBounds()
{
  DCLINE         *line;
  DCBUS          *bus;
  int       loc,gloc=0;


  for(int i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];

    /* Line flow inequality constraints */
    Gl[gloc] = -(line->flowLimit/ps->MVAbase);
    Gu[gloc] =  (line->flowLimit/ps->MVAbase);
       
    gloc += 1;
  }

  for(int i=0; i < ps->Nbus; i++) {
    bus = &ps->bus[i];
    bus->DCBUSGetVariableLocation(&loc);

    /* Bounds on voltage angles and bounds on real power mismatch equality constraints */
    Xl[loc] = -M_PI; Xu[loc] = M_PI;
    Gl[gloc] = 0.0;  Gu[gloc] = 0.0;
    
    for(int k=0; k < bus->ngen; k++) {
      DCGEN *gen;
      bus->DCBUSGetGen(k,&gen);
      loc += 1;
      Xl[loc] = gen->pb/ps->MVAbase; Xu[loc] = gen->pt/ps->MVAbase;
    }
	gloc += 1;
  }
}


/*
  OPFLOWObjGradientFunction - The linear part of the objecive gradient 
  
  Input Parameters:
+ opflow - the OPFLOW object
. X      - the current iterate

  Output Parameters:
. grad -  The linear part of the objecive gradient 
*/
void 
DCOPFLOW::ObjGradient_Lin(double *obj_coef)
{
  DCBUS     *bus;
  DCGEN     *gen;
  int       loc;

  for(int i=0;i<Nvar;i++) obj_coef[i]=0;

  for(int i=0; i < ps->Nbus; i++) {
    bus = &ps->bus[i];
    for(int k=0; k < bus->ngen; k++) {
      bus->DCBUSGetGen(k,&gen);
      obj_coef[loc] = ps->MVAbase*( gen->cost_beta );
    }
  }
}

/*
  OPFLOWObjGradientFunction - The quad part of the objecive gradient 
  
  Input Parameters:
+ opflow - the OPFLOW object
. X      - the current iterate

  Output Parameters:  
  The row and column locations should be such that J(row[i],col[i]) = val
*/
void 
DCOPFLOW::ObjGradient_Quad(int *irow, int *jcol, double *obj_quad)
{
  DCBUS     *bus;
  DCGEN     *gen;
  int       loc;
  int		ctr=0;

  for(int i=0; i < ps->Nbus; i++) {
    bus = &ps->bus[i];
	bus->DCBUSGetVariableLocation(&loc);
    for(int k=0; k < bus->ngen; k++) {
	  loc += 1;
	  irow[ctr] = loc; 
	  jcol[ctr] = loc; 
      bus->DCBUSGetGen(k,&gen);
      obj_quad[ctr] = 2.0 * gen->cost_alpha * ps->MVAbase * ps->MVAbase;
	  ctr += 1;
    }
  }
}


/*
  DCOPFLOWGetConstraintJacobianNonzeros - Gets the number of nonzeros in the constraint jacobian matrix

  Input Paramereters:
. opflow - the optimal power flow application object

  Output Parameters:
. nnz - number of nonzeros in the constraint Jacobian

  Notes:
  DCOPFLOWGetConstraintJacobianNonzeros() must be called before calling this routine.
*/
int 
DCOPFLOW::GetJacNNZ()
{
  int nnz = 0;
  
  DCBUS          *bus;
  DCLINE         *line;
  
  int       nconnlines;
  DCLINE   **connlines;

  if(nnz_jac_g==-1){
    // KVL constraint: flow_ij = -V^2/r * (phase_i - phase_j)
    for(int i=0; i < ps->Nbranch; i++) nnz += 2;

    // KCL constraint:  \sum pgen  +  \sum flow_in -  \sum flow_out = d
    for(int i=0; i < ps->Nbus; i++) {
      bus = &ps->bus[i];

      nnz += bus->ngen;

	  bus->DCBUSGetSupportingLines(&nconnlines,&connlines);

	  nnz += 2*nconnlines;
    }
	// Ref bus constraint: phase angle = 0 for the reference bus
    nnz += 1;
	
	nnz_jac_g = nnz;
  }else{
	nnz = nnz_jac_g;
  }  
  return nnz;
}






/*
  OPFLOWSetConstraintJacobianLocations - Sets the row and column nonzero locations for the
              constraint Jacobian

  Input Parameters:
+ opflow - the OPFLOW object

  Output Parameters:
+ row - the row locations
- col - the col locations

  Notes:
   The row and column locations should be such that J(row[i],col[i]) = val
*/
void
DCOPFLOW::SetJacLocations(int *row, int *col)
{
  int nnz = 0;
  
  DCBUS          *bus;
  DCLINE         *line;
  
  int       nconnlines;
  DCLINE   **connlines;

  DCBUS	  *busf,*bust;
  DCBUS	  **connbuses;  
  
  int   gloc=0,xloc,xlocf,xloct;
  int 	ctr=0;

  
  // KVL constraint: flow_ij = -V^2/r * (phase_i - phase_j)
  for(int i=0; i < ps->Nbranch; i++){
  	line = &ps->line[i];

	line->DCLINEGetConnectedBuses(&connbuses);		
	busf = connbuses[0];
	bust = connbuses[1];	
    busf->DCBUSGetVariableLocation(&xlocf);
    bust->DCBUSGetVariableLocation(&xloct);

    row[ctr] = row[ctr+1] = gloc;
    col[ctr] = xlocf; col[ctr+1] = xloct;
    ctr  += 2;

    gloc += 1;
  }

  // KCL constraint:  \sum pgen  +  \sum flow = d
  for(int i=0; i < ps->Nbus; i++) {
    bus = &ps->bus[i];
    bus->DCBUSGetVariableLocation(&xloc);

    for(int k=0; k < bus->ngen; k++) {
      xloc += 1;

      row[ctr] = gloc;     
	  col[ctr] = xloc;
      ctr += 1;
    }
	
	bus->DCBUSGetSupportingLines(&nconnlines,&connlines);

    for(int k=0; k < nconnlines;k++) {
	  line = connlines[k];
	  line->DCLINEGetConnectedBuses(&connbuses);		
	  busf = connbuses[0];
	  bust = connbuses[1];
	  busf->DCBUSGetVariableLocation(&xlocf);
      busf->DCBUSGetVariableLocation(&xloct);

	  row[ctr] = row[ctr+1] = gloc;
	  col[ctr] = xlocf; col[ctr+1] = xloct;
	  ctr += 2;
    }
	gloc += 1;
  }
  
  // Ref bus constraint: phase angle = 0 for the reference bus
  bus = &ps->bus[refBusID];
  bus->DCBUSGetVariableLocation(&xloc);
  row[ctr] = gloc;
  col[ctr] = xloc;
  ctr  += 1;
  gloc += 1;

  assert(gloc == Ncon);
}



/*
  OPFLOWSetConstraintJacobianValues - Sets the nonzero values for the
              constraint Jacobian

  Input Parameters:
+ opflow - the OPFLOW object
- X      - the current iterate


  Output Parameters:
. values - the nonzero values in the constraint jacobian

  Notes:
   The row and column locations should be such that J(row[i],col[i]) = val
*/
void 
DCOPFLOW::SetJacValues(double *values)
{
  int nnz = 0;
  
  DCBUS          *bus;
  DCLINE         *line;
  
  int       nconnlines;
  DCLINE   **connlines;

  DCBUS	  *busf,*bust;
  DCBUS	  **connbuses;  
  
  int 	ctr=0;

  
  // KVL constraint: flow_ij = -V^2/r * (phase_i - phase_j)
  for(int i=0; i < ps->Nbranch; i++){
  	line = &ps->line[i];
	values[ctr]   =  1.0/line->reactance;
	values[ctr+1] = -1.0/line->reactance;
    ctr  += 2;
  }

  // KCL constraint:  \sum pgen  +  \sum flow = d
  for(int i=0; i < ps->Nbus; i++) {
    bus = &ps->bus[i];

    for(int k=0; k < bus->ngen; k++) {
	  values[ctr] = 1.0;
      ctr += 1;
    }
	
	bus->DCBUSGetSupportingLines(&nconnlines,&connlines);

    for(int k=0; k < nconnlines;k++) {
	  line = connlines[k];
	  busf = connbuses[0];
	  bust = connbuses[1];

	  if(bus == busf) {	  
	    values[ctr]	  = -1.0/line->reactance;
	    values[ctr+1] =  1.0/line->reactance;
	  }else{
		values[ctr]	  =  1.0/line->reactance;
		values[ctr+1] = -1.0/line->reactance;
	  }
	  ctr  += 2;
    }
  }
}




void
DCOPFLOW::DCOPFLOWSetUp_Partition()
{
  DCBUS		   *bus;
  DCLINE       *line;
  
  int       nconnlines;
  DCLINE   **connlines;
  DCBUS	  *busf,*bust;
  DCBUS	  **connbuses;  
  int   busfID,bustID;

  int Nvar_dummy=0;
  Nparts = ps->Nparts;

  Nvar_2nd 			= new int[Nparts];
  Ncon_2nd 			= new int[Nparts];  
  Nconeq_2nd 		= new int[Nparts];
  Nconineq_2nd 		= new int[Nparts];   

  numDummyVar		= new int[Nparts];
  numDummyCon		= new int[Nparts];
  nnz_jac_g_2nd 	= new int[Nparts];
  nnz_jac_g_Link 	= new int[Nparts];    

  busMap_AllTo1st 	= new int[ps->Nbus]; 
  varMap_AllTo2nd 	= new int[ps->Nbus];

  nnz_jac_g_1st = -1;
  
  for(int i=0; i < Nparts; i++) {
	Nvar_2nd[i] 	= 0;
	Ncon_2nd[i] 	= 0;
	Nconeq_2nd[i] 	= 0;
	Nconineq_2nd[i] = 0;
	numDummyVar[i]	= 0;
	numDummyCon[i]	= 0;
	nnz_jac_g_Link[i]  = -1;
	nnz_jac_g_2nd[i]   = -1;
  }  

  for(int i=0; i < ps->Nbus; i++) {
	busMap_AllTo1st[i]=-1;
	varMap_AllTo2nd[i]=-1;
  }  

  // define num of constraints on the cuts
  Nvar_1st  	= 0;
  Nconeq_1st 	= 0;
  if(ps->bus[refBusID].haveCutLine) Nconeq_1st++;    
  Nconineq_1st 	= ps->Ncuts;
  Ncon_1st = Nconineq_1st + Nconeq_1st;

  bool **usedBusID;
  usedBusID = new bool*[Nparts];
  for(int i=0; i<Nparts; i++) usedBusID[i] = new bool[ps->Nbus];
  for(int i=0; i<Nparts; i++) for(int j=0; j<ps->Nbus; j++) usedBusID[i][j] = false;

  // num of var in each partition and the couping part  
  for(int i=0; i < ps->Nbus; i++) {
    bus = &ps->bus[i];
	bus->startloc_Part = Nvar_2nd[bus->partID];

	// local phase angle variable  (for bus with not cuts )
	Nvar_2nd[bus->partID]++;

	// local  node balance constraint. (= #node in each partition)
	Nconeq_2nd[bus->partID]++;

	// ref bus constraint
	if(refBusID == bus->bus_i && !(ps->bus[refBusID].haveCutLine) ) Nconeq_2nd[bus->partID]++;	

	bus->DCBUSGetSupportingLines(&nconnlines,&connlines);
	// find cut lines linked to this bus
	if(bus->haveCutLine){
	  numDummyCon[bus->partID]++;
	  if(busMap_AllTo1st[i]==-1){
		busMap_AllTo1st[i] = Nvar_1st;	  
		Nvar_1st++;
	  }		
	  for(int k=0; k < nconnlines;k++) {
		line = connlines[k];
		if(line->partID == -1){
		  line->DCLINEGetConnectedBuses(&connbuses);		
		  busf = connbuses[0];  busfID = busf->bus_i;
		  bust = connbuses[1];  bustID = bust->bus_i;					
		  // cut line, one node is in this partition.  add phase angle for dummy bus
		  if(busf->partID == bus->partID){
			// from bus is in this part
			if(!usedBusID[bus->partID][bustID]){
			  numDummyVar[bus->partID]++;
			  numDummyCon[bus->partID]++;
			  usedBusID[bus->partID][bustID] = true;
			} 	  
		  }else{
			// to bus is in this part
			if(!usedBusID[bus->partID][busfID]){
			  numDummyVar[bus->partID]++;
			  numDummyCon[bus->partID]++;
			  usedBusID[bus->partID][busfID] = true;
			} 		  
		  }
		}
      }
	}
	// local power generation variable 
	Nvar_2nd[bus->partID] += bus->ngen; 
  }

  for(int i=0; i < ps->Nbranch; i++) {
	line = &ps->line[i];	
	// local line flow constraints	
	if(line->partID!=-1) 
	  Nconineq_2nd[line->partID] += 1;
  }  
  
  // add dummy var/con in the end
  for(int i=0;i<Nparts; i++){  
    Nvar_2nd[i]  += numDummyVar[i]; 
	Nconeq_2nd[i]+= numDummyCon[i]; 	
	Ncon_2nd[i] = Nconeq_2nd[i] + Nconineq_2nd[i];
  }

  // compute global number of var/cons
  Nvar		= Nvar_1st;   
  Nconeq 	= Nconeq_1st;
  Nconineq 	= Nconineq_1st;  
  for(int i=0; i < Nparts; i++) {
	Nvar 	 += Nvar_2nd[i];
	Nconeq	 += Nconeq_2nd[i];
	Nconineq += Nconineq_2nd[i];	
  }  
  Ncon		  = Nconeq + Nconineq;

  assert(Nconineq = ps->Nbranch);

  
  setupcalled_part = true;

  for(int i=0; i<Nparts; i++) 
  	delete [] usedBusID[i];
  	
  delete [] usedBusID; 
}


void
DCOPFLOW::VarAndConBounds_1st_Partition(double *xl,double *xu,double *gl,double *gu)
{
  DCLINE         *line;
  DCBUS          *bus;
  int       loc,gloc=0;
  
  // Line flow inequality constraints for cuts
  for(int i=0; i < ps->Nbranch; i++) {
  	line = &ps->line[i];	
	if(line->partID==-1) {
	  gl[gloc] = -(line->flowLimit/ps->MVAbase);
	  gu[gloc] =  (line->flowLimit/ps->MVAbase);	  
	  gloc += 1;
	}       
  }
  assert(gloc == ps->Ncuts);

  // ref bus constraint
  if(ps->bus[refBusID].haveCutLine) {
	gl[gloc] = 0.0;  
	gu[gloc] = 0.0;
	gloc += 1;
  }

  // Bounds on voltage angles
  for(int i=0; i < Nvar_1st; i++) {
	/*  and bounds on real power mismatch equality constraints */
	xl[i] = -M_PI; xu[i] = M_PI;
  }

}


void
DCOPFLOW::VarAndConBounds_2nd_Partition(const int scen, double *xl,double *xu,double *gl,double *gu)
{
  DCLINE         *line;
  DCBUS          *bus;
  int       loc=0,gloc=0;
  
  // Line flow inequality constraints for cuts
  for(int i=0; i < ps->Nbranch; i++) {
  	line = &ps->line[i];	
	if(line->partID == scen) {
	  gl[gloc] = -(line->flowLimit/ps->MVAbase);
	  gu[gloc] =  (line->flowLimit/ps->MVAbase);
	  gloc += 1;	  
	}       
  }
  assert(gloc == Nconineq_2nd[scen]);

  for(int i=0; i < ps->Nbus; i++) {
  	bus = &ps->bus[i];
	if(bus->partID == scen){
	  /* Bounds on voltage angles and bounds on node balance equality constraints */
	  xl[loc] = -M_PI; xu[loc] = M_PI;
	  gl[gloc] = 0.0;  gu[gloc] = 0.0;
	  loc++;
	  gloc++;

	  /* Bounds on real power generation */
	  for(int k=0; k < bus->ngen; k++) {
		DCGEN *gen;
		bus->DCBUSGetGen(k,&gen);
		xl[loc] = gen->pb/ps->MVAbase; xu[loc] = gen->pt/ps->MVAbase;
		loc++;
	  }

	  // ref bus constraint
	  if(refBusID == bus->bus_i && !(ps->bus[refBusID].haveCutLine) ){ 
	    gl[gloc] = 0.0;  gu[gloc] = 0.0;
	    gloc += 1;	
	  }	
	}
  }

  for(int i=0; i < numDummyVar[scen]; i++) {
	/* Bounds on voltage angles for dummy variables */
	xl[loc] = -M_PI; xu[loc] = M_PI;
	loc++;
  }

  for(int i=0; i < numDummyCon[scen]; i++) {
	/* Bounds on dummy constraint */
	gl[gloc] = 0.0; gu[gloc] = 0.0; 
	gloc++;
  }

  assert(loc  == Nvar_2nd[scen]);	
  assert(gloc == Ncon_2nd[scen]);
}


void 
DCOPFLOW::ObjGradient_Lin_1st_Partition(double *obj_coef)
{
  DCBUS     *bus;
  DCGEN     *gen;
  int       loc;

  for(int i=0;i<Nvar_1st;i++) obj_coef[i]=0;
}


void 
DCOPFLOW::ObjGradient_Lin_2nd_Partition(const int scen, double *obj_coef)
{
  DCBUS     *bus;
  int       loc=0;

  int nVar_loc = Nvar_2nd[scen];
  int nCon_loc = Ncon_2nd[scen];
  
  for(int i=0;i<nVar_loc;i++) obj_coef[i]=0;

  for(int i=0; i < ps->Nbus; i++) {
  	bus = &ps->bus[i];
	if(bus->partID == scen){
	  loc++;
	  /* linear coefficient of real power generation */
	  for(int k=0; k < bus->ngen; k++) {
		DCGEN *gen;
		bus->DCBUSGetGen(k,&gen);
		obj_coef[loc] = ps->MVAbase*( gen->cost_beta );
		loc++;
	  }
	}
  }
}


void 
DCOPFLOW::ObjGradient_Quad_1st_Partition(int *irow, int *jcol, double *obj_quad)
{

}

int  
DCOPFLOW::GetHesNNZ_2nd_Partition(const int scen)
{
  DCBUS     *bus;
  DCGEN     *gen;
  int       nnz_part=0;

  for(int i=0; i < ps->Nbus; i++) {
    bus = &ps->bus[i];
	if(bus->partID == scen){
	  nnz_part += bus->ngen;
	}
  }
  return nnz_part;
}


void 
DCOPFLOW::ObjGradient_Quad_2nd_Partition(const int scen, int *irow, int *jcol, double *obj_quad)
{
  DCBUS     *bus;
  DCGEN     *gen;
  int       loc=0;
  int		ctr=0;

  for(int i=0; i < ps->Nbus; i++) {
    bus = &ps->bus[i];
	if(bus->partID == scen){
	  loc++;
	  /* quadratic coefficient of real power generation */
	  for(int k=0; k < bus->ngen; k++) {
		DCGEN *gen;
		bus->DCBUSGetGen(k,&gen);
	    irow[ctr] = loc; 
	    jcol[ctr] = loc; 		
        obj_quad[ctr] = 2.0 * gen->cost_alpha * ps->MVAbase * ps->MVAbase;
		loc++;
		ctr++;
	  }
	}
  }
  assert(loc == Nvar_2nd[scen]-numDummyVar[scen]);
}




int 
DCOPFLOW::GetJacNNZ_1st_Partition()
{
  int nnz = 0;
  
  DCBUS          *bus;
  DCLINE         *line;
  
  int       nconnlines;
  DCLINE   **connlines;

  if(nnz_jac_g_1st==-1){
    // Line flow inequality constraints for cuts
    for(int i=0; i < ps->Nbranch; i++) {
  	  line = &ps->line[i];
	  if(line->partID==-1) {
		nnz += 2;
	  }
    }
    // ref bus constraint
    if(ps->bus[refBusID].haveCutLine) {
  	  nnz += 1;
    }
	nnz_jac_g_1st = nnz;
  }else{
	nnz = nnz_jac_g_1st;
  }  
  
  return nnz;
}


int 
DCOPFLOW::GetJacNNZ_2nd_Partition(const int scen)
{
  int nnz = 0;
  
  DCBUS          *bus;
  DCLINE         *line;
  
  int       nconnlines;
  DCLINE   **connlines;

  if(nnz_jac_g_2nd[scen]==-1){
    // Local line flow inequality constraints 
    for(int i=0; i < ps->Nbranch; i++) {
  	  line = &ps->line[i];	
	  if(line->partID==scen) {
		nnz += 2;
	  }       
    }
	
    //Local KCL constraint:  \sum pgen  +  \sum flow_in -  \sum flow_out = d
    for(int i=0; i < ps->Nbus; i++) {
      bus = &ps->bus[i];
	  if(bus->partID==scen){
		nnz += bus->ngen;
		bus->DCBUSGetSupportingLines(&nconnlines,&connlines);
		nnz += 2*nconnlines;
	    
	    // ref bus constraint
	    if(refBusID == bus->bus_i && !(ps->bus[refBusID].haveCutLine) ){ 
		  nnz += 1;
	    }	
      }
    }	
	
	for(int i=0; i < numDummyCon[scen]; i++) {
	  /* Bounds on dummy constraint: x_2nd_dummy - x_1st=0 */
	  nnz += 1;
	}
	
	nnz_jac_g_2nd[scen] = nnz;
  }else{
	nnz = nnz_jac_g_2nd[scen];
  }  
  
  return nnz;
}


int 
DCOPFLOW::GetJacNNZ_Link_Partition(const int scen)
{
  int nnz = 0;
  
  DCBUS          *bus;
  DCLINE         *line;
  
  int       nconnlines;
  DCLINE   **connlines;

  if(nnz_jac_g_Link[scen]==-1){	
	nnz = numDummyCon[scen];

	nnz_jac_g_Link[scen] = nnz;
  }else{
	nnz = nnz_jac_g_Link[scen];
  }  
  
  return nnz;
}


void
DCOPFLOW::GetJac_1st_Partition(int *row, int *col, double *ele)
{  
  DCBUS          *bus;
  DCLINE         *line;
  
  int       nconnlines;
  DCLINE   **connlines;

  DCBUS	  *busf,*bust;
  DCBUS	  **connbuses;  

  int 	busfID, bustID;
  int   gloc=0;
  int 	ctr=0;

  // Line flow inequality constraints for cuts:  f_min <= flow_ij = -V^2/r * (phase_i - phase_j) <= f_max
  for(int i=0; i < ps->Nbranch; i++){
  	line = &ps->line[i];

	if(line->partID==-1) {
	  line->DCLINEGetConnectedBuses(&connbuses);		
	  busf = connbuses[0]; busfID = busf->bus_i;
	  bust = connbuses[1]; bustID = bust->bus_i;

      row[ctr] = gloc;						row[ctr+1] = gloc;
      col[ctr] = busMap_AllTo1st[busfID]; 	col[ctr+1] = busMap_AllTo1st[bustID];
	  ele[ctr] =  1.0/line->reactance;		ele[ctr+1] = -1.0/line->reactance;
      ctr  += 2;
      gloc += 1;
    }
  }
  
  // ref bus constraint:  phase = 0;
  if(ps->bus[refBusID].haveCutLine) {
    row[ctr] = gloc;
    col[ctr] = busMap_AllTo1st[refBusID];
	ele[ctr] = 1.0;
    ctr  += 1;
    gloc += 1;
  }

  assert(ctr == nnz_jac_g_1st);
  assert(gloc == Ncon_1st);
}




void 
DCOPFLOW::GetJac_2nd_Link_Partition(	const int scen, int *row, int *col, double *ele, 
  											int *row_link, int *col_link, double *ele_link)
{
  DCBUS          *bus;
  DCLINE         *line;
  
  int       nconnlines;
  DCLINE   **connlines;

  DCBUS	  *busf,*bust;
  DCBUS	  **connbuses;  

  int 	busfID, bustID;
  
  int 	xlocf, xloct;
  int   xloc=0, gloc=0;
  int 	ctr=0;
  int 	ctr_link=0;


  // Local line flow inequality constraints:  f_min <= flow_ij = -V^2/r * (phase_i - phase_j) <= f_max 
  for(int i=0; i < ps->Nbranch; i++) {
    line = &ps->line[i];	
	if(line->partID==scen) {
	  line->DCLINEGetConnectedBuses(&connbuses);		
	  busf = connbuses[0]; 	busfID = busf->bus_i;
	  bust = connbuses[1]; 	bustID = bust->bus_i;
	  busf->DCBUSGetVariableLocation_Part(&xlocf);
	  bust->DCBUSGetVariableLocation_Part(&xloct);

      row[ctr] = gloc; 						row[ctr+1] = gloc;
      col[ctr] = xlocf; 					col[ctr+1] = xloct;	  
	  ele[ctr] = 1.0/line->reactance;		ele[ctr+1] = -1.0/line->reactance;
      ctr  += 2;
      gloc += 1;
    }
  }
  assert(gloc == Nconineq_2nd[scen]);

  bool *usedBusID   	= new bool[ps->Nbus];  
  int  *dummyBusID 		= new int[ps->Nbus];  
  int  *dummyVarID 		= new int[ps->Nbus];  
  
  int  num_usedDummyID	= 0;
  for(int i=0; i < ps->Nbus; i++){ 
  	usedBusID[i]   		= false; 
	dummyBusID[i] 		= -1;
	dummyVarID[i] 		= -1;
  }
  
  for(int i=0; i < ps->Nbus; i++) {
    bus = &ps->bus[i];
	if(bus->partID==scen) {
	  //Local node balance constraint:  \sum pgen	+  \sum flow_in -  \sum flow_out = d
	  bus->DCBUSGetVariableLocation_Part(&xloc);
	  xloc++;
	  for(int k=0; k < bus->ngen; k++) {	  
		row[ctr] = gloc;	 
		col[ctr] = xloc;
		ele[ctr] = 1.0;
		ctr  += 1;
		xloc += 1;
	  }
 
	  bus->DCBUSGetSupportingLines(&nconnlines,&connlines);
	  for(int k=0; k < nconnlines;k++) {
		line = connlines[k];
		line->DCLINEGetConnectedBuses(&connbuses);		
		busf = connbuses[0];  busfID = busf->bus_i;
		bust = connbuses[1];  bustID = bust->bus_i;	
		row[ctr] = gloc;		row[ctr+1] = gloc;

		if(bus == busf) {	
		  ele[ctr]	 = -1.0/line->reactance;
		  ele[ctr+1] =  1.0/line->reactance;
		}else{
		  ele[ctr]	 =  1.0/line->reactance;
		  ele[ctr+1] = -1.0/line->reactance;
		}

		if(line->partID==scen){
		  // local line
		  busf->DCBUSGetVariableLocation_Part(&xlocf);
		  bust->DCBUSGetVariableLocation_Part(&xloct);		  
		}else{
		  assert(line->partID==-1);
		  // cut line, one node is in this partition
		  if(busf->partID == scen){
		  	// from bus is in this part
			busf->DCBUSGetVariableLocation_Part(&xlocf);
			if(!usedBusID[bustID]){
			  xloct = Nvar_2nd[scen] - numDummyVar[scen] + num_usedDummyID;
			  num_usedDummyID++;
			  usedBusID[bustID]  = true;
			  dummyVarID[bustID] = xloct;
			}else{
			  xloct = dummyVarID[bustID];
			}		
		  }else{
			// to bus is in this part
  			bust->DCBUSGetVariableLocation_Part(&xloct);
			if(!usedBusID[busfID]){
			  xlocf = Nvar_2nd[scen] - numDummyVar[scen] + num_usedDummyID;
			  num_usedDummyID++;
			  usedBusID[busfID]  = true;
			  dummyVarID[busfID] = xlocf;
			}else{
			  xlocf = dummyVarID[busfID];
			}			
		  }
		}
		col[ctr] = xlocf; col[ctr+1] = xloct;		
		ctr += 2;
	  }
	  gloc += 1;

	  // ref bus constraint:  phase = 0;
	  if(refBusID == bus->bus_i && !(ps->bus[refBusID].haveCutLine) ){ 
		bus->DCBUSGetVariableLocation_Part(&xloc);
		row[ctr] = gloc;
		col[ctr] = xloc;		
		ele[ctr] = 1.0;
		ctr  += 1;
	    gloc += 1;	
	  }	
  	}
  }	

  assert(gloc == Ncon_2nd[scen]-numDummyCon[scen]);

  // for dummy constraint on bus linked with cuts:  phase angle = global phase angel
  for(int i=0; i < ps->Nbus; i++) {
    bus = &ps->bus[i];
	if(bus->partID==scen && bus->haveCutLine) {
	  bus->DCBUSGetVariableLocation_Part(&xloc);
	  row[ctr] = gloc;
	  col[ctr] = xloc;
	  ele[ctr] = 1.0;

	  row_link[ctr_link] = gloc;
	  col_link[ctr_link] = busMap_AllTo1st[bus->bus_i];
	  ele_link[ctr_link] = -1.0;
	  
	  ctr  		+= 1;
	  ctr_link  += 1;	  
	  gloc 		+= 1;	
	}
  }
  
  // for dummy constraint on dummy bus: phase angle = global phase angel
  for(int i=0; i < ps->Nbus; i++){ 
  	if(usedBusID[i]){
	  row[ctr] = gloc;
	  col[ctr] = dummyVarID[i];	
	  ele[ctr] = 1.0;

	  row_link[ctr_link] = gloc;
	  col_link[ctr_link] = busMap_AllTo1st[bus->bus_i];
	  ele_link[ctr_link] = -1.0;
	  
	  ctr  		+= 1;
	  ctr_link  += 1;	  
	  gloc 		+= 1;
	}
  }

  
  assert(ctr  		== nnz_jac_g_2nd[scen]);
  assert(ctr_link  	== nnz_jac_g_Link[scen]);
  assert(gloc == Ncon_2nd[scen]);

  delete []usedBusID;
  delete []dummyBusID;
  delete []dummyVarID;


}

