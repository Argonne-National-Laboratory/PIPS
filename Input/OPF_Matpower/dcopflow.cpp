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
	  setupcalled(false), setupcalled_part(false),setupcalled_aggregation(false),
	  lambda_g(NULL), lambda_xl(NULL), lambda_xu(0)
{
  if(ps->Nparts==1)
  	DCOPFLOWSetUp();
  else{
  	DCOPFLOWSetUp_Partition();
	SetUp_Aggregation();
  }
}

DCOPFLOW::~DCOPFLOW()
{
  if(setupcalled){
	delete[] (Xl);
	delete[] (Xu);
	delete[] (Gl);
	delete[] (Gu);
	delete[] (gradobj);
	delete[] (lambda_g);
	delete[] (lambda_xl);
	delete[] (lambda_xu);	
  }	
  if(setupcalled_part){
    delete[](busMap_AllTo1st);
    delete[](Nvar_2nd);  
    delete[](Ncon_2nd); 
    delete[](Nconeq_2nd); 
    delete[](Nconineq_2nd); 
	delete[](numDummyVar);
	delete[](numDummyCon);
	delete[] nnz_jac_g_2nd;
	delete[] nnz_jac_g_Link;
    for (int scen=0;scen<Nparts;scen++){
	  delete [] dummyBusVarID[scen];
    }
    delete [] dummyBusVarID;	
  }

  if(setupcalled_aggregation){
    for (int scen=0;scen<Nparts;scen++){
	  delete [] locVarMap_Agg[scen];
	  delete [] locConMap_Agg[scen];
    }
    delete [] locVarMap_Agg;
    delete [] locConMap_Agg;
	
	delete [] firstVarMap_Agg;
	delete [] firstConMap_Agg;	
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
    Gl[gloc] = -(line->flowLimit);
    Gu[gloc] =  (line->flowLimit);
       
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
      Xl[loc] = gen->pb; Xu[loc] = gen->pt;
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
  int xloc=0,xloct,xlocf;
  
  Nparts = ps->Nparts;

  Nvar_2nd 			= new int[Nparts];
  Ncon_2nd 			= new int[Nparts];  
  Nconeq_2nd 		= new int[Nparts];
  Nconineq_2nd 		= new int[Nparts];   

  numDummyVar		= new int[Nparts];
  numDummyCon		= new int[Nparts];
  nnz_jac_g_2nd 	= new int[Nparts];
  nnz_jac_g_Link 	= new int[Nparts];  

  dummyBusVarID 	= new int*[Nparts];  
  

  busMap_AllTo1st 	= new int[ps->Nbus]; 

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
	dummyBusVarID[i] = new int[ps->Nbus]; 
	for(int k=0; k < ps->Nbus; k++) {
	 dummyBusVarID[i][k] = -1;
    } 
  }  

  for(int i=0; i < ps->Nbus; i++) {
	busMap_AllTo1st[i]=-1;
  }  

  // define num of constraints on the cuts, 
  Nvar_1st  	= 0;
  Nslack_1st	= 0;
  Nconeq_1st 	= 0;
  Nconineq_1st 	= 0;

  Nslack_1st	+= ps->Ncuts;	// add slack for inequalities
  if(ps->bus[refBusID].haveCutLine) Nconeq_1st++;    
  Nconeq_1st    += ps->Ncuts;	// add slack for inequalities
  Ncon_1st = Nconineq_1st + Nconeq_1st;

  int  *num_usedDummyID	= new int[ps->Nparts];
  int  *find_cuts		= new int[ps->Nparts];
  bool **usedBusID;  
  usedBusID = new bool*[Nparts];
  for(int i=0; i<Nparts; i++) {usedBusID[i] = new bool[ps->Nbus]; num_usedDummyID[i]=0;find_cuts[i]=0;}
  for(int i=0; i<Nparts; i++) for(int j=0; j<ps->Nbus; j++) usedBusID[i][j] = false;

  // num of var in each partition 
  for(int i=0; i < ps->Nbus; i++) {
    bus = &ps->bus[i];
	bus->startloc_Part = Nvar_2nd[bus->partID];

	// local phase angle variable
	Nvar_2nd[bus->partID]++;
	// local power generation variable 
	Nvar_2nd[bus->partID] += bus->ngen; 

	// local  node balance constraint. (= #node in each partition)
	Nconeq_2nd[bus->partID]++;
	// ref bus constraint
	if(refBusID == bus->bus_i && !(ps->bus[refBusID].haveCutLine) ) Nconeq_2nd[bus->partID]++;
  }
  
  // find num of dummy var/con for cuts
  for(int i=0; i < ps->Nbus; i++) {
	bus = &ps->bus[i];
	int scen = bus->partID;
	bus->DCBUSGetSupportingLines(&nconnlines,&connlines);

	// find cut lines linked to this bus
	if(bus->haveCutLine){
	  numDummyCon[bus->partID]++;
	  if(busMap_AllTo1st[i]==-1){
		busMap_AllTo1st[i] = Nvar_1st;	  
		Nvar_1st++;
	  }		
	  for(int k=0; k < nconnlines && find_cuts[scen]< ps->num_CutInEachPart[scen];k++) {
		line = connlines[k];
		if(line->partID == -1)
		{
		  find_cuts[scen]++;
		  line->DCLINEGetConnectedBuses(&connbuses);		
		  busf = connbuses[0];  busfID = busf->bus_i;
		  bust = connbuses[1];  bustID = bust->bus_i;					
		  // cut line, one node is in this partition.  add phase angle for dummy bus
		  if(busf->partID == bus->partID){
			// from bus is in this part
			if(!usedBusID[scen][bustID]){
			  numDummyVar[scen]++;
			  numDummyCon[scen]++;
			  dummyBusVarID[scen][bustID] = Nvar_2nd[scen] + num_usedDummyID[scen];
			  num_usedDummyID[scen]++;			  
			  usedBusID[scen][bustID] = true;
			} 	  			
		  }else{
			// to bus is in this part
			if(!usedBusID[scen][busfID]){
			  numDummyVar[scen]++;
			  numDummyCon[scen]++;
			  dummyBusVarID[scen][busfID] = Nvar_2nd[scen] + num_usedDummyID[scen];
			  num_usedDummyID[scen]++;
			  usedBusID[scen][busfID] = true;
			} 		  
		  }
		}
      }
	}

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
  Nvar		= Nvar_1st + Nslack_1st;    // add slack var
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

  delete [] num_usedDummyID; 
  delete [] find_cuts; 

}


void
DCOPFLOW::VarAndConBounds_1st_Partition(double *xl,double *xu,double *gl,double *gu)
{
  DCLINE         *line;
  DCBUS          *bus;
  int       loc,gloc=0;

  int findSlack = 0;

  // Bounds on voltage angles
  for(int i=0; i < Nvar_1st; i++) {
	/*  and bounds on real power mismatch equality constraints */
	xl[i] = -M_PI; xu[i] = M_PI;
  }
  
  // Line flow inequality constraints for cuts:  f_min <= flow_ij = -V^2/r * (phase_i - phase_j) <= f_max
  // add slack as :   -V^2/r * (phase_i - phase_j)  - slack = 0,
  //   f_min <= slack <= f_max
  for(int i=0; i < ps->Nbranch; i++) {
  	line = &ps->line[i];	
	if(line->partID==-1) {
	  gl[gloc] = 0.0;
	  gu[gloc] = 0.0; 	  
	  xl[Nvar_1st+findSlack] = -(line->flowLimit); 
	  xu[Nvar_1st+findSlack] =  (line->flowLimit);
	  findSlack += 1;
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
	  gl[gloc] = -(line->flowLimit);
	  gu[gloc] =  (line->flowLimit);
	  gloc += 1;	  
	}       
  }
  assert(gloc == Nconineq_2nd[scen]);

  for(int i=0; i < ps->Nbus; i++) {
  	bus = &ps->bus[i];
	if(bus->partID == scen){
	  /* Bounds on voltage angles and bounds on node balance equality constraints */
	  xl[loc] = -M_PI; xu[loc] = M_PI;
	  gl[gloc] = bus->loads;  
	  gu[gloc] = bus->loads;
	  loc++;
	  gloc++;

	  /* Bounds on real power generation */
	  for(int k=0; k < bus->ngen; k++) {
		DCGEN *gen;
		bus->DCBUSGetGen(k,&gen);
		xl[loc] = gen->pb; xu[loc] = gen->pt;
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
    // Line flow inequality constraints for cuts (adding slack to make it equality constraint)
    for(int i=0; i < ps->Nbranch; i++) {
  	  line = &ps->line[i];
	  if(line->partID==-1) {
		nnz += 2 + 1;  // 2 for th phase angle, 1 for the slack
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

  int findSlack = 0;

  // Line flow inequality constraints for cuts:  f_min <= flow_ij = -V^2/r * (phase_i - phase_j) <= f_max
  // add slack as :   -V^2/r * (phase_i - phase_j)  - slack = 0,
  //   f_min <= slack <= f_max
  for(int i=0; i < ps->Nbranch; i++){
  	line = &ps->line[i];

	if(line->partID==-1) {
	  line->DCLINEGetConnectedBuses(&connbuses);		
	  busf = connbuses[0]; busfID = busf->bus_i;
	  bust = connbuses[1]; bustID = bust->bus_i;

      row[ctr] = gloc;						row[ctr+1] = gloc;						row[ctr+2] = gloc;
      col[ctr] = busMap_AllTo1st[busfID]; 	col[ctr+1] = busMap_AllTo1st[bustID];	col[ctr+2] = Nvar_1st + findSlack;
	  ele[ctr] =  1.0/line->reactance;		ele[ctr+1] = -1.0/line->reactance;		ele[ctr+2] = -1.0;
      ctr  += 3;
      gloc += 1;
	  findSlack += 1;
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
  
  int  num_usedDummyID	= 0;
  for(int i=0; i < ps->Nbus; i++){ 
  	usedBusID[i]   		= false; 
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
//			if(!usedBusID[bustID]){
//			  xloct = Nvar_2nd[scen] - numDummyVar[scen] + num_usedDummyID;
//			  num_usedDummyID++;
//			  usedBusID[bustID]  = true;
//			  dummyBusVarID[scen][bustID] = xloct;
//			}else
			{
			  usedBusID[bustID]  = true;
			  xloct = dummyBusVarID[scen][bustID];
			}		
		  }else{
			// to bus is in this part
  			bust->DCBUSGetVariableLocation_Part(&xloct);
//			if(!usedBusID[busfID]){
//			  xlocf = Nvar_2nd[scen] - numDummyVar[scen] + num_usedDummyID;
//			  num_usedDummyID++;
//			  usedBusID[busfID]  = true;
//			  dummyBusVarID[scen][busfID] = xlocf;
//			}else
			{
			  usedBusID[busfID]  = true;
			  xlocf = dummyBusVarID[scen][busfID];
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
  	bus = &ps->bus[i];
  	if(usedBusID[i]){
	  row[ctr] = gloc;
	  col[ctr] = dummyBusVarID[scen][i];	
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


}

void 
DCOPFLOW::SetUp_Aggregation()
{
  if(!setupcalled_part) DCOPFLOWSetUp_Partition();
  
  DCBUS          *bus;
  DCLINE		 *line;
 
  int   xloc=0;
  nvar_aggregation = 0;
  ncon_aggregation = 0;

  // for second stage
  for(int i=0; i < ps->Nbus; i++) {
    bus = &ps->bus[i];
	if(bus->haveCutLine) nvar_aggregation++;  	//  var: phase angle for the bus with cut
	nvar_aggregation += bus->ngen;				//  var: power generator on the bus
  }
  for (int scen=0;scen<Nparts;scen++){
	nvar_aggregation += numDummyVar[scen];		//  var: phase angle from the dummy bus
	ncon_aggregation += numDummyCon[scen];		//  con: dummy constraint to set phase angle as the global var
	ncon_aggregation += 1;						//  con: power gen and demand for this part
  }

  // for 1st stage
  ncon_aggregation	+= Nslack_1st;				//  con: slack line flow constraint on the cuts
  nvar_aggregation	+= Nvar_1st + Nslack_1st;	//  var: global phase angle variable and slack line flow variable on the cuts


  ncon_aggregation	+= 1;	// reference bus constraint


  firstVarMap_Agg 	= new int[Nvar_1st+Nslack_1st];
  firstConMap_Agg   = new int[Ncon_1st];
  	
  locVarMap_Agg = new int*[Nparts];
  locConMap_Agg = new int*[Nparts];

  for (int scen=0;scen<Nparts;scen++){
	locVarMap_Agg[scen] = new int[Nvar_2nd[scen]];
	locConMap_Agg[scen] = new int[Ncon_2nd[scen]];
	for(int i=0;i<Nvar_2nd[scen];i++){ locVarMap_Agg[scen][i] = -1; locConMap_Agg[scen][i] = -1;}
  }

  // define local var to aggregation var map
  int find_nvar=0;
  for (int scen=0;scen<Nparts;scen++){
    for(int i=0; i < ps->Nbus; i++) {
      bus = &ps->bus[i];	
	  if(bus->partID==scen) {
	    bus->DCBUSGetVariableLocation_Part(&xloc);
	    // phase angle var for node connected with cuts
	    if(bus->haveCutLine){
	  	  locVarMap_Agg[scen][xloc] = find_nvar++;
	    }
	    //power generation at this bus
	    xloc++;
	    for(int k=0; k < bus->ngen; k++) {	  
		  locVarMap_Agg[scen][xloc] = find_nvar++;
		  xloc++;
	    }
	  }
    }
    // local dummy vars
    for(int i=0; i < numDummyVar[scen]; i++) {
	  xloc = Nvar_2nd[scen] - numDummyVar[scen] + i;
	  locVarMap_Agg[scen][xloc] = find_nvar++;
    }
  }
  assert(find_nvar == nvar_aggregation - Nvar_1st - Nslack_1st);

  // attach 1st var in the end
  for(int i=0; i < Nvar_1st+Nslack_1st; i++) {
	firstVarMap_Agg[i] = find_nvar++;
  }  


  // define local con to aggregation var map
  int find_ncon=0,gloc=0;
  bool findRefbus=false;
  for (int scen=0;scen<Nparts;scen++){
  	gloc = 0;
	// SKIP the Line flow inequality constraints for cuts
	for(int i=0; i < ps->Nbranch; i++) {
	  line = &ps->line[i];	  
	  if(line->partID == scen) {
		gloc++ ;		
	  } 	  
	}
	assert(gloc == Nconineq_2nd[scen]);

	for(int i=0; i < ps->num_BusInEachPart[scen]; i++) {
	  bus = &ps->bus[ ps->BusInEachPart[scen][i] ]; 
	  assert(bus->partID==scen); 
	  /* Bounds on voltage angles and bounds on node balance equality constraints */
	  locConMap_Agg[scen][gloc] = find_ncon;
	  gloc++;
	
	  // ref bus constraint
	  if(refBusID == bus->bus_i && !(ps->bus[refBusID].haveCutLine) ){ 
		locConMap_Agg[scen][gloc] = find_ncon+1;
		findRefbus = true;
		gloc++; 
	  } 
	}
	find_ncon++;
	if(findRefbus) find_ncon++;

	/* dummy constraint */
	for(int i=0; i < numDummyCon[scen]; i++) {
	  locConMap_Agg[scen][gloc] = find_ncon;
	  gloc++;	find_ncon++;	  
	}
  }
  assert(find_ncon == ncon_aggregation - Ncon_1st);

  // attach 1st con in the end
  for(int i=0; i < Ncon_1st; i++) {
	firstConMap_Agg[i] = find_ncon++;
  }  





  nnz_jac_Aggregation=-1;

  setupcalled_aggregation = true;
}



// only support equality cons
void
DCOPFLOW::VarAndConBounds_Aggregation(double *xl,double *xu,double *gl,double *gu)
{
  if(!setupcalled_aggregation) SetUp_Aggregation();

  DCLINE         *line;
  DCBUS          *bus;
  int       loc,gloc=0;

  int findSlack = 0;

  bool findRefbus=false;

  // define local var lb and ub
  int find_nvar=0;
  for (int scen=0;scen<Nparts;scen++){
  	//Local generation-demand balance constraint:  \sum pgen	+  \sum flow_in -  \sum flow_out = d
	gl[gloc] = 0.0;  gu[gloc] = 0.0;
  	
    for(int i=0; i < ps->num_BusInEachPart[scen]; i++) {
      bus = &ps->bus[ ps->BusInEachPart[scen][i] ];	
	  assert(bus->partID==scen); 
	  gl[gloc] += bus->loads;  gu[gloc] += bus->loads;
	  
	  // phase angle var for node connected with cuts
	  if(bus->haveCutLine){
		/*  and bounds on phase angle */
//		xl[find_nvar] = -M_PI; xu[find_nvar] = M_PI;	
		find_nvar++;
	  }
	  //power generation at this bus
	  for(int k=0; k < bus->ngen; k++) {	
		DCGEN *gen;
		bus->DCBUSGetGen(k,&gen);			
		xl[find_nvar] = gen->pb; xu[find_nvar] = gen->pt;
		find_nvar++;
	  }
	  if(refBusID == bus->bus_i && !(ps->bus[refBusID].haveCutLine) ) findRefbus = true;
    }
	gloc += 1;
	
    // ref bus constraint in this part: phase = 0;
    if(findRefbus){ 
	  gl[gloc] = 0.0;  gu[gloc] = 0.0;
	  gloc += 1;	
    } 	
	
    // local dummy phase angle vars
    for(int i=0; i < numDummyVar[scen]; i++) {
	  xl[find_nvar+i] = -M_PI; xu[find_nvar+i] = M_PI;	
    }
	find_nvar+=numDummyVar[scen];
  }
  assert(find_nvar == nvar_aggregation - Nvar_1st - Nslack_1st);


  // for 1st stage var/con

  // Bounds on voltage angles
  for(int i=find_nvar; i < find_nvar+Nvar_1st; i++) {
	/*  and bounds on real power mismatch equality constraints */
	xl[i] = -M_PI; xu[i] = M_PI;
  }
  find_nvar += Nvar_1st;

  // Line flow inequality constraints for cuts:  f_min <= flow_ij = -V^2/r * (phase_i - phase_j) <= f_max
  // add slack as :   -V^2/r * (phase_i - phase_j)  - slack = 0,
  //   f_min <= slack <= f_max
  for(int i=0; i < ps->Nbranch; i++) {
  	line = &ps->line[i];	
	if(line->partID==-1) {
	  gl[gloc] = 0.0;
	  gu[gloc] = 0.0; 	  
	  xl[find_nvar+findSlack] = -(line->flowLimit); 
	  xu[find_nvar+findSlack] =  (line->flowLimit);
	  findSlack += 1;
	  gloc += 1;
	}       
  }
  find_nvar += findSlack;

  // ref bus constraint
  if(ps->bus[refBusID].haveCutLine) {
	gl[gloc] = 0.0;  
	gu[gloc] = 0.0;
	gloc += 1;
  }  
  assert(gloc == ncon_aggregation);
  assert(find_nvar == nvar_aggregation);
}


// only for equalities constraint
void
DCOPFLOW::EqConBounds_Aggregation(double *b)
{
  if(!setupcalled_aggregation) SetUp_Aggregation();

  DCLINE         *line;
  DCBUS          *bus;
  int       loc,gloc=0;

  int findSlack = 0;

  bool findRefbus=false;

  // define local constraint bound
  for (int scen=0;scen<Nparts;scen++){
  	//Local generation-demand balance constraint:  \sum pgen	+  \sum flow_in -  \sum flow_out = d
	b[gloc] = 0.0;
  	
    for(int i=0; i < ps->num_BusInEachPart[scen]; i++) {
      bus = &ps->bus[ ps->BusInEachPart[scen][i] ];	
	  assert(bus->partID==scen); 
	  b[gloc] += bus->loads;
	  
	  if(refBusID == bus->bus_i && !(ps->bus[refBusID].haveCutLine) ) findRefbus = true;
    }
	gloc += 1;
	
    // ref bus constraint in this part: phase = 0;
    if(findRefbus){ 
	  b[gloc] = 0.0;
	  gloc += 1;	
    } 	

    for(int i=0; i < numDummyCon[scen]; i++) {
	  /* Bounds on dummy constraint */
	  b[gloc] = 0.0; 
	  gloc++;
    }  
  }

  // Line flow inequality constraints for cuts:  f_min <= flow_ij = -V^2/r * (phase_i - phase_j) <= f_max
  // add slack as :   -V^2/r * (phase_i - phase_j)  - slack = 0,
  //   f_min <= slack <= f_max
  for(int i=0; i < ps->Nbranch; i++) {
  	line = &ps->line[i];	
	if(line->partID==-1) {
	  b[gloc] = 0.0;
	  gloc += 1;
	}       
  }


  // ref bus constraint
  if(ps->bus[refBusID].haveCutLine) {
	b[gloc] = 0.0;  
	gloc += 1;
  }  
  assert(gloc == ncon_aggregation);
}


// only for equalities constraint
void
DCOPFLOW::objLinGrad_Aggregation(double *obj_coef)
{
  if(!setupcalled_aggregation) SetUp_Aggregation();
  DCBUS	  *bus;
  DCGEN	  *gen;
  int 	  loc=0;

  
  int firstStVarStartID = nvar_aggregation - Nvar_1st - Nslack_1st;

  for(int i=firstStVarStartID;i<nvar_aggregation;i++) obj_coef[i]=0.;
  
  for(int scen=0;scen<Nparts;scen++){
	int nVar_loc = Nvar_2nd[scen];
    for(int i=0;i<nVar_loc;i++) 
	  obj_coef[locVarMap_Agg[scen][i]]=0;
  }	
  
  for(int i=0; i < ps->Nbus; i++) {
	bus = &ps->bus[i];
	int scen = bus->partID;
	bus->DCBUSGetVariableLocation_Part(&loc);
	{
	  loc++;
	  /* linear coefficient of real power generation */
	  for(int k=0; k < bus->ngen; k++) {
	    DCGEN *gen;
	    bus->DCBUSGetGen(k,&gen);
	    obj_coef[locVarMap_Agg[scen][loc]] = ps->MVAbase*( gen->cost_beta );
	    loc++;
	  }
	}
  }
}







int 
DCOPFLOW::GetAggregationJacNNZ()
{
  if(!setupcalled_aggregation) SetUp_Aggregation();

  int nnz = 0;
  
  DCBUS          *bus;
  DCLINE         *line;
  
  int       nconnlines;
  DCLINE   **connlines;

  if(nnz_jac_Aggregation==-1){

  	// Jac: total number of generators, defined by node balance
  	for(int i=0; i < ps->Nbus; i++) {
	  bus = &ps->bus[i];
	  nnz += bus->ngen;
  	}   

    for (int scen=0;scen<Nparts;scen++){  
	  // Jac: phase angle in dummy constraints, phase angle = global var  
	  nnz += 2*numDummyCon[scen];
    }
	// Jac: each cut introduce 4 nnz:  dummt phase angle and local phase angle, for two side
	nnz += 4*ps->Ncuts;

	// ref bus constraint
	if(!ps->bus[refBusID].haveCutLine) {
	  nnz += 1;
	}


	// Jac from 1st stage
	nnz += nnz_jac_g_1st;
		
	nnz_jac_Aggregation = nnz;
  }else{
	nnz = nnz_jac_Aggregation;
  }
  return nnz;
}


// find sparsity pattern

void 
DCOPFLOW::GetPrecondMatrixJac_Aggregation(int *row, int *col, double *ele)
{
  if(!setupcalled_aggregation) SetUp_Aggregation();
#if 1
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
  int find_cuts=0;

  bool findRefbus = false;

  int firstStVarStartID = nvar_aggregation - Nvar_1st - Nslack_1st;
 
  for (int scen=0;scen<Nparts;scen++){
  	find_cuts=0;
  	//Local generation-demand balance constraint:  \sum pgen	+  \sum flow_in -  \sum flow_out = d
    for(int i=0; i < ps->num_BusInEachPart[scen]; i++) {
      bus = &ps->bus[ ps->BusInEachPart[scen][i] ];	
	  assert(bus->partID==scen); 
	  bus->DCBUSGetVariableLocation_Part(&xloc);
	  xloc++;
	  //power generation at this bus
	  for(int k=0; k < bus->ngen; k++) {	  
		row[ctr] = gloc;	 
		col[ctr] = locVarMap_Agg[scen][xloc] ;
		ele[ctr] = 1.0;
		ctr  += 1;
		xloc += 1;
	  }
	  // flow in cuts	  
	  if(bus->haveCutLine){
	  	bus->DCBUSGetSupportingLines(&nconnlines,&connlines);
		for(int k=0; k < nconnlines && find_cuts< ps->num_CutInEachPart[scen]; k++) {
		  line = connlines[k];
		  if(line->partID==-1){
			find_cuts++;
			line->DCLINEGetConnectedBuses(&connbuses);		
		 	busf = connbuses[0];  	busfID = busf->bus_i;
			bust = connbuses[1];  	bustID = bust->bus_i;	
			if(bus == busf) {	
			  // from bus is in this part
		      ele[ctr]	 = -1.0/line->reactance;	ele[ctr+1] =  1.0/line->reactance;
			  // phase angle on FROM bus
			  busf->DCBUSGetVariableLocation_Part(&xlocf);
			  // dummy phase angle on TO bus
			  xloct = dummyBusVarID[scen][bustID];
			}else{ 
			  // to bus is in this part
		  	  ele[ctr]	 =  1.0/line->reactance;	ele[ctr+1] = -1.0/line->reactance;
			  // phase angle on TO bus
			  bust->DCBUSGetVariableLocation_Part(&xloct);
			  // dummy phase angle on FROM bus
			  xlocf = dummyBusVarID[scen][busfID];			  
			}			
			row[ctr] = gloc;							row[ctr+1] = gloc;
			col[ctr] = locVarMap_Agg[scen][xlocf]; 		col[ctr+1] = locVarMap_Agg[scen][xloct];		
			ctr += 2;
		  }
		}
	  }
	  if(refBusID == bus->bus_i && !(ps->bus[refBusID].haveCutLine) ) findRefbus = true;
    }
	gloc += 1;
	
	// ref bus constraint in this part:	phase = 0;
	if(findRefbus){ 
	  bus = &ps->bus[ refBusID ];
	  bus->DCBUSGetVariableLocation_Part(&xloc);
	  row[ctr] = gloc;
	  col[ctr] = locVarMap_Agg[scen][xloc];		  
	  ele[ctr] = 1.0;
	  ctr  += 1;
	  gloc += 1;  
	} 
	
	// for dummy constraint on bus linked with cuts:  phase angle = global phase angel
	for(int i=0; i < ps->num_BusInEachPart[scen]; i++) {
	  bus = &ps->bus[ ps->BusInEachPart[scen][i] ];
	  assert(bus->partID==scen); 
	  if(bus->haveCutLine) {
		bus->DCBUSGetVariableLocation_Part(&xloc);
		row[ctr] = gloc;						row[ctr+1] = gloc;
		col[ctr] = locVarMap_Agg[scen][xloc];	col[ctr+1] = firstStVarStartID+busMap_AllTo1st[bus->bus_i];
		ele[ctr] = 1.0;							ele[ctr+1] = -1.0;
		ctr 	  += 2;
		gloc	  += 1;   
	  }
	}
	
	// for dummy constraint on dummy bus: phase angle = global phase angel
	for(int i=0; i < ps->Nbus; i++){ 
	  int locID_ori = dummyBusVarID[scen][i];
	  if(locID_ori>=0){
		row[ctr] = gloc;							row[ctr+1] = gloc;
		col[ctr] = locVarMap_Agg[scen][locID_ori];	col[ctr+1] = firstStVarStartID+busMap_AllTo1st[i];
		ele[ctr] = 1.0;								ele[ctr+1] = -1.0;
		ctr 	  += 2;
		gloc	  += 1;
	  }
	}
  }


  // Line flow inequality constraints for cuts:  f_min <= flow_ij = -V^2/r * (phase_i - phase_j) <= f_max
  // add slack as :   -V^2/r * (phase_i - phase_j)  - slack = 0,
  //   f_min <= slack <= f_max
  int findSlack=0;
  for(int i=0; i < ps->Nbranch; i++){
  	line = &ps->line[i];

	if(line->partID==-1) {
	  line->DCLINEGetConnectedBuses(&connbuses);		
	  busf = connbuses[0]; busfID = busf->bus_i;
	  bust = connbuses[1]; bustID = bust->bus_i;

      row[ctr] = gloc;						row[ctr+1] = gloc;						row[ctr+2] = gloc;
      col[ctr] 	 = firstStVarStartID + busMap_AllTo1st[busfID]; 	
	  col[ctr+1] = firstStVarStartID + busMap_AllTo1st[bustID];	
	  col[ctr+2] = firstStVarStartID + Nvar_1st + findSlack;
	  ele[ctr] =  1.0/line->reactance;		ele[ctr+1] = -1.0/line->reactance;		ele[ctr+2] = -1.0;
      ctr  += 3;
      gloc += 1;
	  findSlack += 1;
    }
  }

  // ref bus constraint:  phase = 0;
  if(ps->bus[refBusID].haveCutLine) {
    row[ctr] = gloc;
    col[ctr] = firstStVarStartID + busMap_AllTo1st[refBusID];
	ele[ctr] = 1.0;
    ctr  += 1;
    gloc += 1;
  }

#endif


}



int 
DCOPFLOW::GetAggregationHesNNZ()
{
  return ps->Ngen;
}


// find sparsity pattern

void 
DCOPFLOW::GetPrecondMatrixHes_Aggregation(int *row, int *col, double *ele)
{

  if(!setupcalled_aggregation) SetUp_Aggregation();
#if 1
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

  for (int scen=0;scen<Nparts;scen++){
  	//Local generation-demand balance constraint:  \sum pgen	+  \sum flow_in -  \sum flow_out = d
    for(int i=0; i < ps->num_BusInEachPart[scen]; i++) {
      bus = &ps->bus[ ps->BusInEachPart[scen][i] ];	
	  assert(bus->partID==scen); 
	  bus->DCBUSGetVariableLocation_Part(&xloc);
	  xloc++;
	  //power generation at this bus
	  for(int k=0; k < bus->ngen; k++) {	
		DCGEN *gen;
		bus->DCBUSGetGen(k,&gen);	  	
		row[ctr] = locVarMap_Agg[scen][xloc] ; 
		col[ctr] = locVarMap_Agg[scen][xloc] ;
		ele[ctr] = 2.0 * gen->cost_alpha * ps->MVAbase * ps->MVAbase;
		ctr  += 1;
		xloc += 1;
	  }
    }
  }
  
  assert(ctr == ps->Ngen);
#endif

}




#if 0
void 
DCOPFLOW::GetAggregationGrad()
{
  return ps->Ngen;
}

void 
DCOPFLOW::GetAggregationConRHS()
{
  return ps->Ngen;
}
#endif
