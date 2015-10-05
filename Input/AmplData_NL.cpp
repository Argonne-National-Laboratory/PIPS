/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "AmplData_NL.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>

#include "asl.h"
#include "asl_pfgh.h"
#include "getstub.h"

using namespace std;

/* static variables to be available to all routines in this file */
static fint nslacks = 0;
static int setLinear = 0;

int nobj = 0;
int hes_obj = 1;
int hes_con = 1;
int hes_tri = 1;
double OW[1] = {1};

static Long lprintLevel = 0;
const int nk = 1;
static keyword keywds[] = {      /* must be in alphabetical order */
  KW((char *) "print_level",    L_val, &lprintLevel,  
	 (char *) "Amount of output" )
};

// ASL options
Option_Info Oinfo = { (char*)"PIPS-NLP", (char*) "PIPS-NLP", (char*)"PIPS-NLP_options", keywds, nk };

void AmplData_NL::initialize(int nvar, int ncons) 
{
  collb.resize(nvar);
  colub.resize(nvar);
  objGrad.resize(nvar);
  rowlb.resize(ncons);
  rowub.resize(ncons);
  rownames.resize(ncons);
  colnames.resize(nvar);
  VarsStatus.resize(nvar);
  ConsStatus.resize(ncons);	
  didLoad = true;
}


//the number of nonzeros in the sparse Hessian W of the Lagrangian (if uptri = 0) or its upper triangle (if uptri = 1)
int AmplData_NL::Ampl_nnz_Hessian_Tri()
{
  ASL_pfgh *asl = locASL;
  int nzH = sphsetup(-1, hes_obj, hes_con, hes_tri); 
  return nzH;
}

ASL_pfgh* AmplData_NL::initASL(char *nlFileName[], AmplSuffix* amplSuffix) 
{
  // init asl structure
  ASL_pfgh *asl;
  FILE *nlFile;
  char *stub;
  int i;
	
  // set up option info
  Oinfo.options = keywds;
  Oinfo.n_options = nk;
	
  // initialize Ampl Solver Library
  asl = (ASL_pfgh *)ASL_alloc(ASL_read_pfgh);
  asl->i.want_xpi0_ = 3;

  stub = getstub(&nlFileName, &Oinfo);
	
  if (!stub) { // Couldn't find the name of the stub file
	usage_ASL(&Oinfo, 1);
	exit(0);
  }
	
  if (getopts(nlFileName, &Oinfo)){
	printf("Reach optimal solution at the initial point from input file. \n");	
	exit(0);
  }

  // Get dimension information and .nl file
  nlFile = jac0dim(stub, (fint)strlen(stub));
  
  // Set up ASL structure to store information
  #define Int(n) (int *)M1alloc((n)*sizeof(int));
  #define Double(n) (real *)M1alloc((n)*sizeof(real));

  asl->i.LUv_ 	= Double(asl->i.n_var_); //variable lower bound
  asl->i.Uvx_ 	= Double(asl->i.n_var_); //variable upper bound
  asl->i.LUrhs_ = Double(asl->i.n_con_); //constraint lower bound 
  asl->i.Urhsx_ = Double(asl->i.n_con_); //constraint upper bound
	
  asl->i.X0_ 		= Double(asl->i.n_var_);
  asl->i.havex0_ 	= new char[asl->i.n_var_];
  asl->i.pi0_ 		= Double(asl->i.n_con_);
  asl->i.Fortran_ 	= 0;

  if(amplSuffix)
	amplSuffix->InitSpaceForSuffixes(asl);

  int err3 = pfgh_read(nlFile, 0);
  if (err3!=0){
	printf("pfgh_read returns err=%d\n",err3);
	exit(1);
  }

  /* congrad should fill in according to goff [hooking, p.12] */
  asl->i.congrd_mode = 2;

  for(i=0;i<asl->i.n_con_;i++){
	if (fabs(asl->i.LUrhs_[i]-asl->i.Urhsx_[i])>1e-8){
	  nslacks++;
	  if (asl->i.LUrhs_[i]<-1e20 && asl->i.Urhsx_[i]>1e20){
	  	printf("Bounds on constraint %d: %f %f => free\n",i,asl->i.LUrhs_[i], asl->i.Urhsx_[i]);
	  	printf("NOT SUPPORTED!\n");
	  	exit(0);
	  }
	}
  }

  didLoad=true;
  locASL = asl;
  return locASL;	
}	





AmplSuffix::AmplSuffix()
	:suftab_ (NULL)
{}

AmplSuffix::~AmplSuffix()
{
  if (suftab_) {
	int n = (int)Suf_Idx.size();
	for (int i=0; i<n; i++) {
	  delete [] suftab_[i].name;
	  suftab_[i].name = NULL;
	}
  }
  delete [] suftab_;
  suftab_ = NULL;
}

void AmplSuffix::InitSpaceForSuffixes(ASL_pfgh* asl_)
{
  ASL_pfgh* asl = asl_;
  
  int n = (int) Suf_Idx.size();
  suftab_ = new SufDecl[n];
  
  for (int i=0; i<n; i++) {
	int id_len = (int)strlen(Suf_Idx[i].c_str());
	suftab_[i].name = new char[id_len + 1];
	strcpy(suftab_[i].name, Suf_Idx[i].c_str());

	suftab_[i].table = 0;

	if (Suf_Type[i] == Suffix_Var) {
	  suftab_[i].kind = ASL_Sufkind_var;
	}
	else if (Suf_Type[i]  == Suffix_Con) {
	  suftab_[i].kind = ASL_Sufkind_con;
	}
	else {
	  assert(false && "Unknown suffix source");
	}

	if (Suf_NumType[i] == Suffix_Double) {
	  suftab_[i].kind = suftab_[i].kind | ASL_Sufkind_real;
	}

	suftab_[i].nextra = 0;
  }

  suf_declare(suftab_, n);
}


int*
AmplSuffix::GetSuffixVal_Int(ASL_pfgh* asl_, std::string suffix_string, Suffix_Type type)
{
  ASL_pfgh* asl = asl_;

  int kind;
  if (type == Suffix_Var) {
	kind = ASL_Sufkind_var;
  }
  else if (type == Suffix_Con) {
	kind = ASL_Sufkind_con;
  }
  else {
	kind = 0;
  }
  SufDesc* dp = suf_get(suffix_string.c_str(), kind);
  return dp->u.i;
}


double*
AmplSuffix::GetSuffixVal_Double( ASL_pfgh* asl_, std::string suffix_string, Suffix_Type type)
{
  ASL_pfgh* asl = asl_;

  int kind;
  if (type == Suffix_Var) {
	kind = ASL_Sufkind_var;
  }
  else if (type == Suffix_Con) {
	kind = ASL_Sufkind_con;
  }
  else {
	kind = 0;
  }
  SufDesc* dp = suf_get(suffix_string.c_str(), kind);
  return dp->u.r;
}


