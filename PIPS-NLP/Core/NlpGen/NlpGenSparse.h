/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */
 
 /* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef SPSEQNLPGENFACTORY
#define SPSEQNLPGENFACTORY

#include "NlpGen.h"
class NlpGenData;
class NlpGenVars;
class NlpInfo;

class NlpGenSparse : public NlpGen{
protected:
  int nnzQ;
  int nnzA;
  int nnzC;
public:

  NlpGenSparse( int nx_, int my_, int mz_,
		  int nnzQ_, int nnzA_, int nnzC_ )
    : NlpGen( nx_, my_, mz_ ),
      nnzQ(nnzQ_), nnzA(nnzA_), nnzC(nnzC_)
  {}

  Data  * makeData();
  Data  * makeData( double    c[],
			 int    krowQ[],  int  jcolQ[],  double dQ[],
			 double  xlow[],  char ixlow[],
			 double  xupp[],  char ixupp[],
			 int    krowA[],  int  jcolA[],  double dA[],
			 double     b[],
			 int    krowC[],  int  jcolC[],  double dC[],
			 double  clow[],  char iclow[],
			 double  cupp[],  char icupp[] );

  Data   *
  copyDataFromSparseTriple( double c[],
			    int irowQ[], int nnzQ,  int jcolQ[],  double dQ[],
			    double xlow[],  char ixlow[],
			    double xupp[],  char ixupp[],
			    int irowA[], int nnzA,  int jcolA[],  double dA[],
			    double   bA[],
			    int irowC[],  int nnzC,  int jcolC[], double dC[],
			    double clow[], char iclow[],
			    double cupp[], char icupp[], int *rowMap,NlpInfo * updateNlp);

  virtual void joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			OoqpVector& rhs2_in, OoqpVector& rhs3_in);

  virtual void separateVars( OoqpVector& x_in, OoqpVector& y_in,
			     OoqpVector& z_in, OoqpVector& vars_in);




  Data   *
  copyDataFromSparseTriple( double c[],
			  int irowQ[], int nnzQ,  int jcolQ[],	double dQ[],
			  double xlow[],  char ixlow[],
			  double xupp[],  char ixupp[],
			  int irowA[], int nnzA,  int jcolA[],	double dA[],
			  double   bA[],
			  int irowC[],	int nnzC,  int jcolC[], double dC[],
			  double clow[], char iclow[],
			  double cupp[], char icupp[], int *rowMap,
			  int nxL,int nxU,int nsL,int nsU,NlpInfo * updateNlp);


  virtual void joinRHSXSYZ( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			OoqpVector& rhs2_in, OoqpVector& rhs3_in, OoqpVector& rhs4_in );

  virtual void separateVarsXSYZ( OoqpVector& x_in, OoqpVector& s_in, 
  			OoqpVector& y_in, OoqpVector& z_in, OoqpVector& vars_in);  

  virtual void copyXSYZ_fromArray( OoqpVector& vec_xsyz, double* array_in, const int nb_col);
  virtual void copyXSYZ_toArray( OoqpVector& vec_xsyz, double* array_in, const int nb_col);  


};

#endif
