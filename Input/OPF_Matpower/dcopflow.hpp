/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/
 

#ifndef OPFLOW_H
#define OPFLOW_H

#include "ps.h"

typedef struct _p_OPFLOW *OPFLOW;



class DCOPFLOW{
public:
  const DCPS   *ps;   /* Power system context */

  double obj; 		/* Objective function */

  double *gradobj; 	/* Gradient of the objective function */

  double *Xl; 		/* Lower bound on solution */
  double *Xu; 		/* Upper bound on solution */

  double *Gl; 		/* Lower bound on G */
  double *Gu; 		/* Upper bound on G */


  int Nvar; 		/* Number of variables */

  int Nconeq; 		/* Number of equality constraints */
  int Nconineq; 	/* Number of inequality constraints */
  int Ncon;     	/* Total number of constraints (equality + inequality) */

  int n; 			/* Number of variables */
  int m; 			/* Number of constraints */
  int nnz_jac_g; 	/* Number of nonzeros in the jacobian of the constraints */

  int refBusID;		/* reference bus id*/


  /* Lagrange multipliers */
  double *lambda_g;
  double *lambda_xl;
  double *lambda_xu;

  
  bool setupcalled; 		/* OPFLOWSetUp called? */
  bool setupcalled_part; 	/* OPFLOWSetUp_part called? */


  int Nvar_1st;
  int Ncon_1st;
  int Nconeq_1st;
  int Nconineq_1st;

  int *Nvar_2nd;
  int *Ncon_2nd;
  int *Nconeq_2nd;
  int *Nconineq_2nd;

  int nnz_jac_g_1st; 	/* Number of nonzeros in the jacobian of the constraints */
  int *nnz_jac_g_2nd; 	/* Number of nonzeros in the jacobian of the constraints */
  int *nnz_jac_g_Link; 	/* Number of nonzeros in the jacobian of the constraints */




  int *busMap_AllTo1st;
  int *varMap_AllTo2nd;


  int *numDummyVar;
  int *numDummyCon;


  int Nparts;

  DCOPFLOW(const DCPS* ps_);  
  ~DCOPFLOW();

  
  virtual void DCOPFLOWSetUp();
  virtual void SetVariableandConstraintBounds();
  virtual void ObjGradient_Lin(double *obj_coef);
  virtual void ObjGradient_Quad(int *irow, int *jcol, double *obj_quad);
  
  virtual int  GetJacNNZ();
  virtual void SetJacLocations(int *row, int *col);
  virtual void SetJacValues(double *values);


  virtual void DCOPFLOWSetUp_Partition();
  virtual void VarAndConBounds_1st_Partition(double *xl,double *xu,double *gl,double *gu);
  virtual void VarAndConBounds_2nd_Partition(const int scen, double *xl,double *xu,double *gl,double *gu);
  virtual void ObjGradient_Lin_1st_Partition(double *obj_coef);
  virtual void ObjGradient_Lin_2nd_Partition(const int scen, double *obj_coef);

  
  virtual int  GetHesNNZ_1st_Partition(){return 0;};
  virtual void ObjGradient_Quad_1st_Partition(int *irow, int *jcol, double *obj_quad);  
  virtual int  GetHesNNZ_2nd_Partition(const int scen);
  virtual void ObjGradient_Quad_2nd_Partition(const int scen,int *irow, int *jcol, double *obj_quad);

  virtual int  GetJacNNZ_1st_Partition();
  virtual int  GetJacNNZ_2nd_Partition(const int scen);
  virtual int  GetJacNNZ_Link_Partition(const int scen);

  virtual void GetJac_1st_Partition(int *row, int *col, double *ele);
  virtual void GetJac_2nd_Link_Partition(	const int scen, int *row, int *col, double *ele, 
  												int *row_link, int *col_link, double *ele_link);
  
};



#endif

