/* This file is part of OOPS.
 *
 * OOPS is (c) 2003-2009 Jacek Gondzio and Andreas Grothey, 
 *                       University of Edinburgh
 *
 * OOPS is distributed in a restricted form in the hope that it will be a useful
 * example of what can be done with SML, however it is NOT released under a free
 * software license.
 *
 * You may only redistribute this version of OOPS with a version of SML. You
 * may not link OOPS with code which is not part of SML licensed under the
 * LGPL v3.
 *
 * You may NOT modify, disassemble, or otherwise reverse engineer OOPS.
 *
 * OOPS is distributed WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef HOPDM_H
#define HOPDM_H

#include "oops/DenseVector.h"
#include "oops/Vector.h"
#include "oops/LogiVector.h"
#include "oops/Algebra.h"
#include "oops/GlobalOpt.h"

#define hopdm_opt_type HopdmOptions
#define hopdm_prt_type PrintOptions
#define primal_dual_pb PDProblem

class HopdmOptions {

 public:

  HopdmOptions();
  ~HopdmOptions();
  int readGlobalOptions();
  int use_start_point;
  int use_wsh;
  int use_presolve;
  int nb_beg_center_it;
  int nb_end_center_it;
  int ret_adv_center;
  int nb_ret_adv_center_it;
  double ret_adv_center_gap;
  double ret_target_mu;
  int ret_path;
  GlobalOpt *glopt;

};

enum PrintLevelValues {

  PRINT_NONE,
  PRINT_ITER,
  PRINT_INFO,
  PRINT_VERBOSE,
  LAST_LEVEL
};


class PrintOptions {

 public:

  PrintOptions(const int PrintLevel = PRINT_ITER);

   int PrtFact;
   int PrtIRSlv;
   int PrtInit;
   int PrtDTheta;
   int PrtPushV;
   int PrtResid;
   int PrtComPr;
   int Prthopdm;
   int PrthoDir;
   int PrtMaxS;
   int PrtMakeS;
   int PrtTT;

 private:

  void setPrintLevel(const int PrintLevel);

};

typedef struct {

  double val;
  double dVal;
  int ifail;
  int iters; 
  int found_adv_center;
  double err_b, err_c, err_u;
  double gap;

} hopdm_ret;

class PDProblem {

 public:

  PDProblem(Algebra *AlgAug = NULL,
            Vector *b = NULL, Vector *c = NULL, Vector *u = NULL,
            Vector *x = NULL, Vector *y = NULL, Vector *z = NULL);

  ~PDProblem();

  Algebra *AlgAug;
  Vector *b;
  Vector *c;
  Vector *u;
  Vector *l;
      Vector *x, *cx;
      Vector *y, *cy;
      Vector *z, *cz;
      Vector *s, *cs;
      Vector *w, *cw;
  double adv_mu;
  Vector **path;
  Vector *target;
      int is_feas;
      double obj_fact;
      int set_by_constr; 

};

typedef struct
{
  Algebra *AlgAug;
  Vector *vtheta;
  Vector *vthetay;
  Vector *vpdRegTh;
  Vector *vpdRegx;
  Vector *vpdRegy;

} InverseRep;

InverseRep *
NewInverseRep(Algebra *AlgAug, Vector *vtheta, Vector *vthetay, 
	      Vector *pdRegTh, Vector *pdRegx, Vector *pdRegy);

typedef struct
{
  Vector *x;
  Vector *y;
  Vector *z;
  Vector *s;
  Vector *w;
  Vector *xy;

} PDPoint;

void
MaxStep(FILE *out, PDPoint *pdPoint, PDPoint *pdDir,
	double *alphaX, double *alphaZ, double *alphaS, double *alphaW,
	Algebra *AlgA, LogiVector *vwhere_u, LogiVector *vwhere_l, Vector *vstep,
        PrintOptions *Prt);

void
CompPDRes(FILE *out, Algebra * AlgA, Algebra *AlgQ,
	  Vector *vb, Vector *vc, Vector *vu, PDPoint *pdPoint,
	  Vector *vxib, Vector *vxic, Vector *vxiu,
	  double *err_b, double *err_c, double *err_u, 
          LogiVector *vwhere_u, PrintOptions *Prt);

int
DefTheta(FILE *out, Vector *vx, Vector *vs, Vector *vz, Vector *vw,
         Vector *vtheta, LogiVector *vwhere_u, LogiVector *vwhere_l,
         PrintOptions *Prt);

void
pdFactor (FILE *out, Algebra *AlgAug, Vector *vtheta, Vector *vthetay, 
	  Vector *vpdRegTh, Vector *vpdRegx, Vector *vpdRegy, 
          PrintOptions *Prt);

void
IterRefSolveNew (FILE *out, InverseRep *ivr, int *Alarm,
		 Vector *vrhs_x, Vector *vrhs_y, 
		 Vector *vdel_xy, Vector *vdel_x, Vector *vdel_y,
		 PrintOptions *Prt, Vector *vNwrk3, Vector *vMwrk3,
		 double ireps);

void
HopdmDir(FILE *out, InverseRep *ivr, const int iDir, const int onlyCent,
	 const double oldbarr, const double barr, double AlphaP, double AlphaD,
	 Vector *vxib, Vector *vxic, Vector *vxiu,
	 Vector *vXrhs_x, Vector *target,
	  PDPoint *pdPoint, PDPoint *pdPredDir, PDPoint *pdNewDir,
         LogiVector *vwhere_u, LogiVector *vwhere_l, PrintOptions *Prt,
	 Vector *vNw1, Vector *vMw1, const double ireps);

hopdm_ret*
hopdm(FILE *out, PDProblem *P, HopdmOptions *opt, PrintOptions *Prt);

hopdm_ret*
hopdm_std(PDProblem *P, HopdmOptions *opt, PrintOptions *Prt);

hopdm_ret*
SolveOops(PDProblem* pb);

PDProblem*
NewPDProblem(Algebra *AlgAug, Vector *b, Vector *c, Vector *u,
	     Vector *x, Vector *y, Vector *z);

void
FreePDProblem(PDProblem *P);

PDPoint*
NewPDPoint(Vector *x, Vector *y, Vector *z, Vector *s, Vector *w, Vector *xy);

void
FreePDPoint(PDPoint *P);

HopdmOptions*
NewHopdmOptions(void);

void
FreeHopdmOptions(HopdmOptions *opt);

PrintOptions*
NewHopdmPrt(const int level);

#endif
