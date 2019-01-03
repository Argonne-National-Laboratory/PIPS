/* global constants (symbol dimensions, symbol layout, etc.)
 * that might ultimately come from global
 * created: Steve Dirkse, July 2007
 */

#include "gclgms.h"

const char *gmsGdxTypeText[GMS_DT_MAX] =
  {"Set","Parameter","Variable","Equation","Alias"};
const char *gmsVarTypeText[GMS_VARTYPE_MAX] = {"Unknown","Binary","Integer","Positive","Negative","Free","Sos1","Sos2","Semicont","Semiint"};
const char *gmsValTypeText[GMS_VAL_MAX] = {".l",".m",".lo",".up",".scale"};
const char *gmsSVText[GMS_SVIDX_MAX] = {"UNdef","NA","+Inf","-Inf","Eps","0","AcroN"};

const char * const rcStat[STAT_MAX] = {"     "," NOPT","INFES","UNBND"," ****","UNKNW","REDEF","DEPND","REDIR"};

const char * const solveStatusTxt[SS_MAX] = {"NA                      ",
                                             "Normal Completion       ",
                                             "Iteration Interrupt     ",
                                             "Resource Interrupt      ",
                                             "Terminated By Solver    ",
                                             "Evaluation Interrupt    ",
                                             "Capability Problems     ",
                                             "Licensing Problems      ",
                                             "User Interrupt          ",
                                             "Setup Failure           ",
                                             "Solver Failure          ",
                                             "Internal Solver Failure ",
                                             "Solve Processing Skipped",
                                             "System Failure          "};

const char * const modelStatusTxt[MS_MAX] = {"NA                      ",
                                             "Optimal                 ",
                                             "Locally Optimal         ",
                                             "Unbounded               ",
                                             "Infeasible              ",
                                             "Locally Infeasible      ",
                                             "Intermediate Infeasible ",
                                             "Feasible Solution       ",
                                             "Integer Solution        ",
                                             "Intermediate Non-Integer",
                                             "Integer Infeasible      ",
                                             "Licensing Problem       ",
                                             "Error Unknown           ",
                                             "Error No Solution       ",
                                             "No Solution Returned    ",
                                             "Solved Unique           ",
                                             "Solved                  ",
                                             "Solved Singular         ",
                                             "Unbounded - No Solution ",
                                             "Infeasible - No Solution"};


const double gmsDefRecVar[GMS_VARTYPE_MAX][GMS_VAL_MAX] = {
 /* .l   .m           .lo         .ub  .scale */
  { 0.0, 0.0,         0.0,         0.0, 1.0},    /* unknown */
  { 0.0, 0.0,         0.0,         1.0, 1.0},    /* binary */
  { 0.0, 0.0,         0.0,       100.0, 1.0},    /* integer */
  { 0.0, 0.0,         0.0, GMS_SV_PINF, 1.0},    /* positive */
  { 0.0, 0.0, GMS_SV_MINF,         0.0, 1.0},    /* negative */
  { 0.0, 0.0, GMS_SV_MINF, GMS_SV_PINF, 1.0},    /* free */
  { 0.0, 0.0,         0.0, GMS_SV_PINF, 1.0},    /* sos1 */
  { 0.0, 0.0,         0.0, GMS_SV_PINF, 1.0},    /* sos2 */
  { 0.0, 0.0,         1.0, GMS_SV_PINF, 1.0},    /* semicont */
  { 0.0, 0.0,         1.0,       100.0, 1.0}     /* semiint */
};

const double gmsDefRecEqu[GMS_EQUTYPE_MAX][GMS_VAL_MAX] = {
 /* .l   .m           .lo         .ub  .scale */
  { 0.0, 0.0,         0.0,         0.0, 1.0},    /* =e= */
  { 0.0, 0.0,         0.0, GMS_SV_PINF, 1.0},    /* =g= */
  { 0.0, 0.0, GMS_SV_MINF,         0.0, 1.0},    /* =l= */
  { 0.0, 0.0, GMS_SV_MINF, GMS_SV_PINF, 1.0},    /* =n= */
  { 0.0, 0.0,         0.0,         0.0, 1.0},    /* =x= */
  { 0.0, 0.0,         0.0, GMS_SV_PINF, 1.0}     /* =c= */
};

/* extract an equation type in [GMS_EQUTYPE_E,GMS_EQUTYPE_MAX)
 * from the userInfo value stored for an equation symbol,
 * or a negative value if the equation type was not stored
 */
int
gmsFixEquType (int userInfo)
{
  int equType = userInfo - GMS_EQU_USERINFO_BASE;

  if (equType >= GMS_EQUTYPE_MAX)
    equType = -1;

  return equType;
} /* gmsFixEquType */

/* extract a variable type in [GMS_VARTYPE_UNKNOWN,GMS_VARTYPE_MAX)
 * from the userInfo value stored for a variable symbol,
 * or a negative value if the variable type was not stored
 */
int
gmsFixVarType (int userInfo)
{
  int varType = userInfo;

  if (varType >= GMS_VARTYPE_MAX)
    varType = -1;

  return varType;
} /* gmsFixVarType */
