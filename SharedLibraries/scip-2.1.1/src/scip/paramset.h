/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   paramset.h
 * @brief  internal methods for handling parameter settings
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PARAMSET_H__
#define __SCIP_PARAMSET_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_paramset.h"
#include "scip/pub_paramset.h"
#include "scip/pub_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates parameter set */
extern
SCIP_RETCODE SCIPparamsetCreate(
   SCIP_PARAMSET**       paramset,           /**< pointer to store the parameter set */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** frees parameter set */
extern
void SCIPparamsetFree(
   SCIP_PARAMSET**       paramset,           /**< pointer to the parameter set */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** creates a bool parameter, sets it to its default value, and adds it to the parameter set */
extern
SCIP_RETCODE SCIPparamsetAddBool(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Bool             defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a int parameter, sets it to its default value, and adds it to the parameter set */
extern
SCIP_RETCODE SCIPparamsetAddInt(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   int*                  valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   int                   defaultvalue,       /**< default value of the parameter */
   int                   minvalue,           /**< minimum value for parameter */
   int                   maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a SCIP_Longint parameter, sets it to its default value, and adds it to the parameter set */
extern
SCIP_RETCODE SCIPparamsetAddLongint(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Longint          defaultvalue,       /**< default value of the parameter */
   SCIP_Longint          minvalue,           /**< minimum value for parameter */
   SCIP_Longint          maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a SCIP_Real parameter, sets it to its default value, and adds it to the parameter set */
extern
SCIP_RETCODE SCIPparamsetAddReal(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Real             defaultvalue,       /**< default value of the parameter */
   SCIP_Real             minvalue,           /**< minimum value for parameter */
   SCIP_Real             maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
extern
SCIP_RETCODE SCIPparamsetAddChar(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char*                 valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   char                  defaultvalue,       /**< default value of the parameter */
   const char*           allowedvalues,      /**< array with possible parameter values, or NULL if not restricted */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
extern
SCIP_RETCODE SCIPparamsetAddString(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char**                valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   const char*           defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   );

/** gets the value of an existing SCIP_Bool parameter */
extern
SCIP_RETCODE SCIPparamsetGetBool(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing int parameter */
extern
SCIP_RETCODE SCIPparamsetGetInt(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           name,               /**< name of the parameter */
   int*                  value               /**< pointer to store the parameter */
   );

/** gets the value of an existing SCIP_Longint parameter */
extern
SCIP_RETCODE SCIPparamsetGetLongint(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint*         value               /**< pointer to store the parameter */
   );

/** gets the value of an existing SCIP_Real parameter */
extern
SCIP_RETCODE SCIPparamsetGetReal(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           name,               /**< name of the parameter */
   SCIP_Real*            value               /**< pointer to store the parameter */
   );

/** gets the value of an existing char parameter */
extern
SCIP_RETCODE SCIPparamsetGetChar(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           name,               /**< name of the parameter */
   char*                 value               /**< pointer to store the parameter */
   );

/** gets the value of an existing string parameter */
extern
SCIP_RETCODE SCIPparamsetGetString(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           name,               /**< name of the parameter */
   char**                value               /**< pointer to store the parameter */
   );

/** changes the value of an existing parameter */
extern
SCIP_RETCODE SCIPparamsetSet(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   void*                 value               /**< new value of the parameter */
   );

/** changes the value of an existing SCIP_Bool parameter */
extern
SCIP_RETCODE SCIPparamsetSetBool(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool             value               /**< new value of the parameter */
   );

/** changes the value of an existing int parameter */
extern
SCIP_RETCODE SCIPparamsetSetInt(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   int                   value               /**< new value of the parameter */
   );

/** changes the value of an existing SCIP_Longint parameter */
extern
SCIP_RETCODE SCIPparamsetSetLongint(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint          value               /**< new value of the parameter */
   );

/** changes the value of an existing SCIP_Real parameter */
extern
SCIP_RETCODE SCIPparamsetSetReal(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   SCIP_Real             value               /**< new value of the parameter */
   );

/** changes the value of an existing char parameter */
extern
SCIP_RETCODE SCIPparamsetSetChar(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   char                  value               /**< new value of the parameter */
   );

/** changes the value of an existing string parameter */
extern
SCIP_RETCODE SCIPparamsetSetString(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of the parameter */
   const char*           value               /**< new value of the parameter */
   );

/** reads parameters from a file */
SCIP_RETCODE SCIPparamsetRead(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           filename            /**< file name */
   );

/** writes all parameters in the parameter set to a file */
SCIP_RETCODE SCIPparamsetWrite(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   const char*           filename,           /**< file name, or NULL for stdout */
   SCIP_Bool             comments,           /**< should parameter descriptions be written as comments? */
   SCIP_Bool             onlychanged         /**< should only the parameters been written, that are changed from default? */
   );

/** installs default values for all parameters */
extern
SCIP_RETCODE SCIPparamsetSetToDefaults(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip                /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   );

/** installs default value for a single parameter */
extern
SCIP_RETCODE SCIPparamsetSetToDefault(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure, or NULL if paramchgd method should not be called */   
   const char*           paramname           /**< name of the parameter */
   );

/** sets parameters to 
 *  - SCIP_PARAMSETTING_DEFAULT to use default values (see also SCIPparamsetSetToDefault())
 *  - SCIP_PARAMSETTING_COUNTER to get feasible and "fast" counting process
 *  - SCIP_PARAMSETTING_CPSOLVER to get CP like search (e.g. no LP relaxation)
 *  - SCIP_PARAMSETTING_EASYCIP to solve easy problems fast
 *  - SCIP_PARAMSETTING_FEASIBILITY to detect feasibility fast 
 *  - SCIP_PARAMSETTING_HARDLP to be capable to handle hard LPs
 *  - SCIP_PARAMSETTING_OPTIMALITY to prove optimality fast
 */
extern
SCIP_RETCODE SCIPparamsetSetEmphasis(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMEMPHASIS    paramemphasis,      /**< parameter emphasis */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** sets parameters to deactivate separators and heuristics that use auxiliary SCIP instances; should be called for
 *  auxiliary SCIP instances to avoid recursion
 */
extern
SCIP_RETCODE SCIPparamsetSetToSubscipsOff(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** sets heuristic parameters values to 
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all heuristic parameters 
 *  - SCIP_PARAMSETTING_FAST such that the time spend for heuristic is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the heuristic are called more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all heuristics
 */
extern
SCIP_RETCODE SCIPparamsetSetHeuristics(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** sets presolving parameters to 
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all presolving parameters 
 *  - SCIP_PARAMSETTING_FAST such that the time spend for presolving is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the presolving is more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all presolving
 */
extern
SCIP_RETCODE SCIPparamsetSetPresolving(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** sets separating parameters to 
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all separating parameters 
 *  - SCIP_PARAMSETTING_FAST such that the time spend for separating is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the separating is done more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all separating
 */
extern
SCIP_RETCODE SCIPparamsetSetSeparating(
   SCIP_PARAMSET*        paramset,           /**< parameter set */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   );

/** returns the array of parameters */
extern
SCIP_PARAM** SCIPparamsetGetParams(
   SCIP_PARAMSET*        paramset            /**< parameter set */
   );

/** returns the number of parameters in the parameter set */
extern
int SCIPparamsetGetNParams(
   SCIP_PARAMSET*        paramset            /**< parameter set */
   );

/** copies all parameter values of the source parameter set to the corresponding parameters in the target set */
extern
SCIP_RETCODE SCIPparamsetCopyParams(
 SCIP_PARAMSET*          sourceparamset,     /**< source parameter set */
 SCIP_PARAMSET*          targetparamset,     /**< target parameter set */
 SCIP*                   targetscip          /**< target SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
