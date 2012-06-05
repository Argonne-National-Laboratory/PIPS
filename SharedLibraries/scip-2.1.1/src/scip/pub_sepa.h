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

/**@file   pub_sepa.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_SEPA_H__
#define __SCIP_PUB_SEPA_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_sepa.h"

#ifdef __cplusplus
extern "C" {
#endif

/** compares two separators w. r. to their priority */
extern
SCIP_DECL_SORTPTRCOMP(SCIPsepaComp);

/** gets user data of separator */
extern
SCIP_SEPADATA* SCIPsepaGetData(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** sets user data of separator; user has to free old data in advance! */
extern
void SCIPsepaSetData(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata            /**< new separator user data */
   );

/** gets name of separator */
extern
const char* SCIPsepaGetName(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets description of separator */
extern
const char* SCIPsepaGetDesc(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets priority of separator */
extern
int SCIPsepaGetPriority(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets frequency of separator */
extern
int SCIPsepaGetFreq(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** sets frequency of separator */
extern
void SCIPsepaSetFreq(
   SCIP_SEPA*            sepa,               /**< separator */
   int                   freq                /**< new frequency of separator */
   );

/** get maximal bound distance at which the separator is called */
extern
SCIP_Real SCIPsepaGetMaxbounddist(
   SCIP_SEPA*            sepa                /**< separator */
   );

/**< does the separator use a secondary SCIP instance? */
extern
SCIP_Bool SCIPsepaUsesSubscip(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets time in seconds used in this separator */
extern
SCIP_Real SCIPsepaGetTime(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the total number of times, the separator was called */
extern
SCIP_Longint SCIPsepaGetNCalls(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the number of times, the separator was called at the current node */
extern
int SCIPsepaGetNCallsAtNode(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets total number of times, the separator detected a cutoff */
extern
SCIP_Longint SCIPsepaGetNCutoffs(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the total number of cutting planes found by this separator */
extern
SCIP_Longint SCIPsepaGetNCutsFound(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets the number of cutting planes found by this separator at the current node */
extern
SCIP_Longint SCIPsepaGetNCutsFoundAtNode(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets total number of additional constraints added by this separator */
extern
SCIP_Longint SCIPsepaGetNConssFound(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** gets total number of domain reductions found by this separator */
extern
SCIP_Longint SCIPsepaGetNDomredsFound(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** should separator be delayed, if other separators found cuts? */
extern
SCIP_Bool SCIPsepaIsDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** was separation of the LP solution delayed at the last call? */
extern
SCIP_Bool SCIPsepaWasLPDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** was separation of the primal solution delayed at the last call? */
extern
SCIP_Bool SCIPsepaWasSolDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   );

/** is separator initialized? */
extern
SCIP_Bool SCIPsepaIsInitialized(
   SCIP_SEPA*            sepa                /**< separator */
   );

#ifdef __cplusplus
}
#endif

#endif
