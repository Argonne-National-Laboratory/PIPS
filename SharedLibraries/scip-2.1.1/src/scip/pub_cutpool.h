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

/**@file   pub_cutpool.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for storing cuts in a cut pool
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_CUTPOOL_H__
#define __SCIP_PUB_CUTPOOL_H__


#include "scip/def.h"
#include "scip/type_cutpool.h"

#ifdef __cplusplus
extern "C" {
#endif

/** gets the row of the cut */
extern
SCIP_ROW* SCIPcutGetRow(
   SCIP_CUT*             cut                 /**< cut */
   );

/** gets the age of the cut: the number of consecutive cut pool separation rounds where the cut was neither in the LP nor violated */
extern
int SCIPcutGetAge(
   SCIP_CUT*             cut                 /**< cut */
   );

/** gets array of cuts in the cut pool */
extern
SCIP_CUT** SCIPcutpoolGetCuts(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );

/** get number of cuts in the cut pool */
extern
int SCIPcutpoolGetNCuts(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );

/** get maximum number of cuts that were stored in the cut pool at the same time */
extern
int SCIPcutpoolGetMaxNCuts(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );

/** gets time in seconds used for separating cuts from the pool */
extern
SCIP_Real SCIPcutpoolGetTime(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );

/** get number of times, the cut pool was separated */
extern
SCIP_Longint SCIPcutpoolGetNCalls(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );

/** get total number of cuts that were separated from the cut pool */
extern
SCIP_Longint SCIPcutpoolGetNCutsFound(
   SCIP_CUTPOOL*         cutpool             /**< cut pool */
   );

#ifdef __cplusplus
}
#endif

#endif
