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

/**@file   relax_xyz.c
 * @brief  xyz relaxator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/relax_xyz.h"


#define RELAX_NAME             "xyz"
#define RELAX_DESC             "relaxator template"
#define RELAX_PRIORITY         0
#define RELAX_FREQ             1




/*
 * Data structures
 */

/* TODO: fill in the necessary relaxator data */

/** relaxator data */
struct SCIP_RelaxData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of relaxator
 */

/* TODO: Implement all necessary relaxator methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for relaxator plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_RELAXCOPY(relaxCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxCopyXyz NULL
#endif

/** destructor of relaxator to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_RELAXFREE(relaxFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxFreeXyz NULL
#endif


/** initialization method of relaxator (called after problem was transformed) */
#if 0
static
SCIP_DECL_RELAXINIT(relaxInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxInitXyz NULL
#endif


/** deinitialization method of relaxator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_RELAXEXIT(relaxExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxExitXyz NULL
#endif


/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_RELAXINITSOL(relaxInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxInitsolXyz NULL
#endif


/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_RELAXEXITSOL(relaxExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define relaxExitsolXyz NULL
#endif


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz relaxator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}





/*
 * relaxator specific interface methods
 */

/** creates the xyz relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;

   /* create xyz relaxator data */
   relaxdata = NULL;
   /* TODO: (optional) create relaxator specific data here */

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelax(scip, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ, 
         relaxCopyXyz,
         relaxFreeXyz, relaxInitXyz, relaxExitXyz, relaxInitsolXyz, relaxExitsolXyz, relaxExecXyz,
         relaxdata) );

   /* add xyz relaxator parameters */
   /* TODO: (optional) add relaxator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
