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

/**@file   nodesel_xyz.c
 * @brief  xyz node selector
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/nodesel_xyz.h"


#define NODESEL_NAME            "xyz"
#define NODESEL_DESC            "node selector template"
#define NODESEL_STDPRIORITY     0
#define NODESEL_MEMSAVEPRIORITY 0




/*
 * Data structures
 */

/* TODO: fill in the necessary node selector data */

/** node selector data */
struct SCIP_NodeselData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of node selector
 */

/* TODO: Implement all necessary node selector methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for node selector plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_NODESELCOPY(nodeselCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
 
   return SCIP_OKAY;
}
#else
#define nodeselCopyXyz NULL
#endif

/** destructor of node selector to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_NODESELFREE(nodeselFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselFreeXyz NULL
#endif


/** initialization method of node selector (called after problem was transformed) */
#if 0
static
SCIP_DECL_NODESELINIT(nodeselInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselInitXyz NULL
#endif


/** deinitialization method of node selector (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_NODESELEXIT(nodeselExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselExitXyz NULL
#endif


/** solving process initialization method of node selector (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_NODESELINITSOL(nodeselInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselInitsolXyz NULL
#endif


/** solving process deinitialization method of node selector (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_NODESELEXITSOL(nodeselExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselExitsolXyz NULL
#endif


/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** node comparison method of node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return 0;
}




/*
 * node selector specific interface methods
 */

/** creates the xyz node selector and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESELDATA* nodeseldata;

   /* create xyz node selector data */
   nodeseldata = NULL;
   /* TODO: (optional) create node selector specific data here */

   /* include node selector */
   SCIP_CALL( SCIPincludeNodesel(scip, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
         nodeselCopyXyz,
         nodeselFreeXyz, nodeselInitXyz, nodeselExitXyz, 
         nodeselInitsolXyz, nodeselExitsolXyz, nodeselSelectXyz, nodeselCompXyz,
         nodeseldata) );

   /* add xyz node selector parameters */
   /* TODO: (optional) add node selector specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
