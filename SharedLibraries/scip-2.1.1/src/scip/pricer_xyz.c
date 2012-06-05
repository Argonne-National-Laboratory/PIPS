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

/**@file   pricer_xyz.c
 * @brief  xyz variable pricer
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/pricer_xyz.h"


#define PRICER_NAME            "xyz"
#define PRICER_DESC            "variable pricer template"
#define PRICER_PRIORITY        0
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */




/*
 * Data structures
 */

/* TODO: fill in the necessary variable pricer data */

/** variable pricer data */
struct SCIP_PricerData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of variable pricer
 */

/* TODO: Implement all necessary variable pricer methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for pricer plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_PRICERCOPY(pricerCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
 
   return SCIP_OKAY;
}
#else
#define pricerCopyXyz NULL
#endif

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_PRICERFREE(pricerFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerFreeXyz NULL
#endif


/** initialization method of variable pricer (called after problem was transformed) */
#if 0
static
SCIP_DECL_PRICERINIT(pricerInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerInitXyz NULL
#endif


/** deinitialization method of variable pricer (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PRICEREXIT(pricerExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerExitXyz NULL
#endif


/** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_PRICERINITSOL(pricerInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerInitsolXyz NULL
#endif


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerExitsolXyz NULL
#endif


/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


#if 0
/** Farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerFarkasXyz NULL
#endif




/*
 * variable pricer specific interface methods
 */

/** creates the xyz variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICERDATA* pricerdata;

   /* create xyz variable pricer data */
   pricerdata = NULL;
   /* TODO: (optional) create variable pricer specific data here */

   /* include variable pricer */
   SCIP_CALL( SCIPincludePricer(scip, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerCopyXyz,
         pricerFreeXyz, pricerInitXyz, pricerExitXyz, 
         pricerInitsolXyz, pricerExitsolXyz, pricerRedcostXyz, pricerFarkasXyz,
         pricerdata) );

   /* add xyz variable pricer parameters */
   /* TODO: (optional) add variable pricer specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
