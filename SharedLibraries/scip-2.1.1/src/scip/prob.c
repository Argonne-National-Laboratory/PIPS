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

/**@file   prob.c
 * @brief  Methods and datastructures for storing and manipulating the main problem
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/message.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/pub_misc.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/primal.h"
#include "scip/tree.h"
#include "scip/branch.h"
#include "scip/cons.h"


#define OBJSCALE_MAXDNOM          1000000LL  /**< maximal denominator in objective integral scaling */
#define OBJSCALE_MAXSCALE         1000000.0  /**< maximal scalar to reach objective integrality */
#define OBJSCALE_MAXFINALSCALE       1000.0  /**< maximal final value to apply as scaling */



/*
 * dymanic memory arrays
 */

/** resizes vars array to be able to store at least num entries */
static
SCIP_RETCODE probEnsureVarsMem(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(prob != NULL);
   assert(set != NULL);

   if( num > prob->varssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&prob->vars, newsize) );
      prob->varssize = newsize;
   }
   assert(num <= prob->varssize);

   return SCIP_OKAY;
}

/** resizes fixedvars array to be able to store at least num entries */
static
SCIP_RETCODE probEnsureFixedvarsMem(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(prob != NULL);
   assert(set != NULL);

   if( num > prob->fixedvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&prob->fixedvars, newsize) );
      prob->fixedvarssize = newsize;
   }
   assert(num <= prob->fixedvarssize);

   return SCIP_OKAY;
}

/** resizes deletedvars array to be able to store at least num entries */
static
SCIP_RETCODE probEnsureDeletedvarsMem(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(prob != NULL);
   assert(set != NULL);

   if( num > prob->deletedvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&prob->deletedvars, newsize) );
      prob->deletedvarssize = newsize;
   }
   assert(num <= prob->deletedvarssize);

   return SCIP_OKAY;
}

/** resizes conss array to be able to store at least num entries */
static
SCIP_RETCODE probEnsureConssMem(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(prob != NULL);
   assert(set != NULL);

   if( num > prob->consssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&prob->conss, newsize) );
      prob->consssize = newsize;
   }
   assert(num <= prob->consssize);

   return SCIP_OKAY;
}

/** returns whether the constraint has a name */
static
SCIP_Bool consHasName(
   SCIP_CONS*            cons                /**< constraint */
   )
{
   const char* name;

   name = SCIPconsGetName(cons);

   return (name != NULL && name[0] != '\0');
}

/** returns whether the variable has a name */
static
SCIP_Bool varHasName(
   SCIP_VAR*             var                 /**< variable */
   )
{
   const char* name;

   name = SCIPvarGetName(var);

   return (name != NULL && name[0] != '\0');
}



/*
 * problem creation
 */

/** creates problem data structure by copying the source problem; 
 *  If the problem type requires the use of variable pricers, these pricers should be activated with calls
 *  to SCIPactivatePricer(). These pricers are automatically deactivated, when the problem is freed.
 */
SCIP_RETCODE SCIPprobCopy(
   SCIP_PROB**           prob,               /**< pointer to problem data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< problem name */
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP_PROB*            sourceprob,         /**< source problem structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints */
   SCIP_Bool             global              /**< create a global or a local copy? */
   )
{
   SCIP_PROBDATA* targetdata;
   SCIP_RESULT result;

   assert(prob != NULL);
   assert(set != NULL);
   assert(blkmem != NULL);
   assert(sourcescip != NULL);
   assert(sourceprob != NULL);
   assert(varmap != NULL);
   assert(consmap != NULL);

   result = SCIP_DIDNOTRUN;
   targetdata = NULL;

   /* create problem and initialize callbacks with NULL */
   SCIP_CALL( SCIPprobCreate(prob, blkmem, set, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL, FALSE) );
   
   /* call user copy callback method */
   if( sourceprob->probdata != NULL && sourceprob->probcopy != NULL )
   {
      SCIP_CALL( sourceprob->probcopy(set->scip, sourcescip, sourceprob->probdata, varmap, consmap, &targetdata, global, &result) );

      /* evaluate result */
      if( result != SCIP_DIDNOTRUN && result != SCIP_SUCCESS )
      {
         SCIPerrorMessage("probdata copying method returned invalid result <%d>\n", result);
         return SCIP_INVALIDRESULT;
      }

      assert(targetdata == NULL || result == SCIP_SUCCESS);
   }

   /* if copying was successful, add data and callbacks */
   if( result == SCIP_SUCCESS )
   {
      assert( targetdata != NULL );
      (*prob)->probdelorig = sourceprob->probdelorig;
      (*prob)->probtrans = sourceprob->probtrans;
      (*prob)->probdeltrans = sourceprob->probdeltrans;
      (*prob)->probinitsol = sourceprob->probinitsol;
      (*prob)->probexitsol = sourceprob->probexitsol;
      (*prob)->probcopy = sourceprob->probcopy;
      (*prob)->probdata = targetdata;
   }

   return SCIP_OKAY;
}

/** creates problem data structure
 *  If the problem type requires the use of variable pricers, these pricers should be activated with calls
 *  to SCIPactivatePricer(). These pricers are automatically deactivated, when the problem is freed.
 */
SCIP_RETCODE SCIPprobCreate(
   SCIP_PROB**           prob,               /**< pointer to problem data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< problem name */
   SCIP_DECL_PROBDELORIG ((*probdelorig)),   /**< frees user data of original problem */
   SCIP_DECL_PROBTRANS   ((*probtrans)),     /**< creates user data of transformed problem by transforming original user data */
   SCIP_DECL_PROBDELTRANS((*probdeltrans)),  /**< frees user data of transformed problem */
   SCIP_DECL_PROBINITSOL ((*probinitsol)),   /**< solving process initialization method of transformed data */
   SCIP_DECL_PROBEXITSOL ((*probexitsol)),   /**< solving process deinitialization method of transformed data */
   SCIP_DECL_PROBCOPY    ((*probcopy)),      /**< copies user data if you want to copy it to a subscip, or NULL */
   SCIP_PROBDATA*        probdata,           /**< user problem data set by the reader */
   SCIP_Bool             transformed         /**< is this the transformed problem? */
   )
{
   assert(prob != NULL);

   SCIP_ALLOC( BMSallocMemory(prob) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*prob)->name, name, strlen(name)+1) );

   (*prob)->probdata = probdata;
   (*prob)->probcopy = probcopy;
   (*prob)->probdelorig = probdelorig;
   (*prob)->probtrans = probtrans;
   (*prob)->probdeltrans = probdeltrans;
   (*prob)->probinitsol = probinitsol;
   (*prob)->probexitsol = probexitsol;
   if( set->misc_usevartable )
   {
      SCIP_CALL( SCIPhashtableCreate(&(*prob)->varnames, blkmem, 
            (set->misc_usesmalltables ? SCIP_HASHSIZE_NAMES_SMALL : SCIP_HASHSIZE_NAMES),
            SCIPhashGetKeyVar, SCIPhashKeyEqString, SCIPhashKeyValString, NULL) );
   }
   else
      (*prob)->varnames = NULL;
   (*prob)->vars = NULL;
   (*prob)->varssize = 0;
   (*prob)->nvars = 0;
   (*prob)->nbinvars = 0;
   (*prob)->nintvars = 0;
   (*prob)->nimplvars = 0;
   (*prob)->ncontvars = 0;
   (*prob)->ncolvars = 0;
   (*prob)->fixedvars = NULL;
   (*prob)->fixedvarssize = 0;
   (*prob)->nfixedvars = 0;
   (*prob)->deletedvars = NULL;
   (*prob)->deletedvarssize = 0;
   (*prob)->ndeletedvars = 0;
   if( set->misc_useconstable )
   {
      SCIP_CALL( SCIPhashtableCreate(&(*prob)->consnames, blkmem,
            (set->misc_usesmalltables ? SCIP_HASHSIZE_NAMES_SMALL : SCIP_HASHSIZE_NAMES),
            SCIPhashGetKeyCons, SCIPhashKeyEqString, SCIPhashKeyValString, NULL) );
   }
   else
      (*prob)->consnames = NULL;
   (*prob)->conss = NULL;
   (*prob)->consssize = 0;
   (*prob)->nconss = 0;
   (*prob)->maxnconss = 0;
   (*prob)->startnvars = 0;
   (*prob)->startnconss = 0;
   (*prob)->objsense = SCIP_OBJSENSE_MINIMIZE;
   (*prob)->objoffset = 0.0;
   (*prob)->objscale = 1.0;
   (*prob)->objlim = SCIP_INVALID;
   (*prob)->objisintegral = FALSE;
   (*prob)->transformed = transformed;

   return SCIP_OKAY;
}

/** frees problem data structure */
SCIP_RETCODE SCIPprobFree(
   SCIP_PROB**           prob,               /**< pointer to problem data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp                  /**< current LP data (or NULL, if it's the original problem) */
   )
{
   int v;

   assert(prob != NULL);
   assert(*prob != NULL);
   assert(set != NULL);
   
   /* remove all constraints from the problem */
   while( (*prob)->nconss > 0 )
   {
      assert((*prob)->conss != NULL);
      SCIP_CALL( SCIPprobDelCons(*prob, blkmem, set, stat, (*prob)->conss[0]) );
   }

   if( (*prob)->transformed )
   {
      int h;

      /* unlock variables for all constraint handlers that don't need constraints */
      for( h = 0; h < set->nconshdlrs; ++h )
      {
         if( !SCIPconshdlrNeedsCons(set->conshdlrs[h]) )
         {
            SCIP_CALL( SCIPconshdlrUnlockVars(set->conshdlrs[h], set) );
         }
      }
   }

   /* free constraint array */
   BMSfreeMemoryArrayNull(&(*prob)->conss);
   
   /* free user problem data */
   if( (*prob)->transformed )
   {
      if( (*prob)->probdeltrans != NULL )
      {
         SCIP_CALL( (*prob)->probdeltrans(set->scip, &(*prob)->probdata) );
      }
   }
   else
   {
      if( (*prob)->probdelorig != NULL )
      {
         SCIP_CALL( (*prob)->probdelorig(set->scip, &(*prob)->probdata) );
      }
   }

   /* release problem variables */
   for( v = 0; v < (*prob)->nvars; ++v )
   {
      assert(SCIPvarGetProbindex((*prob)->vars[v]) >= 0);
      SCIP_CALL( SCIPvarRemove((*prob)->vars[v], blkmem, set, TRUE) );
      SCIP_CALL( SCIPvarRelease(&(*prob)->vars[v], blkmem, set, eventqueue, lp) );
   }
   BMSfreeMemoryArrayNull(&(*prob)->vars);

   /* release fixed problem variables */
   for( v = 0; v < (*prob)->nfixedvars; ++v )
   {
      assert(SCIPvarGetProbindex((*prob)->fixedvars[v]) == -1);
      SCIP_CALL( SCIPvarRelease(&(*prob)->fixedvars[v], blkmem, set, eventqueue, lp) );
   }
   BMSfreeMemoryArrayNull(&(*prob)->fixedvars);

   /* free deleted problem variables array */
   BMSfreeMemoryArrayNull(&(*prob)->deletedvars);

   /* free hash tables for names */
   if( (*prob)->varnames != NULL )
   {
      SCIPhashtableFree(&(*prob)->varnames);
   }
   if( (*prob)->consnames != NULL )
   {
      SCIPhashtableFree(&(*prob)->consnames);
   }
   BMSfreeMemoryArray(&(*prob)->name);
   BMSfreeMemory(prob);
   
   return SCIP_OKAY;
}

/** transform problem data into normalized form */
SCIP_RETCODE SCIPprobTransform(
   SCIP_PROB*            source,             /**< problem to transform */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_PROB**           target              /**< pointer to target problem data structure */
   )
{
   SCIP_VAR* targetvar;
   SCIP_CONS* targetcons;
   char transname[SCIP_MAXSTRLEN];
   int v;
   int c;
   int h;

   assert(set != NULL);
   assert(source != NULL);
   assert(blkmem != NULL);
   assert(target != NULL);

   SCIPdebugMessage("transform problem: original has %d variables\n", source->nvars);

   /* create target problem data (probdelorig and probtrans are not needed, probdata is set later) */
   (void) SCIPsnprintf(transname, SCIP_MAXSTRLEN, "t_%s", source->name);
   SCIP_CALL( SCIPprobCreate(target, blkmem, set, transname, source->probdelorig, source->probtrans, source->probdeltrans, 
         source->probinitsol, source->probexitsol, source->probcopy, NULL, TRUE) );
   SCIPprobSetObjsense(*target, source->objsense);

   /* transform objective limit */
   if( source->objlim < SCIP_INVALID )
      SCIPprobSetObjlim(*target, source->objlim);

   /* transform and copy all variables to target problem */
   SCIP_CALL( probEnsureVarsMem(*target, set, source->nvars) );
   for( v = 0; v < source->nvars; ++v )
   {
      SCIP_CALL( SCIPvarTransform(source->vars[v], blkmem, set, stat, source->objsense, &targetvar) );
      SCIP_CALL( SCIPprobAddVar(*target, blkmem, set, lp, branchcand, eventfilter, eventqueue, targetvar) );
      SCIP_CALL( SCIPvarRelease(&targetvar, blkmem, set, eventqueue, NULL) );
   }
   assert((*target)->nvars == source->nvars);

   /* call user data transformation */
   if( source->probtrans != NULL )
   {
      SCIP_CALL( source->probtrans(set->scip, source->probdata, &(*target)->probdata) );
   }
   else
      (*target)->probdata = source->probdata;

   /* transform and copy all constraints to target problem */
   for( c = 0; c < source->nconss; ++c )
   {
      SCIP_CALL( SCIPconsTransform(source->conss[c], blkmem, set, &targetcons) );
      SCIP_CALL( SCIPprobAddCons(*target, set, stat, targetcons) );
      SCIP_CALL( SCIPconsRelease(&targetcons, blkmem, set) );
   }

   /* lock variables for all constraint handlers that don't need constraints */
   for( h = 0; h < set->nconshdlrs; ++h )
   {
      if( !SCIPconshdlrNeedsCons(set->conshdlrs[h]) )
      {
         SCIP_CALL( SCIPconshdlrLockVars(set->conshdlrs[h], set) );
      }
   }
   
   /* objective value is always integral, iff original objective value is always integral and shift is integral */
   (*target)->objisintegral = source->objisintegral && SCIPsetIsIntegral(set, (*target)->objoffset);

   /* check, whether objective value is always integral by inspecting the problem, if it is the case adjust the
    * cutoff bound if primal solution is already known 
    */
   SCIP_CALL( SCIPprobCheckObjIntegral(*target, blkmem, set, stat, primal, tree, lp, eventqueue) );

   return SCIP_OKAY;
}

/** resets the global and local bounds of original variables in original problem to their original values */
SCIP_RETCODE SCIPprobResetBounds(
   SCIP_PROB*            prob,               /**< original problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< problem statistics */
   )
{
   int v;

   assert(prob != NULL);
   assert(!prob->transformed);
   assert(prob->nfixedvars == 0);

   for( v = 0; v < prob->nvars; ++v )
   {
      SCIP_CALL( SCIPvarResetBounds(prob->vars[v], blkmem, set, stat) );
   }

   return SCIP_OKAY;
}




/*
 * problem modification
 */

/** sets user problem data */
void SCIPprobSetData(
   SCIP_PROB*            prob,               /**< problem */
   SCIP_PROBDATA*        probdata            /**< user problem data to use */
   )
{
   assert(prob != NULL);

   prob->probdata = probdata;
}

/** inserts variable at the correct position in vars array, depending on its type */
static
void probInsertVar(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_VAR*             var                 /**< variable to insert */
   )
{
   int insertpos;
   int intstart;
   int implstart;
   int contstart;

   assert(prob != NULL);
   assert(prob->vars != NULL);
   assert(prob->nvars < prob->varssize);
   assert(var != NULL);
   assert(SCIPvarGetProbindex(var) == -1);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   /* original variables cannot go into transformed problem and transformed variables cannot go into original problem */
   assert((SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL) == prob->transformed);

   /* insert variable in array */
   insertpos = prob->nvars;
   intstart = prob->nbinvars;
   implstart = intstart + prob->nintvars;
   contstart = implstart + prob->nimplvars;

   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      prob->ncontvars++;
   else
   {
      if( insertpos > contstart )
      {
         prob->vars[insertpos] = prob->vars[contstart];
         SCIPvarSetProbindex(prob->vars[insertpos], insertpos);
         insertpos = contstart;
      }
      assert(insertpos == contstart);

      if( SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
         prob->nimplvars++;
      else
      {
         if( insertpos > implstart )
         {
            prob->vars[insertpos] = prob->vars[implstart];
            SCIPvarSetProbindex(prob->vars[insertpos], insertpos);
            insertpos = implstart;
         }
         assert(insertpos == implstart);

         if( SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER )
            prob->nintvars++;
         else
         {
            assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
            if( insertpos > intstart )
            {
               prob->vars[insertpos] = prob->vars[intstart];
               SCIPvarSetProbindex(prob->vars[insertpos], insertpos);
               insertpos = intstart;
            }
            assert(insertpos == intstart);

            prob->nbinvars++;
         }
      }
   }
   prob->nvars++;

   assert(prob->nvars == prob->nbinvars + prob->nintvars + prob->nimplvars + prob->ncontvars);
   assert((SCIPvarGetType(var) == SCIP_VARTYPE_BINARY && insertpos == prob->nbinvars - 1)
      || (SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER && insertpos == prob->nbinvars + prob->nintvars - 1)
      || (SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT && insertpos == prob->nbinvars + prob->nintvars + prob->nimplvars - 1)
      || (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS
         && insertpos == prob->nbinvars + prob->nintvars + prob->nimplvars + prob->ncontvars - 1));

   prob->vars[insertpos] = var;
   SCIPvarSetProbindex(var, insertpos);

   /* update number of column variables in problem */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      prob->ncolvars++;
   assert(0 <= prob->ncolvars && prob->ncolvars <= prob->nvars);
}

/** removes variable from vars array */
static
SCIP_RETCODE probRemoveVar(
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var                 /**< variable to remove */
   )
{
   int freepos;
   int intstart;
   int implstart;
   int contstart;

   assert(prob != NULL);
   assert(var != NULL);
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(prob->vars != NULL);
   assert(prob->vars[SCIPvarGetProbindex(var)] == var);

   intstart = prob->nbinvars;
   implstart = intstart + prob->nintvars;
   contstart = implstart + prob->nimplvars;

   switch( SCIPvarGetType(var) )
   {
   case SCIP_VARTYPE_BINARY:
      assert(0 <= SCIPvarGetProbindex(var) && SCIPvarGetProbindex(var) < intstart);
      prob->nbinvars--;
      break;
   case SCIP_VARTYPE_INTEGER:
      assert(intstart <= SCIPvarGetProbindex(var) && SCIPvarGetProbindex(var) < implstart);
      prob->nintvars--;
      break;
   case SCIP_VARTYPE_IMPLINT:
      assert(implstart <= SCIPvarGetProbindex(var) && SCIPvarGetProbindex(var) < contstart);
      prob->nimplvars--;
      break;
   case SCIP_VARTYPE_CONTINUOUS:
      assert(contstart <= SCIPvarGetProbindex(var) && SCIPvarGetProbindex(var) < prob->nvars);
      prob->ncontvars--;
      break;
   default:
      SCIPerrorMessage("unknown variable type\n");
      SCIPABORT();
   }

   /* move last binary, last integer, last implicit, and last continuous variable forward to fill the free slot */
   freepos = SCIPvarGetProbindex(var);
   if( freepos < intstart-1 )
   {
      /* move last binary variable to free slot */
      prob->vars[freepos] = prob->vars[intstart-1];
      SCIPvarSetProbindex(prob->vars[freepos], freepos);
      freepos = intstart-1;
   }
   if( freepos < implstart-1 )
   {
      /* move last integer variable to free slot */
      prob->vars[freepos] = prob->vars[implstart-1];
      SCIPvarSetProbindex(prob->vars[freepos], freepos);
      freepos = implstart-1;
   }
   if( freepos < contstart-1 )
   {
      /* move last implicit integer variable to free slot */
      prob->vars[freepos] = prob->vars[contstart-1];
      SCIPvarSetProbindex(prob->vars[freepos], freepos);
      freepos = contstart-1;
   }
   if( freepos < prob->nvars-1 )
   {
      /* move last implicit integer variable to free slot */
      prob->vars[freepos] = prob->vars[prob->nvars-1];
      SCIPvarSetProbindex(prob->vars[freepos], freepos);
      freepos = prob->nvars-1;
   }
   assert(freepos == prob->nvars-1);

   prob->nvars--;
   assert(prob->nvars == prob->nbinvars + prob->nintvars + prob->nimplvars + prob->ncontvars);

   /* update number of column variables in problem */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      prob->ncolvars--;
   assert(0 <= prob->ncolvars && prob->ncolvars <= prob->nvars);

   /* inform the variable that it is no longer in the problem; if necessary, delete it from the implication graph */
   SCIP_CALL( SCIPvarRemove(var, blkmem, set, FALSE) );

   return SCIP_OKAY;
}

/** adds variable to the problem and captures it */
SCIP_RETCODE SCIPprobAddVar(
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var                 /**< variable to add */
   )
{
   assert(prob != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(SCIPvarGetProbindex(var) == -1);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   /* original variables cannot go into transformed problem and transformed variables cannot go into original problem */
   assert((SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL) == prob->transformed);

#ifndef NDEBUG
   /* check if we add this variables to the same scip, where we created it */
   if( var->scip != set->scip )
   {
      SCIPerrorMessage("variable belongs to a different scip instance\n");
      return SCIP_INVALIDDATA;
   }
#endif

   /* capture variable */
   SCIPvarCapture(var);

   /* allocate additional memory */
   SCIP_CALL( probEnsureVarsMem(prob, set, prob->nvars+1) );

   /* insert variable in vars array and mark it to be in problem */
   probInsertVar(prob, var);

   /* add variable's name to the namespace */
   if( varHasName(var) && prob->varnames != NULL )
   {
      SCIP_CALL( SCIPhashtableInsert(prob->varnames, (void*)var) );
   }

   /* update branching candidates and pseudo and loose objective value in the LP */
   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL )
   {
      SCIP_CALL( SCIPbranchcandUpdateVar(branchcand, set, var) );
      SCIP_CALL( SCIPlpUpdateAddVar(lp, set, var) );
   }

   SCIPdebugMessage("added variable <%s> to problem (%d variables: %d binary, %d integer, %d implicit, %d continuous)\n",
      SCIPvarGetName(var), prob->nvars, prob->nbinvars, prob->nintvars, prob->nimplvars, prob->ncontvars);

   if( prob->transformed )
   {
      SCIP_EVENT* event;

      /* issue VARADDED event */
      SCIP_CALL( SCIPeventCreateVarAdded(&event, blkmem, var) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, eventfilter, &event) );
   }

   return SCIP_OKAY;
}

/** marks variable to be removed from the problem; however, the variable is NOT removed from the constraints */
SCIP_RETCODE SCIPprobDelVar(
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool*            deleted             /**< pointer to store whether marking variable to be deleted was successful */
   )
{
   assert(prob != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(eventqueue != NULL);
   assert(deleted != NULL);
   assert(SCIPvarGetProbindex(var) != -1);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   *deleted = FALSE;

   /* don't remove variables that are not in the problem */
   /**@todo what about negated variables? should the negation variable be removed instead? */
   if( SCIPvarGetProbindex(var) == -1 )
      return SCIP_OKAY;

   /* don't remove the direct counterpart of an original variable from the transformed problem, because otherwise
    * operations on the original variables would be applied to a NULL pointer
    */
   if( SCIPvarIsTransformedOrigvar(var) )
      return SCIP_OKAY;

   assert(SCIPvarGetNegatedVar(var) == NULL);

   SCIPdebugMessage("deleting variable <%s> from problem (%d variables: %d binary, %d integer, %d implicit, %d continuous)\n",
      SCIPvarGetName(var), prob->nvars, prob->nbinvars, prob->nintvars, prob->nimplvars, prob->ncontvars);

   /* mark variable to be deleted from the problem */
   SCIPvarMarkDeleted(var);

   if( prob->transformed )
   {
      SCIP_EVENT* event;

      /* issue VARDELETED event */
      SCIP_CALL( SCIPeventCreateVarDeleted(&event, blkmem, var) );
      SCIP_CALL( SCIPeventqueueAdd(eventqueue, blkmem, set, NULL, NULL, NULL, NULL, &event) );
   }

   /* remember that the variable should be deleted from the problem in SCIPprobPerformVarDeletions() */
   SCIP_CALL( probEnsureDeletedvarsMem(prob, set, prob->ndeletedvars+1) );
   prob->deletedvars[prob->ndeletedvars] = var;
   prob->ndeletedvars++;

   *deleted = TRUE;

   return SCIP_OKAY;
}

/** actually removes the deleted variables from the problem and releases them */
SCIP_RETCODE SCIPprobPerformVarDeletions(
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data (may be NULL) */
   SCIP_BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   int i;

   assert(prob != NULL);
   assert(set != NULL);

   /* delete variables from the constraints;
    * do this only in solving stage, in presolving, it is already handled by the constraint handlers
    */
   if( SCIPsetGetStage(set) == SCIP_STAGE_SOLVING )
   {
      for( i = 0; i < set->nconshdlrs; ++i )
      {
         SCIP_CALL( SCIPconshdlrDelVars(set->conshdlrs[i], blkmem, set, stat) );
      }
   }

   for( i = 0; i < prob->ndeletedvars; ++i )
   {
      SCIP_VAR* var;

      var = prob->deletedvars[i];

      /* don't delete the variable, if it was fixed or aggregated in the meantime */
      if( SCIPvarGetProbindex(var) >= 0 )
      {
         SCIPdebugMessage("perform deletion of <%s> [%p]\n", SCIPvarGetName(var), (void*)var);

         /* convert column variable back into loose variable, free LP column */
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_CALL( SCIPvarLoose(var, blkmem, set, eventqueue, prob, lp) );
         }
         
         /* update branching candidates and pseudo and loose objective value in the LP */
         if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL )
         {
            SCIP_CALL( SCIPlpUpdateDelVar(lp, set, var) );
            SCIP_CALL( SCIPbranchcandRemoveVar(branchcand, var) );
         }
         
         /* remove variable's name from the namespace */
         if( varHasName(var) && prob->varnames != NULL )
         {
            assert(SCIPhashtableExists(prob->varnames, (void*)var));
            SCIP_CALL( SCIPhashtableRemove(prob->varnames, (void*)var) );
         }

         /* remove variable from vars array and mark it to be not in problem */
         SCIP_CALL( probRemoveVar(prob, blkmem, set, var) );

         /* release variable */
         SCIP_CALL( SCIPvarRelease(&prob->deletedvars[i], blkmem, set, eventqueue, lp) );
      }
   }
   prob->ndeletedvars = 0;

   return SCIP_OKAY;
}

/** changes the type of a variable in the problem */
SCIP_RETCODE SCIPprobChgVarType(
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_VAR*             var,                /**< variable to add */
   SCIP_VARTYPE          vartype             /**< new type of variable */
   )
{
   assert(prob != NULL);
   assert(var != NULL);
   assert(SCIPvarGetProbindex(var) >= 0);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
   assert(branchcand != NULL || SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL);

   if( SCIPvarGetType(var) == vartype )
      return SCIP_OKAY;

   /* temporarily remove variable from branching candidates */
   if( branchcand != NULL )
   {
      SCIP_CALL( SCIPbranchcandRemoveVar(branchcand, var) );
   }

   /* temporarily remove variable from problem */
   SCIP_CALL( probRemoveVar(prob, blkmem, set, var) );

   /* change the type of the variable */
   SCIP_CALL( SCIPvarChgType(var, vartype) );

   /* reinsert variable into problem */
   probInsertVar(prob, var);

   /* update branching candidates */
   if( branchcand != NULL )
   {
      SCIP_CALL( SCIPbranchcandUpdateVar(branchcand, set, var) );
   }

   return SCIP_OKAY;
}

/** informs problem, that the given loose problem variable changed its status */
SCIP_RETCODE SCIPprobVarChangedStatus(
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert(prob != NULL);
   assert(var != NULL);
   assert(SCIPvarGetProbindex(var) != -1);

   /* get current status of variable */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      SCIPerrorMessage("variables cannot switch to ORIGINAL status\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_LOOSE:
      /* variable switched from column to loose */
      prob->ncolvars--;
      break;

   case SCIP_VARSTATUS_COLUMN:
      /* variable switched from non-column to column */
      prob->ncolvars++;
      break;

   case SCIP_VARSTATUS_FIXED:
   case SCIP_VARSTATUS_AGGREGATED:
   case SCIP_VARSTATUS_MULTAGGR:
   case SCIP_VARSTATUS_NEGATED:
      /* variable switched from unfixed to fixed (if it was fixed before, probindex would have been -1) */

      /* remove variable from problem */
      SCIP_CALL( probRemoveVar(prob, blkmem, set, var) );
      
      /* insert variable in fixedvars array */
      SCIP_CALL( probEnsureFixedvarsMem(prob, set, prob->nfixedvars+1) );
      prob->fixedvars[prob->nfixedvars] = var;
      prob->nfixedvars++;

      /* update branching candidates */
      SCIP_CALL( SCIPbranchcandUpdateVar(branchcand, set, var) );
      break;

   default:
      SCIPerrorMessage("invalid variable status <%d>\n", SCIPvarGetStatus(var));
      return SCIP_INVALIDDATA;
   }
   assert(0 <= prob->ncolvars && prob->ncolvars <= prob->nvars);

   return SCIP_OKAY;
}

/** adds constraint to the problem and captures it;
 *  a local constraint is automatically upgraded into a global constraint
 */
SCIP_RETCODE SCIPprobAddCons(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   assert(prob != NULL);
   assert(cons != NULL);
   assert(cons->addconssetchg == NULL);
   assert(cons->addarraypos == -1);

#ifndef NDEBUG
   /* check if we add this constraint to the same scip, where we create the constraint */
   if( cons->scip != set->scip )
   {
      SCIPerrorMessage("constraint belongs to different scip instance\n");
      return SCIP_INVALIDDATA;
   }
#endif

   SCIPdebugMessage("adding constraint <%s> to global problem -> %d constraints\n",
      SCIPconsGetName(cons), prob->nconss+1);

   /* mark the constraint as problem constraint, and remember the constraint's position */
   cons->addconssetchg = NULL;
   cons->addarraypos = prob->nconss;

   /* add the constraint to the problem's constraint array */
   SCIP_CALL( probEnsureConssMem(prob, set, prob->nconss+1) );
   prob->conss[prob->nconss] = cons;
   prob->nconss++;
   prob->maxnconss = MAX(prob->maxnconss, prob->nconss);

   /* undelete constraint, if it was globally deleted in the past */
   cons->deleted = FALSE;

   /* mark constraint to be globally valid */
   cons->local = FALSE;

   /* capture constraint */
   SCIPconsCapture(cons);

   /* add constraint's name to the namespace */
   if( consHasName(cons) && prob->consnames != NULL )
   {
      SCIP_CALL( SCIPhashtableInsert(prob->consnames, (void*)cons) );
   }

   /* if the problem is the transformed problem, activate and lock constraint */
   if( prob->transformed )
   {
      /* activate constraint */
      if( !SCIPconsIsActive(cons) )
      {
         SCIP_CALL( SCIPconsActivate(cons, set, stat, -1, (stat->nnodes <= 1)) );
      }

      /* if constraint is a check-constraint, lock roundings of constraint's variables */
      if( SCIPconsIsChecked(cons) )
      {
         SCIP_CALL( SCIPconsAddLocks(cons, set, +1, 0) );
      }
   }

   return SCIP_OKAY;
}

/** releases and removes constraint from the problem; if the user has not captured the constraint for his own use, the
 *  constraint may be invalid after the call
 */
SCIP_RETCODE SCIPprobDelCons(
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_CONS*            cons                /**< constraint to remove */
   )
{
   int arraypos;

   assert(prob != NULL);
   assert(blkmem != NULL);
   assert(cons != NULL);
   assert(cons->addconssetchg == NULL);
   assert(0 <= cons->addarraypos && cons->addarraypos < prob->nconss);
   assert(prob->conss != NULL);
   assert(prob->conss[cons->addarraypos] == cons);

   /* if the problem is the transformed problem, deactivate and unlock constraint */
   if( prob->transformed )
   {
      /* if constraint is a check-constraint, unlock roundings of constraint's variables */
      if( SCIPconsIsChecked(cons) )
      {
         SCIP_CALL( SCIPconsAddLocks(cons, set, -1, 0) );
      }

      /* deactivate constraint, if it is currently active */
      if( cons->active && !cons->updatedeactivate )
      {
         SCIP_CALL( SCIPconsDeactivate(cons, set, stat) );
      }
   }
   assert(!cons->active || cons->updatedeactivate);
   assert(!cons->enabled || cons->updatedeactivate);

   /* remove constraint's name from the namespace */
   if( consHasName(cons) && prob->consnames != NULL )
   {
      assert(SCIPhashtableExists(prob->consnames, (void*)cons));
      SCIP_CALL( SCIPhashtableRemove(prob->consnames, (void*)cons) );
   }

   /* remove the constraint from the problem's constraint array */
   arraypos = cons->addarraypos;
   prob->conss[arraypos] = prob->conss[prob->nconss-1];
   assert(prob->conss[arraypos] != NULL);
   assert(prob->conss[arraypos]->addconssetchg == NULL);
   prob->conss[arraypos]->addarraypos = arraypos;
   prob->nconss--;

   /* mark the constraint to be no longer in the problem */
   cons->addarraypos = -1;

   /* release constraint */
   SCIP_CALL( SCIPconsRelease(&cons, blkmem, set) );

   return SCIP_OKAY;
}

/** remembers the current number of constraints in the problem's internal data structure
 *  - resets maximum number of constraints to current number of constraints
 *  - remembers current number of constraints as starting number of constraints
 */
void SCIPprobMarkNConss(
   SCIP_PROB*            prob                /**< problem data */
   )
{
   assert(prob != NULL);

   /* remember number of constraints for statistic */
   prob->maxnconss = prob->nconss;
   prob->startnvars = prob->nvars;
   prob->startnconss = prob->nconss;
}

/** sets objective sense: minimization or maximization */
void SCIPprobSetObjsense(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_OBJSENSE         objsense            /**< new objective sense */
   )
{
   assert(prob != NULL);
   assert(prob->objsense == SCIP_OBJSENSE_MAXIMIZE || prob->objsense == SCIP_OBJSENSE_MINIMIZE);
   assert(objsense == SCIP_OBJSENSE_MAXIMIZE || objsense == SCIP_OBJSENSE_MINIMIZE);

   prob->objsense = objsense;
}

/** adds value to objective offset */
void SCIPprobAddObjoffset(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             addval              /**< value to add to objective offset */
   )
{
   assert(prob != NULL);
   assert(prob->transformed);

   SCIPdebugMessage("adding %g to objective offset %g: new offset = %g\n", addval, prob->objoffset, prob->objoffset + addval);
   prob->objoffset += addval;
}

/** sets limit on objective function, such that only solutions better than this limit are accepted */
void SCIPprobSetObjlim(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             objlim              /**< external objective limit */
   )
{
   assert(prob != NULL);

   prob->objlim = objlim;
}

/** informs the problem, that its objective value is always integral in every feasible solution */
void SCIPprobSetObjIntegral(
   SCIP_PROB*            prob                /**< problem data */
   )
{
   assert(prob != NULL);
   
   prob->objisintegral = TRUE;
}

/** sets integral objective value flag, if all variables with non-zero objective values are integral and have 
 *  integral objective value and also updates the cutoff bound if primal solution is already known
 */
SCIP_RETCODE SCIPprobCheckObjIntegral(
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   SCIP_Real obj;
   int v;

   assert(prob != NULL);
   
   /* if we know already, that the objective value is integral, nothing has to be done */
   if( prob->objisintegral )
      return SCIP_OKAY;
   
   /* if there exist unknown variables, we cannot conclude that the objective value is always integral */
   if( set->nactivepricers != 0 )
      return SCIP_OKAY;

   /* if the objective value offset is fractional, the value itself is possibly fractional */
   if( !SCIPsetIsIntegral(set, prob->objoffset) )
      return SCIP_OKAY;

   /* scan through the variables */
   for( v = 0; v < prob->nvars; ++v )
   {
      /* get objective value of variable */
      obj = SCIPvarGetObj(prob->vars[v]);

      /* check, if objective value is non-zero */
      if( !SCIPsetIsZero(set, obj) )
      {
         /* if variable's objective value is fractional, the problem's objective value may also be fractional */
         if( !SCIPsetIsIntegral(set, obj) )
            break;
         
         /* if variable with non-zero objective value is continuous, the problem's objective value may be fractional */
         if( SCIPvarGetType(prob->vars[v]) == SCIP_VARTYPE_CONTINUOUS )
            break;
      }
   }

   /* objective value is integral, if the variable loop scanned all variables */
   if( v == prob->nvars )
   {
      prob->objisintegral = TRUE;

      /* update upper bound and cutoff bound in primal data structure due to new internality information */
      SCIP_CALL( SCIPprimalUpdateObjoffset(primal, blkmem, set, stat, eventqueue, prob, tree, lp) );
   }

   return SCIP_OKAY;
}

/** if possible, scales objective function such that it is integral with gcd = 1 */
SCIP_RETCODE SCIPprobScaleObj(
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PRIMAL*          primal,             /**< primal data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue          /**< event queue */
   )
{
   int v;
   int nints;

   assert(prob != NULL);
   assert(set != NULL);

   /* do not change objective if there are pricers involved */
   if( set->nactivepricers != 0 )
      return SCIP_OKAY;

   nints = prob->nvars - prob->ncontvars;

   /* scan through the continuous variables */
   for( v = nints; v < prob->nvars; ++v )
   {
      SCIP_Real obj;

      /* get objective value of variable; it it is non-zero, no scaling can be applied */
      obj = SCIPvarGetObj(prob->vars[v]);
      if( !SCIPsetIsZero(set, obj) )
         break;
   }

   /* only continue if all continuous variables have obj = 0 */
   if( v == prob->nvars )
   {
      SCIP_Real* objvals;
      SCIP_Real intscalar;
      SCIP_Bool success;

      /* get temporary memory */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &objvals, nints) );

      /* get objective values of integer variables */
      for( v = 0; v < nints; ++v )
         objvals[v] = SCIPvarGetObj(prob->vars[v]);

      /* calculate integral scalar */
      SCIP_CALL( SCIPcalcIntegralScalar(objvals, nints, -SCIPsetEpsilon(set), +SCIPsetEpsilon(set), OBJSCALE_MAXDNOM, OBJSCALE_MAXSCALE,
         &intscalar, &success) );

      SCIPdebugMessage("integral objective scalar: success=%u, intscalar=%g\n", success, intscalar);

      if( success )
      {
         SCIP_Longint gcd;

         assert(intscalar > 0.0);

         /* calculate gcd of resulting integral coefficients */
         gcd = 0;
         for( v = 0; v < nints && gcd != 1; ++v )
         {
            SCIP_Longint absobj;

            absobj = (SCIP_Longint)(REALABS(objvals[v]) * intscalar + 0.5);
            if( gcd == 0 )
               gcd = absobj;
            else if( absobj > 0 )
               gcd = SCIPcalcGreComDiv(gcd, absobj);
         }
         if( gcd != 0 )
            intscalar /= gcd;
         SCIPdebugMessage("integral objective scalar: gcd=%"SCIP_LONGINT_FORMAT", intscalar=%g\n", gcd, intscalar);

         /* only apply scaling if the final scalar is small enough */
         if( intscalar <= OBJSCALE_MAXFINALSCALE )
         {
            /* apply scaling */
            if( !SCIPsetIsEQ(set, intscalar, 1.0) )
            {
               /* calculate scaled objective values */
               for( v = 0; v < nints; ++v )
               {
                  SCIP_Real newobj;
                  
                  /* check if new obj is really integral */
                  newobj = intscalar * SCIPvarGetObj(prob->vars[v]);
                  if( !SCIPsetIsFeasIntegral(set, newobj) )
                     break;
                  objvals[v] = SCIPsetFeasFloor(set, newobj);
               }

               /* change the variables' objective values and adjust objscale and objoffset */
               if( v == nints )
               {
                  for( v = 0; v < nints; ++v )
                  {
                     SCIPdebugMessage(" -> var <%s>: newobj = %.6f\n", SCIPvarGetName(prob->vars[v]), objvals[v]);
                     SCIP_CALL( SCIPvarChgObj(prob->vars[v], blkmem, set, primal, lp, eventqueue, objvals[v]) );
                  }
                  prob->objoffset *= intscalar;
                  prob->objscale /= intscalar;
                  prob->objisintegral = TRUE;
                  SCIPdebugMessage("integral objective scalar: objscale=%g\n", prob->objscale);

                  /* update upperbound and cutoffbound in primal data structure */
                  SCIP_CALL( SCIPprimalUpdateObjoffset(primal, blkmem, set, stat, eventqueue, prob, tree, lp) );
               }
            }
         }
      }

      /* free temporary memory */
      SCIPsetFreeBufferArray(set, &objvals);
   }

   return SCIP_OKAY;
}

/** remembers the current solution as root solution in the problem variables */
void SCIPprobStoreRootSol(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             roothaslp           /**< is the root solution from LP? */
   )
{
   int v;

   assert(prob != NULL);
   assert(prob->transformed);

   if( roothaslp )
   {
      for( v = 0; v < prob->nvars; ++v )
         SCIPvarStoreRootSol(prob->vars[v], stat, lp, roothaslp);

      SCIPlpSetRootLPIsRelax(lp, SCIPlpIsRelax(lp));
      SCIPlpStoreRootObjval(lp, set);
   }
}

/** informs problem, that the presolving process was finished, and updates all internal data structures */
SCIP_RETCODE SCIPprobExitPresolve(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** initializes problem for branch and bound process and resets all constraint's ages and histories of current run */
SCIP_RETCODE SCIPprobInitSolve(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   int c;
   int v;

   assert(prob != NULL);
   assert(prob->transformed);
   assert(set != NULL);

   /* reset constraint's ages */
   for( c = 0; c < prob->nconss; ++c )
   {
      SCIP_CALL( SCIPconsResetAge(prob->conss[c], set) );
   }

   /* initialize variables for solving */
   for( v = 0; v < prob->nvars; ++v )
      SCIPvarInitSolve(prob->vars[v]);

   /* call user data function */
   if( prob->probinitsol != NULL )
   {
      SCIP_CALL( prob->probinitsol(set->scip, prob->probdata) );
   }

   return SCIP_OKAY;
}

/** deinitializes problem after branch and bound process, and converts all COLUMN variables back into LOOSE variables */
SCIP_RETCODE SCIPprobExitSolve(
   SCIP_PROB*            prob,               /**< problem data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_Bool             restart             /**< was this exit solve call triggered by a restart? */
   )
{
   SCIP_VAR* var;
   int v;

   assert(prob != NULL);
   assert(prob->transformed);
   assert(set != NULL);

   /* call user data function */
   if( prob->probexitsol != NULL )
   {
      SCIP_CALL( prob->probexitsol(set->scip, prob->probdata, restart) );
   }

   /* convert all COLUMN variables back into LOOSE variables */
   if( prob->ncolvars > 0 )
   {
      for( v = 0; v < prob->nvars; ++v )
      {
         var = prob->vars[v];
         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_CALL( SCIPvarLoose(var, blkmem, set, eventqueue, prob, lp) );
         }
      }
   }
   assert(prob->ncolvars == 0);

   return SCIP_OKAY;
}




/*
 * problem information
 */

/** gets problem name */
const char* SCIPprobGetName(
   SCIP_PROB*            prob                /**< problem data */
   )
{
   assert(prob != NULL);
   return prob->name;
}

/** sets problem name */
SCIP_RETCODE SCIPprobSetName(
   SCIP_PROB*            prob,               /**< problem data */
   const char*           name                /**< name to be set */
   )
{
   assert(prob != NULL);

   BMSfreeMemoryArray(&(prob->name));
   SCIP_ALLOC( BMSduplicateMemoryArray(&(prob->name), name, strlen(name)+1) );

   return SCIP_OKAY;
}

/** gets user problem data */
SCIP_PROBDATA* SCIPprobGetData(
   SCIP_PROB*            prob                /**< problem */
   )
{
   assert(prob != NULL);

   return prob->probdata;
}

/** returns the external value of the given internal objective value */
SCIP_Real SCIPprobExternObjval(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             objval              /**< internal objective value */
   )
{
   assert(prob != NULL);
   assert(prob->transformed);
   assert(prob->objscale > 0.0);

   if( SCIPsetIsInfinity(set, objval) )
      return (SCIP_Real)prob->objsense * SCIPsetInfinity(set);
   else if( SCIPsetIsInfinity(set, -objval) )
      return -(SCIP_Real)prob->objsense * SCIPsetInfinity(set);
   else
      return (SCIP_Real)prob->objsense * prob->objscale * (objval + prob->objoffset);
}

/** returns the internal value of the given external objective value */
SCIP_Real SCIPprobInternObjval(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Real             objval              /**< external objective value */
   )
{
   assert(prob != NULL);
   assert(prob->transformed);
   assert(prob->objscale > 0.0);

   if( SCIPsetIsInfinity(set, objval) )
      return (SCIP_Real)prob->objsense * SCIPsetInfinity(set);
   else if( SCIPsetIsInfinity(set, -objval) )
      return -(SCIP_Real)prob->objsense * SCIPsetInfinity(set);
   else
      return (SCIP_Real)prob->objsense * objval/prob->objscale - prob->objoffset;
}

/** gets limit on objective function in external space */
SCIP_Real SCIPprobGetObjlim(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(prob != NULL);
   assert(set != NULL);

   return prob->objlim >= SCIP_INVALID ? (SCIP_Real)(prob->objsense) * SCIPsetInfinity(set) : prob->objlim;
}

/** returns whether the objective value is known to be integral in every feasible solution */
SCIP_Bool SCIPprobIsObjIntegral(
   SCIP_PROB*            prob                /**< problem data */
   )
{
   assert(prob != NULL);

   return prob->objisintegral;
}

/** returns variable of the problem with given name */
SCIP_VAR* SCIPprobFindVar(
   SCIP_PROB*            prob,               /**< problem data */
   const char*           name                /**< name of variable to find */
   )
{
   assert(prob != NULL);
   assert(name != NULL);

   if( prob->varnames == NULL )
   {
      SCIPerrorMessage("Cannot find variable if variable-names hashtable was disabled (due to parameter <misc/usevartable>)\n");
      SCIPABORT();/*lint --e{527}*/ /* only in debug mode */
      return NULL;
   }

   return (SCIP_VAR*)(SCIPhashtableRetrieve(prob->varnames, (char*)name));
}

/** returns constraint of the problem with given name */
SCIP_CONS* SCIPprobFindCons(
   SCIP_PROB*            prob,               /**< problem data */
   const char*           name                /**< name of variable to find */
   )
{
   assert(prob != NULL);
   assert(name != NULL);

   if( prob->consnames == NULL )
   {
      SCIPerrorMessage("Cannot find constraint if constraint-names hashtable was disabled (due to parameter <misc/useconstable>)\n");
      SCIPABORT();/*lint --e{527}*/ /* only in debug mode */
      return NULL;
   }

   return (SCIP_CONS*)(SCIPhashtableRetrieve(prob->consnames, (char*)name));
}

/** returns TRUE iff all columns, i.e. every variable with non-empty column w.r.t. all ever created rows, are present
 *  in the LP, and FALSE, if there are additional already existing columns, that may be added to the LP in pricing
 */
SCIP_Bool SCIPprobAllColsInLP(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_LP*              lp                  /**< current LP data */
   )
{
   assert(SCIPlpGetNCols(lp) <= prob->ncolvars && prob->ncolvars <= prob->nvars);

   return (SCIPlpGetNCols(lp) == prob->ncolvars && set->nactivepricers == 0);
}

/** displays current pseudo solution */
void SCIPprobPrintPseudoSol(
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_VAR* var;
   SCIP_Real solval;
   int v;
   
   for( v = 0; v < prob->nvars; ++v )
   {
      var = prob->vars[v];
      assert(var != NULL);
      solval = SCIPvarGetPseudoSol(var);
      if( !SCIPsetIsZero(set, solval) )
         SCIPmessagePrintInfo(" <%s>=%.15g", SCIPvarGetName(var), solval);
   }
   SCIPmessagePrintInfo("\n");
}

/** outputs problem statistics */
void SCIPprobPrintStatistics(
   SCIP_PROB*            prob,               /**< problem data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(prob != NULL);

   SCIPmessageFPrintInfo(file, "  Problem name     : %s\n", prob->name);
   SCIPmessageFPrintInfo(file, "  Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      prob->nvars, prob->nbinvars, prob->nintvars, prob->nimplvars, prob->ncontvars);
   SCIPmessageFPrintInfo(file, "  Constraints      : %d initial, %d maximal\n", prob->startnconss, prob->maxnconss);
   if( ! prob->transformed )
      SCIPmessageFPrintInfo(file, "  Objective sense  : %s\n", prob->objsense == SCIP_OBJSENSE_MINIMIZE ? "minimize" : "maximize");
}
