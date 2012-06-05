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

/**@file   cons_xor.c
 * @brief  Constraint handler for "xor" constraints,  \f$rhs = x_1 \oplus x_2 \oplus \dots  \oplus x_n\f$
 * @author Tobias Achterberg
 * @author Stefan Heinz
 * @author Michael Winkler
 *
 * This constraint handler deals with "xor" constraint. These are constraint of the form:
 *
 * \f[
 *    rhs = x_1 \oplus x_2 \oplus \dots  \oplus x_n
 * \f]
 *
 * where \f$x_i\f$ is a binary variable for all \f$i\f$ and \f$rhs\f$ is bool. The variables \f$x\f$'s are called
 * operators. This constraint is satisfied if \f$rhs\f$ is TRUE and an odd number of the operators are TRUE or if the
 * \f$rhs\f$ is FALSE and a even number of operators are TRUE. Hence, if the sum of \f$rhs\f$ and operators is even.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/pub_misc.h"
#include "scip/cons_xor.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "xor"
#define CONSHDLR_DESC          "constraint handler for xor constraints: r = xor(x1, ..., xn)"
#define CONSHDLR_SEPAPRIORITY   +850200 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -850200 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -850200 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP

#define EVENTHDLR_NAME         "xor"
#define EVENTHDLR_DESC         "event handler for xor constraints"

#define DEFAULT_PRESOLPAIRWISE     TRUE /**< should pairwise constraint comparison be performed in presolving? */
#define HASHSIZE_XORCONS         131101 /**< minimal size of hash table in logicor constraint tables */
#define DEFAULT_PRESOLUSEHASHING   TRUE /**< should hash table be used for detecting redundant constraints in advance */
#define NMINCOMPARISONS          200000 /**< number for minimal pairwise presolving comparisons */
#define MINGAINPERNMINCOMPARISONS 1e-06 /**< minimal gain per minimal pairwise presolving comparisons to repeat pairwise comparison round */

#define NROWS 4



/*
 * Data structures
 */

/** constraint data for xor constraints */
struct SCIP_ConsData
{
   SCIP_VAR**            vars;               /**< variables (including resultant) in the xor operation */
   SCIP_VAR*             intvar;             /**< internal variable for LP relaxation */
   SCIP_ROW*             rows[NROWS];        /**< rows for linear relaxation of xor constraint */
   int                   nvars;              /**< number of variables (including resultant) in xor operation */
   int                   varssize;           /**< size of vars array */
   int                   watchedvar1;        /**< position of first watched operator variable */
   int                   watchedvar2;        /**< position of second watched operator variable */
   int                   filterpos1;         /**< event filter position of first watched operator variable */
   int                   filterpos2;         /**< event filter position of second watched operator variable */
   SCIP_Bool             rhs;                /**< right hand side of the constraint */
   unsigned int          propagated:1;       /**< is constraint already preprocessed/propagated? */
   unsigned int          sorted:1;           /**< are the constraint's variables sorted? */
   unsigned int          changed:1;          /**< was constraint changed since last pair preprocessing round? */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for events on watched variables */
   SCIP_Bool             presolpairwise;     /**< should pairwise constraint comparison be performed in presolving? */
   SCIP_Bool             presolusehashing;   /**< should hash table be used for detecting redundant constraints in advance */
};




/*
 * Propagation rules
 */

enum Proprule
{
   PROPRULE_1,                          /**< all except one variable fixed  =>  fix remaining variable */
   PROPRULE_INVALID                     /**< propagation was applied without a specific propagation rule */
};
typedef enum Proprule PROPRULE;




/*
 * Local methods
 */

/** installs rounding locks for the given variable in the given xor constraint */
static
SCIP_RETCODE lockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding in both directions may violate the constraint */
   SCIP_CALL( SCIPlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** removes rounding locks for the given variable in the given xor constraint */
static
SCIP_RETCODE unlockRounding(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
   SCIP_VAR*             var                 /**< variable of constraint entry */
   )
{
   /* rounding in both directions may violate the constraint */
   SCIP_CALL( SCIPunlockVarCons(scip, var, cons, TRUE, TRUE) );

   return SCIP_OKAY;
}

/** creates constraint handler data */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );

   /* get event handler for catching events on watched variables */
   (*conshdlrdata)->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( (*conshdlrdata)->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for xor constraints not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   
   return SCIP_OKAY;
}

/** frees constraint handler data */
static
SCIP_RETCODE conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   SCIPfreeMemory(scip, conshdlrdata);

   return SCIP_OKAY;
}

/** stores the given variable numbers as watched variables, and updates the event processing */
static
SCIP_RETCODE consdataSwitchWatchedvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< xor constraint data */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   watchedvar1,        /**< new first watched variable */
   int                   watchedvar2         /**< new second watched variable */
   )
{
   assert(consdata != NULL);
   assert(watchedvar1 == -1 || watchedvar1 != watchedvar2);
   assert(watchedvar1 != -1 || watchedvar2 == -1);
   assert(watchedvar1 == -1 || (0 <= watchedvar1 && watchedvar1 < consdata->nvars));
   assert(watchedvar2 == -1 || (0 <= watchedvar2 && watchedvar2 < consdata->nvars));

   /* if one watched variable is equal to the old other watched variable, just switch positions */
   if( watchedvar1 == consdata->watchedvar2 || watchedvar2 == consdata->watchedvar1 )
   {
      int tmp;
      
      tmp = consdata->watchedvar1;
      consdata->watchedvar1 = consdata->watchedvar2;
      consdata->watchedvar2 = tmp;
      tmp = consdata->filterpos1;
      consdata->filterpos1 = consdata->filterpos2;
      consdata->filterpos2 = tmp;
   }
   assert(watchedvar1 == -1 || watchedvar1 != consdata->watchedvar2);
   assert(watchedvar2 == -1 || watchedvar2 != consdata->watchedvar1);

   /* drop events on old watched variables */
   if( consdata->watchedvar1 != -1 && consdata->watchedvar1 != watchedvar1 )
   {
      assert(consdata->filterpos1 != -1);
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar1], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (SCIP_EVENTDATA*)consdata, consdata->filterpos1) );
   }
   if( consdata->watchedvar2 != -1 && consdata->watchedvar2 != watchedvar2 )
   {
      assert(consdata->filterpos2 != -1);
      SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[consdata->watchedvar2], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (SCIP_EVENTDATA*)consdata, consdata->filterpos2) );
   }

   /* catch events on new watched variables */
   if( watchedvar1 != -1 && watchedvar1 != consdata->watchedvar1 )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[watchedvar1], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (SCIP_EVENTDATA*)consdata, &consdata->filterpos1) );
   }
   if( watchedvar2 != -1 && watchedvar2 != consdata->watchedvar2 )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[watchedvar2], SCIP_EVENTTYPE_BOUNDCHANGED, eventhdlr,
            (SCIP_EVENTDATA*)consdata, &consdata->filterpos2) );
   }

   /* set the new watched variables */
   consdata->watchedvar1 = watchedvar1;
   consdata->watchedvar2 = watchedvar2;

   return SCIP_OKAY;
}

/** ensures, that the vars array can store at least num entries */
static
SCIP_RETCODE consdataEnsureVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< linear constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(consdata != NULL);
   assert(consdata->nvars <= consdata->varssize);
   
   if( num > consdata->varssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, consdata->varssize, newsize) );
      consdata->varssize = newsize;
   }
   assert(num <= consdata->varssize);

   return SCIP_OKAY;
}

/** creates constraint data for xor constraint */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the constraint data */
   SCIP_Bool             rhs,                /**< right hand side of the constraint */
   int                   nvars,              /**< number of variables in the xor operation */
   SCIP_VAR**            vars                /**< variables in xor operation */
   )
{
   int r;

   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   /* store the resultant variable as first variable in the array */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );

   (*consdata)->rhs = rhs;
   (*consdata)->intvar = NULL;
   for( r = 0; r < NROWS; ++r )
      (*consdata)->rows[r] = NULL;
   (*consdata)->nvars = nvars;
   (*consdata)->varssize = nvars;
   (*consdata)->watchedvar1 = -1;
   (*consdata)->watchedvar2 = -1;
   (*consdata)->filterpos1 = -1;
   (*consdata)->filterpos2 = -1;
   (*consdata)->propagated = FALSE;
   (*consdata)->sorted = FALSE;
   (*consdata)->changed = TRUE;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
   }

   return SCIP_OKAY;
}

/** releases LP row of constraint data */
static
SCIP_RETCODE consdataFreeRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   int r;

   assert(consdata != NULL);

   for( r = 0; r < NROWS; ++r )
   {
      if( consdata->rows[r] != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->rows[r]) );
      }
   }

   return SCIP_OKAY;
}

/** frees constraint data for xor constraint */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to the constraint data */
   SCIP_EVENTHDLR*       eventhdlr           /**< event handler to call for the event processing */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   if( SCIPisTransformed(scip) )
   {
      /* drop events for watched variables */
      SCIP_CALL( consdataSwitchWatchedvars(scip, *consdata, eventhdlr, -1, -1) );
   }
   else
   {
      assert((*consdata)->watchedvar1 == -1);
      assert((*consdata)->watchedvar2 == -1);
   }

   /* release LP row */
   SCIP_CALL( consdataFreeRows(scip, *consdata) );

   /* release internal variable */
   if( (*consdata)->intvar != NULL )
   {
      /* if internal variable is not defined anymore by any checked constraint, delete it from the problem;
       * otherwise, it could be fixed to a wrong value in dual presolving
       */
      if( SCIPvarMayRoundDown((*consdata)->intvar) && SCIPvarMayRoundUp((*consdata)->intvar) )
      {
         SCIP_Bool deleted;

         SCIP_CALL( SCIPdelVar(scip, (*consdata)->intvar, &deleted) );
      }
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->intvar) );
   }

   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->varssize);
   SCIPfreeBlockMemory(scip, consdata);
 
   return SCIP_OKAY;
}

/** prints xor constraint to file stream */
static
SCIP_RETCODE consdataPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< xor constraint data */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             endline             /**< should an endline be set? */
   )
{
   assert(consdata != NULL);

   /* start variable list */
   SCIPinfoMessage(scip, file, "xor(");

   /* print variable list */
   SCIP_CALL( SCIPwriteVarsList(scip, file, consdata->vars, consdata->nvars, TRUE, ',') );

   /* close variable list and write right hand side */
   SCIPinfoMessage(scip, file, ") = %d", consdata->rhs);
   
   if( endline )
      SCIPinfoMessage(scip, file, "\n");

   return SCIP_OKAY;
}

/** adds coefficient to xor constraint */
static
SCIP_RETCODE addCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_VAR*             var                 /**< variable to add to the constraint */
   )
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(var != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows[0] == NULL);

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert(var != NULL);
   assert(transformed == SCIPvarIsTransformed(var));

   SCIP_CALL( consdataEnsureVarsSize(scip, consdata, consdata->nvars+1) );
   consdata->vars[consdata->nvars] = var;
   consdata->nvars++;
   consdata->sorted = (consdata->nvars == 1);
   consdata->changed = TRUE;

   /* install the rounding locks for the new variable */
   SCIP_CALL( lockRounding(scip, cons, var) );

   /**@todo update LP rows */
   if( consdata->rows[0] != NULL )
   {
      SCIPerrorMessage("cannot add coefficients to xor constraint after LP relaxation was created\n");
      return SCIP_INVALIDCALL;
   }

   return SCIP_OKAY;
}

/** deletes coefficient at given position from xor constraint data */
static
SCIP_RETCODE delCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;

   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nvars);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(consdata->vars[pos]));

   /* remove the rounding locks of the deleted variable */
   SCIP_CALL( unlockRounding(scip, cons, consdata->vars[pos]) );

   if( SCIPconsIsTransformed(cons) )
   {
      /* if the position is watched, stop watching the position */
      if( consdata->watchedvar1 == pos )
      {
         SCIP_CALL( consdataSwitchWatchedvars(scip, consdata, eventhdlr, consdata->watchedvar2, -1) );
      }
      if( consdata->watchedvar2 == pos )
      {
         SCIP_CALL( consdataSwitchWatchedvars(scip, consdata, eventhdlr, consdata->watchedvar1, -1) );
      }
   }
   assert(pos != consdata->watchedvar1);
   assert(pos != consdata->watchedvar2);

   /* move the last variable to the free slot */
   consdata->vars[pos] = consdata->vars[consdata->nvars-1];
   consdata->nvars--;

   /* if the last variable (that moved) was watched, update the watched position */
   if( consdata->watchedvar1 == consdata->nvars )
      consdata->watchedvar1 = pos;
   if( consdata->watchedvar2 == consdata->nvars )
      consdata->watchedvar2 = pos;

   consdata->propagated = FALSE;
   consdata->sorted = FALSE;
   consdata->changed = TRUE;

   return SCIP_OKAY;
}

/** index comparison method of and constraints: compares two indices of the variable set in the xor constraint */
static
SCIP_DECL_SORTINDCOMP(consdataCompVar)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata = (SCIP_CONSDATA*)dataptr;

   assert(consdata != NULL);
   assert(0 <= ind1 && ind1 < consdata->nvars);
   assert(0 <= ind2 && ind2 < consdata->nvars);

   return SCIPvarCompareActiveAndNegated(consdata->vars[ind1], consdata->vars[ind2]);
}

/** sorts and constraint's variables by non-decreasing variable index */
static
SCIP_RETCODE consdataSort(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->nvars <= 1 )
      consdata->sorted = TRUE;
   else if( !consdata->sorted )
   {
      SCIP_VAR* varv;

      int* perm;
      int v;
      int i;
      int nexti;

      /* get temporary memory to store the sorted permutation */
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, consdata->nvars) );

      /* calculate sorted permutation of coefficients */
      SCIPsort(perm, consdataCompVar, (void*)consdata, consdata->nvars);

      /* permute the variables in the constraint according to the resulting permutation */
      for( v = 0; v < consdata->nvars; ++v )
      {
         if( perm[v] != v )
         {
            SCIP_Bool iswatchedvar1;
            SCIP_Bool iswatchedvar2;

            varv = consdata->vars[v];
            iswatchedvar1 = (consdata->watchedvar1 == v);
            iswatchedvar2 = (consdata->watchedvar2 == v);
            i = v;
            do
            {
               assert(0 <= perm[i] && perm[i] < consdata->nvars);
               assert(perm[i] != i);
               consdata->vars[i] = consdata->vars[perm[i]];
               if( consdata->watchedvar1 == perm[i] )
                  consdata->watchedvar1 = i;
               if( consdata->watchedvar2 == perm[i] )
                  consdata->watchedvar2 = i;
               nexti = perm[i];
               perm[i] = i;
               i = nexti;
            }
            while( perm[i] != v );
            consdata->vars[i] = varv;
            if( iswatchedvar1 )
               consdata->watchedvar1 = i;
            if( iswatchedvar2 )
               consdata->watchedvar2 = i;
            perm[i] = i;
         }
      }
      consdata->sorted = TRUE;

#ifdef SCIP_DEBUG
      /* check sorting */
      for( v = 0; v < consdata->nvars; ++v )
      {
         assert(v == consdata->nvars-1 || SCIPvarCompareActiveAndNegated(consdata->vars[v], consdata->vars[v+1]) <= 0);
         assert(perm[v] == v);
      }
#endif

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &perm);
   }
   assert(consdata->sorted);

   return SCIP_OKAY;
}


/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyXorcons)
{  /*lint --e{715}*/
   /* the key is the element itself */ 
   return elem;
}

/** returns TRUE iff both keys are equal; two constraints are equal if they have the same variables */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqXorcons)
{
   SCIP* scip;
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;
   int i;

   consdata1 = SCIPconsGetData((SCIP_CONS*)key1);
   consdata2 = SCIPconsGetData((SCIP_CONS*)key2);
   scip = (SCIP*)userptr; 
   assert(scip != NULL);
   
   /* checks trivial case */
   if( consdata1->nvars != consdata2->nvars )
      return FALSE;

   /* sorts the constraints */
   SCIP_CALL_ABORT( consdataSort(scip, consdata1) );
   SCIP_CALL_ABORT( consdataSort(scip, consdata2) );

   for( i = 0; i < consdata1->nvars ; ++i )
   {
      /* tests if variables are equal */
      if( consdata1->vars[i] != consdata2->vars[i] )
      {
         assert(SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == 1 || 
            SCIPvarCompare(consdata1->vars[i], consdata2->vars[i]) == -1);
         return FALSE;
      }
      assert(SCIPvarCompareActiveAndNegated(consdata1->vars[i], consdata2->vars[i]) == 0);
   } 
   
   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValXorcons)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   unsigned int hashval;
   int minidx;
   int mididx;
   int maxidx;
   
   consdata = SCIPconsGetData((SCIP_CONS*)key);
   assert(consdata != NULL);
   assert(consdata->sorted);
   assert(consdata->nvars > 0);

   /* only active, fixed or negated variables are allowed */
   assert(consdata->vars[0] != NULL);
   assert(consdata->vars[consdata->nvars / 2] != NULL);
   assert(consdata->vars[consdata->nvars - 1] != NULL);
   assert(SCIPvarIsActive(consdata->vars[0]) || SCIPvarGetStatus(consdata->vars[0]) == SCIP_VARSTATUS_NEGATED || SCIPvarGetStatus(consdata->vars[0]) == SCIP_VARSTATUS_FIXED);
   assert(SCIPvarIsActive(consdata->vars[consdata->nvars / 2]) || SCIPvarGetStatus(consdata->vars[consdata->nvars / 2]) == SCIP_VARSTATUS_NEGATED || SCIPvarGetStatus(consdata->vars[consdata->nvars / 2]) == SCIP_VARSTATUS_FIXED);
   assert(SCIPvarIsActive(consdata->vars[consdata->nvars - 1]) || SCIPvarGetStatus(consdata->vars[consdata->nvars - 1]) == SCIP_VARSTATUS_NEGATED || SCIPvarGetStatus(consdata->vars[consdata->nvars - 1]) == SCIP_VARSTATUS_FIXED);

   minidx = SCIPvarGetIndex(consdata->vars[0]);
   mididx = SCIPvarGetIndex(consdata->vars[consdata->nvars / 2]);
   maxidx = SCIPvarGetIndex(consdata->vars[consdata->nvars - 1]);
   /* note that for all indices it does not hold that they are sorted, because variables are sorted with
    * SCIPvarCompareActiveAndNegated (see var.c)
    */

   hashval = (consdata->nvars << 29) + (minidx << 22) + (mididx << 11) + maxidx; /*lint !e701*/

   return hashval;
}

/** updates the flags of the first constraint according to the ones of the second constraint */
static
SCIP_RETCODE updateFlags(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons0,              /**< constraint that should stay */
   SCIP_CONS*            cons1               /**< constraint that should be deleted */
   )
{
   if( SCIPconsIsInitial(cons1) )
   {
      SCIP_CALL( SCIPsetConsInitial(scip, cons0, TRUE) );
   }
   if( SCIPconsIsSeparated(cons1) )
   {
      SCIP_CALL( SCIPsetConsSeparated(scip, cons0, TRUE) );
   }
   if( SCIPconsIsEnforced(cons1) )
   {
      SCIP_CALL( SCIPsetConsEnforced(scip, cons0, TRUE) );
   }
   if( SCIPconsIsChecked(cons1) )
   {
      SCIP_CALL( SCIPsetConsChecked(scip, cons0, TRUE) );
   }
   if( SCIPconsIsPropagated(cons1) )
   {
      SCIP_CALL( SCIPsetConsPropagated(scip, cons0, TRUE) );
   }
   if( !SCIPconsIsDynamic(cons1) )
   {
      SCIP_CALL( SCIPsetConsDynamic(scip, cons0, FALSE) );
   }
   if( !SCIPconsIsRemovable(cons1) )
   {
      SCIP_CALL( SCIPsetConsRemovable(scip, cons0, FALSE) );
   }
   if( SCIPconsIsStickingAtNode(cons1) )
   {
      SCIP_CALL( SCIPsetConsStickingAtNode(scip, cons0, TRUE) );
   }

   return SCIP_OKAY;
}

/** deletes all fixed variables and all pairs equal variables variables */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   int*                  nchgcoefs           /**< pointer to add up the number of changed coefficients */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->nvars == 0 || consdata->vars != NULL);
   assert(nchgcoefs != NULL);

   SCIPdebugMessage("before fixings: ");
   SCIPdebug( SCIP_CALL(consdataPrint(scip, consdata, NULL, TRUE)) );

   v = 0;
   while( v < consdata->nvars )
   {
      SCIP_VAR* var;

      var = consdata->vars[v];
      assert(SCIPvarIsBinary(var));

      if( SCIPvarGetUbGlobal(var) < 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0));
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
         (*nchgcoefs)++;
      }
      else if( SCIPvarGetLbGlobal(var) > 0.5 )
      {
         assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(var), 1.0));
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
         consdata->rhs = !consdata->rhs;
         (*nchgcoefs)++;
      }
      else
      {
         SCIP_VAR* repvar;
         SCIP_Bool negated;
         
         /* get binary representative of variable */
         SCIP_CALL( SCIPgetBinvarRepresentative(scip, var, &repvar, &negated) );

         /* check, if the variable should be replaced with the representative */
         if( repvar != var )
         {
            /* delete old (aggregated) variable */
            SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );

            /* add representative instead */
            SCIP_CALL( addCoef(scip, cons, eventhdlr, repvar) );
         }
         else
            ++v;
      }
   }

   /* sort the variables in the constraint */
   SCIP_CALL( consdataSort(scip, consdata) );

   SCIPdebugMessage("after sort    : ");
   SCIPdebug( SCIP_CALL(consdataPrint(scip, consdata, NULL, TRUE)) );

   /* delete pairs of equal or negated variables; scan from back to front because deletion doesn't affect the
    * order of the front variables
    */
   v = consdata->nvars-2;
   while ( v >= 0 )
   {
      if( consdata->vars[v] == consdata->vars[v+1] )
      {
         /* delete both variables */
         SCIPdebugMessage("xor constraint <%s>: deleting pair of equal variables <%s>\n",
            SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[v]));
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v+1) );
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
         (*nchgcoefs) += 2;
         v = MIN(v, consdata->nvars-1);
      }
      else if( consdata->vars[v] == SCIPvarGetNegatedVar(consdata->vars[v+1]) )
      {
         /* delete both variables and negate the rhs */
         SCIPdebugMessage("xor constraint <%s>: deleting pair of negated variables <%s> and <%s>\n",
            SCIPconsGetName(cons), SCIPvarGetName(consdata->vars[v]), SCIPvarGetName(consdata->vars[v+1]));
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v+1) );
         SCIP_CALL( delCoefPos(scip, cons, eventhdlr, v) );
         (*nchgcoefs) += 2;
         consdata->rhs = !consdata->rhs;
         v = MIN(v, consdata->nvars-1);
      }
      else
         assert(SCIPvarGetProbvar(consdata->vars[v]) != SCIPvarGetProbvar(consdata->vars[v+1]));
      --v;
   }

   SCIPdebugMessage("after fixings : ");
   SCIPdebug( SCIP_CALL(consdataPrint(scip, consdata, NULL, TRUE)) );

   return SCIP_OKAY;
}

/** creates LP row corresponding to xor constraint: 
 *    x1 + ... + xn - 2q == rhs
 *  with internal integer variable q;
 *  in the special case of 3 variables and c = 0, the following linear system is created:
 *    + x - y - z <= 0
 *    - x + y - z <= 0
 *    - x - y + z <= 0
 *    + x + y + z <= 2
 *  in the special case of 3 variables and c = 1, the following linear system is created:
 *    - x + y + z <= 1
 *    + x - y + z <= 1
 *    + x + y - z <= 1
 *    - x - y - z <= -1
 */
static 
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_CONSDATA* consdata;
   char varname[SCIP_MAXSTRLEN];

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->rows[0] == NULL);

   if( SCIPconsIsModifiable(cons) || consdata->nvars != 3 )  /* resultant is member of vars array */
   {
      SCIP_Real rhsval;

      /* create internal variable, if not yet existing */
      if( consdata->intvar == NULL )
      {
         int ub;

         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "%s_int", SCIPconsGetName(cons));
         ub = consdata->nvars/2;
         SCIP_CALL( SCIPcreateVar(scip, &consdata->intvar, varname, 0.0, (SCIP_Real)ub, 0.0,
               consdata->nvars >= 4 ? SCIP_VARTYPE_INTEGER : SCIP_VARTYPE_BINARY,
               SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, consdata->intvar) );

         /* install the rounding locks for the internal variable */
         SCIP_CALL( lockRounding(scip, cons, consdata->intvar) );
      }

      /* create LP row (resultant variable is also stored in vars array) */
      rhsval = (consdata->rhs ? 1.0 : 0.0);
      SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[0], SCIPconsGetName(cons), rhsval, rhsval,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
      SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[0], consdata->intvar, -2.0) );
      SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[0], consdata->nvars, consdata->vars, 1.0) );
   }
   else if( !consdata->rhs )
   {
      char rowname[SCIP_MAXSTRLEN];
      int r;

      /* create the <= 0 rows with one positive sign */
      for( r = 0; r < 3; ++r )
      {
         int v;

         (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_%d", SCIPconsGetName(cons), r);
         SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[r], rowname, -SCIPinfinity(scip), 0.0,
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
         for( v = 0; v < 3; ++v )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[r], consdata->vars[v], v == r ? +1.0 : -1.0) );
         }
      }

      /* create the <= 2 row with all positive signs */
      (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_3", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[3], rowname, -SCIPinfinity(scip), 2.0,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
      SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[3], consdata->nvars, consdata->vars, 1.0) );
   }
   else
   {
      char rowname[SCIP_MAXSTRLEN];
      int r;

      /* create the <= 1 rows with one negative sign */
      for( r = 0; r < 3; ++r )
      {
         int v;

         (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_%d", SCIPconsGetName(cons), r);
         SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[r], rowname, -SCIPinfinity(scip), 1.0,
               SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
         for( v = 0; v < 3; ++v )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, consdata->rows[r], consdata->vars[v], v == r ? -1.0 : +1.0) );
         }
      }

      /* create the <= -1 row with all negative signs */
      (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "%s_3", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->rows[3], rowname, -SCIPinfinity(scip), -1.0,
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
      SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, consdata->rows[3], consdata->nvars, consdata->vars, -1.0) );
   }

   return SCIP_OKAY;
}  

/** adds linear relaxation of or constraint to the LP */
static 
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to check */
   )
{
   SCIP_CONSDATA* consdata;
   int r;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->rows[0] == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->rows[0] != NULL);
   for( r = 0; r < NROWS; ++r )
   {
      if( consdata->rows[r] != NULL && !SCIProwIsInLP(consdata->rows[r]) )
      {
         SCIP_CALL( SCIPaddCut(scip, NULL, consdata->rows[r], FALSE) );
      }
   }

   return SCIP_OKAY;
}

/** returns whether all rows of the LP relaxation are in the current LP */
static
SCIP_Bool allRowsInLP(
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->rows[0] == NULL )      /* LP relaxation does not exist */
      return FALSE;
   else
   {
      int r;
      for( r = 0; r < NROWS; ++r )
      {
         if( consdata->rows[r] != NULL && !SCIProwIsInLP(consdata->rows[r]) )
            return FALSE;
      }
      return TRUE;
   }
}

/** checks xor constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< solution to check, NULL for current solution */
   SCIP_Bool             checklprows,        /**< should LP rows be checked? */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;

   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violated = FALSE;

   /* check feasibility of constraint if necessary */
   if( checklprows || !allRowsInLP(consdata) )
   {
      SCIP_Real solval;
      int i;
      SCIP_Bool odd;

      /* increase age of constraint; age is reset to zero, if a violation was found */
      SCIP_CALL( SCIPincConsAge(scip, cons) );
      
      /* check, if all variables (including resultant) and the rhs sum up to an even value */
      odd = consdata->rhs;
      for( i = 0; i < consdata->nvars; ++i )
      {
         solval = SCIPgetSolVal(scip, sol, consdata->vars[i]);
         assert(SCIPisFeasIntegral(scip, solval));
         odd = (odd != (solval > 0.5));
      }
      if( odd )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
         *violated = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** separates current LP solution */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool*            separated           /**< pointer to store whether a cut was found */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real feasibility;
   int r;

   assert(separated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *separated = FALSE;

   /* create row for the linear relaxation */
   if( consdata->rows[0] == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->rows[0] != NULL);

   /* test rows for feasibility and add it, if it is infeasible */
   for( r = 0; r < NROWS; ++r )
   {
      if( sol != NULL || (consdata->rows[r] != NULL && !SCIProwIsInLP(consdata->rows[r])) )
      {
         feasibility = SCIPgetRowSolFeasibility(scip, consdata->rows[r], sol);
         if( SCIPisFeasNegative(scip, feasibility) )
         {
            SCIP_CALL( SCIPaddCut(scip, sol, consdata->rows[r], FALSE) );
            *separated = TRUE;
         }
      }
   }

   return SCIP_OKAY;
}

/** analyzes conflicting assignment on given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< xor constraint that detected the conflict */
   )
{
   SCIP_CONSDATA* consdata;
   int v;

   /* conflict analysis can only be applied in solving stage */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip) );
   for( v = 0; v < consdata->nvars; ++v )
   {
      SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->vars[v]) );
   }

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** propagates constraint with the following rule:
 *   (1) all except one variable fixed  =>  fix remaining variable
 */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< xor constraint to be processed */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   int watchedvar1;
   int watchedvar2;
   int i;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   SCIP_Bool odd;

   assert(cutoff != NULL);
   assert(nfixedvars != NULL);

   /* propagation can only be applied, if we know all operator variables */
   if( SCIPconsIsModifiable(cons) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   vars = consdata->vars;
   nvars = consdata->nvars;

   /* don't process the constraint, if the watched variables weren't fixed to any value since last propagation call */
   if( consdata->propagated )
      return SCIP_OKAY;

   /* increase age of constraint; age is reset to zero, if a conflict or a propagation was found */
   if( !SCIPinRepropagation(scip) )
   {
      SCIP_CALL( SCIPincConsAge(scip, cons) );
   }

   /* propagation cannot be applied, if we have at least two unfixed variables left;
    * that means, we only have to watch (i.e. capture events) of two variables, and switch to other variables
    * if these ones get fixed
    */
   watchedvar1 = consdata->watchedvar1;
   watchedvar2 = consdata->watchedvar2;

   /* check, if watched variables are still unfixed */
   if( watchedvar1 != -1 )
   {
      if( SCIPvarGetLbLocal(vars[watchedvar1]) > 0.5 || SCIPvarGetUbLocal(vars[watchedvar1]) < 0.5 )
         watchedvar1 = -1;
   }
   if( watchedvar2 != -1 )
   {
      if( SCIPvarGetLbLocal(vars[watchedvar2]) > 0.5 || SCIPvarGetUbLocal(vars[watchedvar2]) < 0.5 )
         watchedvar2 = -1;
   }
   
   /* if only one watched variable is still unfixed, make it the first one */
   if( watchedvar1 == -1 )
   {
      watchedvar1 = watchedvar2;
      watchedvar2 = -1;
   }
   assert(watchedvar1 != -1 || watchedvar2 == -1);

   /* if the watched variables are invalid (fixed), find new ones if existing; count the parity */
   odd = consdata->rhs;
   if( watchedvar2 == -1 )
   {
      for( i = 0; i < nvars; ++i )
      {
         if( SCIPvarGetLbLocal(vars[i]) > 0.5 )
            odd = !odd;
         else if( SCIPvarGetUbLocal(vars[i]) > 0.5 )
         {
            if( watchedvar1 == -1 )
            {
               assert(watchedvar2 == -1);
               watchedvar1 = i;
            }
            else if( watchedvar1 != i )
            {
               watchedvar2 = i;
               break;
            }
         }
      }
   }
   assert(watchedvar1 != -1 || watchedvar2 == -1);

   /* if all variables are fixed, we can decide the feasibility of the constraint */
   if( watchedvar1 == -1 )
   {
      assert(watchedvar2 == -1);
      
      SCIPdebugMessage("constraint <%s>: all vars fixed -> constraint is %s\n",
         SCIPconsGetName(cons), odd ? "infeasible" : "feasible");
      if( odd )
      {
         /* use conflict analysis to get a conflict constraint out of the conflicting assignment */
         SCIP_CALL( analyzeConflict(scip, cons) );

         *cutoff = TRUE;
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );

      return SCIP_OKAY;
   }

   /* if only one variable is not fixed, this variable can be deduced */
   if( watchedvar2 == -1 )
   {
      assert(watchedvar1 != -1);
      
      SCIPdebugMessage("constraint <%s>: only one unfixed variable -> fix <%s> to %u\n",
         SCIPconsGetName(cons), SCIPvarGetName(vars[watchedvar1]), odd);
      SCIP_CALL( SCIPinferBinvarCons(scip, vars[watchedvar1], odd, cons, (int)PROPRULE_1, &infeasible, &tightened) );
      assert(!infeasible);
      assert(tightened);
      (*nfixedvars)++;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );

      return SCIP_OKAY;
   }

   /* switch to the new watched variables */
   SCIP_CALL( consdataSwitchWatchedvars(scip, consdata, eventhdlr, watchedvar1, watchedvar2) );

   /* mark the constraint propagated */
   consdata->propagated = TRUE;

   return SCIP_OKAY;
}

/** resolves a conflict on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rule (see propagateCons()):
 *   (1) all except one variable fixed  =>  fix remaining variable
 */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced */
   PROPRULE              proprule,           /**< propagation rule that deduced the value */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_RESULT*          result              /**< pointer to store the result of the propagation conflict resolving call */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   vars = consdata->vars;
   nvars = consdata->nvars;

   switch( proprule )
   {
   case PROPRULE_1:
      /* the variable was infered, because all other variables were fixed */
      for( i = 0; i < nvars; ++i )
      {
         if( vars[i] != infervar )
         {
            assert(SCIPisEQ(scip, SCIPvarGetLbAtIndex(vars[i], bdchgidx, FALSE), 
                  SCIPvarGetUbAtIndex(vars[i], bdchgidx, FALSE)));
            SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i]) );
         }
         else
            assert(SCIPisEQ(scip, SCIPvarGetLbAtIndex(vars[i], bdchgidx, TRUE), 
                  SCIPvarGetUbAtIndex(vars[i], bdchgidx, TRUE)));
      }
      *result = SCIP_SUCCESS;
      break;

   case PROPRULE_INVALID:
   default:
      SCIPerrorMessage("invalid inference information %d in xor constraint <%s>\n", proprule, SCIPconsGetName(cons));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** compares each constraint with all other constraints for possible redundancy and removes or changes constraint 
 *  accordingly; in contrast to preprocessConstraintPairs(), it uses a hash table 
 */
static
SCIP_RETCODE detectRedundantConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   int*                  firstchange,        /**< pointer to store first changed constraint */
   int*                  nchgcoefs,          /**< pointer to add up the number of changed coefficients */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   SCIP_Bool*            cutoff              /**< pointer to store TRUE, if a cutoff was found */
)
{
   SCIP_HASHTABLE* hashtable;
   int hashtablesize;
   int c;

   assert(conss != NULL);
   assert(ndelconss != NULL);

   /* create a hash table for the constraint set */
   hashtablesize = SCIPcalcHashtableSize(10*nconss);
   hashtablesize = MAX(hashtablesize, HASHSIZE_XORCONS);
   SCIP_CALL( SCIPhashtableCreate(&hashtable, blkmem, hashtablesize,
         hashGetKeyXorcons, hashKeyEqXorcons, hashKeyValXorcons, (void*) scip) );

   /* check all constraints in the given set for redundancy */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONS* cons0;
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata0;
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      cons0 = conss[c];

      if( !SCIPconsIsActive(cons0) || SCIPconsIsModifiable(cons0) )
         continue;

      /* get constraint handler data */
      conshdlr = SCIPconsGetHdlr(cons0);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      /* it can happen that during preprocessing some variables got aggregated and a constraint now has not active
       * variables inside so we need to remove them for sorting
       */
      /* remove all variables that are fixed to zero and all pairs of variables fixed to one;
       * merge multiple entries of the same or negated variables
       */
      SCIP_CALL( applyFixings(scip, cons0, conshdlrdata->eventhdlr, nchgcoefs) );

      consdata0 = SCIPconsGetData(cons0);
      /* sort the constraint */
      SCIP_CALL( consdataSort(scip, consdata0) );

      /* get constraint from current hash table with same variables as cons0 */
      cons1 = (SCIP_CONS*)(SCIPhashtableRetrieve(hashtable, (void*)cons0));
 
      if( cons1 != NULL )
      {
         SCIP_CONSDATA* consdata1;

         assert(SCIPconsIsActive(cons1));
         assert(!SCIPconsIsModifiable(cons1));
      
         consdata1 = SCIPconsGetData(cons1);
         
         assert(consdata0 != NULL && consdata1 != NULL);
         assert(consdata0->nvars >= 1 && consdata0->nvars == consdata1->nvars);
         
         assert(consdata0->sorted && consdata1->sorted);
         assert(consdata0->vars[0] == consdata1->vars[0]);

         if( consdata0->rhs != consdata1->rhs )
         {
            *cutoff = TRUE;
            goto TERMINATE;
         }

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( updateFlags(scip, cons1, cons0) ); 

         /* delete consdel */
         SCIP_CALL( SCIPdelCons(scip, cons0) );
         (*ndelconss)++;

         /* update the first changed constraint to begin the next aggregation round with */
         if( consdata0->changed && SCIPconsGetPos(cons1) < *firstchange )
            *firstchange = SCIPconsGetPos(cons1);

         assert(SCIPconsIsActive(cons1));
      }
      else
      {
         /* no such constraint in current hash table: insert cons0 into hash table */  
         SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) cons0) );
      }
   }

 TERMINATE:
   /* free hash table */
   SCIPhashtableFree(&hashtable);

   return SCIP_OKAY;
}

/** compares constraint with all prior constraints for possible redundancy or aggregation,
 *  and removes or changes constraint accordingly
 */
static
SCIP_RETCODE preprocessConstraintPairs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   firstchange,        /**< first constraint that changed since last pair preprocessing round */
   int                   chkind,             /**< index of constraint to check against all prior indices upto startind */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nfixedvars,         /**< pointer to add up the number of found domain reductions */
   int*                  naggrvars,          /**< pointer to count number of aggregated variables */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  nchgcoefs           /**< pointer to add up the number of changed coefficients */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons0;
   SCIP_CONSDATA* consdata0;
   SCIP_Bool cons0changed;
   int c;

   assert(conss != NULL);
   assert(firstchange <= chkind);
   assert(cutoff != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(ndelconss != NULL);
   assert(nchgcoefs != NULL);

   /* get the constraint to be checked against all prior constraints */
   cons0 = conss[chkind];
   assert(SCIPconsIsActive(cons0));
   assert(!SCIPconsIsModifiable(cons0));

   consdata0 = SCIPconsGetData(cons0);
   assert(consdata0 != NULL);
   assert(consdata0->nvars >= 1);

   /* get constraint handler data */
   conshdlr = SCIPconsGetHdlr(cons0);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* it can happen that during preprocessing some variables got aggregated and a constraint now has not active
    * variables inside so we need to remove them for sorting
    */
   /* remove all variables that are fixed to zero and all pairs of variables fixed to one;
    * merge multiple entries of the same or negated variables
    */
   SCIP_CALL( applyFixings(scip, cons0, conshdlrdata->eventhdlr, nchgcoefs) );

   /* sort cons0 */
   SCIP_CALL( consdataSort(scip, consdata0) );

   /* check constraint against all prior constraints */
   cons0changed = consdata0->changed;
   consdata0->changed = FALSE;
   for( c = (cons0changed ? 0 : firstchange); c < chkind && !(*cutoff) && SCIPconsIsActive(cons0) && !SCIPisStopped(scip);
        ++c )
   {
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata1;
      SCIP_VAR* singlevar0;
      SCIP_VAR* singlevar1;
      SCIP_Bool parity;
      SCIP_Bool cons0hastwoothervars;
      SCIP_Bool cons1hastwoothervars;
      SCIP_Bool aborted;
      int v0;
      int v1;

      cons1 = conss[c];

      /* ignore inactive and modifiable constraints */
      if( !SCIPconsIsActive(cons1) || SCIPconsIsModifiable(cons1) )
         continue;

      consdata1 = SCIPconsGetData(cons1);
      assert(consdata1 != NULL);

      /* it can happen that during preprocessing some variables got aggregated and a constraint now has not active
       * variables inside so we need to remove them for sorting
       */
      /* remove all variables that are fixed to zero and all pairs of variables fixed to one;
       * merge multiple entries of the same or negated variables
       */
      SCIP_CALL( applyFixings(scip, cons1, conshdlrdata->eventhdlr, nchgcoefs) );

      SCIPdebugMessage("preprocess xor constraint pair <%s>[chg:%u] and <%s>[chg:%u]\n",
         SCIPconsGetName(cons0), cons0changed, SCIPconsGetName(cons1), consdata1->changed);

      /* if both constraints were not changed since last round, we can ignore the pair */
      if( !cons0changed && !consdata1->changed )
         continue;

      assert(consdata1->nvars >= 1);

      /* sort cons1 */
      SCIP_CALL( consdataSort(scip, consdata1) );

      /* check whether
       *  (a) one problem variable set is a subset of the other, or
       *  (b) the problem variable sets are almost equal with only one variable in each constraint that is not
       *      member of the other
       */
      assert(consdata0->sorted);
      assert(consdata1->sorted);
      aborted = FALSE;
      parity = (consdata0->rhs ^ consdata1->rhs);
      cons0hastwoothervars = FALSE;
      cons1hastwoothervars = FALSE;
      singlevar0 = NULL;
      singlevar1 = NULL;
      v0 = 0;
      v1 = 0;
      while( (v0 < consdata0->nvars || v1 < consdata1->nvars) && !aborted )
      {
         int cmp;

         assert(v0 <= consdata0->nvars);
         assert(v1 <= consdata1->nvars);

         if( v0 == consdata0->nvars )
            cmp = +1;
         else if( v1 == consdata1->nvars )
            cmp = -1;
         else
            cmp = SCIPvarCompareActiveAndNegated(consdata0->vars[v0], consdata1->vars[v1]);

         switch( cmp )
         {
         case -1:
            /* variable doesn't appear in cons1 */
            assert(v0 < consdata0->nvars);
            if( singlevar0 == NULL )
            {
               singlevar0 = consdata0->vars[v0];
               if( cons1hastwoothervars )
                  aborted = TRUE;
            }
            else
            {
               cons0hastwoothervars = TRUE;
               if( singlevar1 != NULL )
                  aborted = TRUE;
            }
            v0++;
            break;

         case +1:
            /* variable doesn't appear in cons0 */
            assert(v1 < consdata1->nvars);
            if( singlevar1 == NULL )
            {
               singlevar1 = consdata1->vars[v1];
               if( cons0hastwoothervars )
                  aborted = TRUE;
            }
            else
            {
               cons1hastwoothervars = TRUE;
               if( singlevar0 != NULL )
                  aborted = TRUE;
            }
            v1++;
            break;

         case 0:
            /* variable appears in both constraints */
            assert(v0 < consdata0->nvars);
            assert(v1 < consdata1->nvars);
            assert(SCIPvarGetProbvar(consdata0->vars[v0]) == SCIPvarGetProbvar(consdata1->vars[v1]));
            if( consdata0->vars[v0] != consdata1->vars[v1] )
            {
               assert(SCIPvarGetNegatedVar(consdata0->vars[v0]) == consdata1->vars[v1]);
               parity = !parity;
            }
            v0++;
            v1++;
            break;

         default:
            SCIPerrorMessage("invalid comparison result\n");
            SCIPABORT();
         }
      }

      /* check if a useful presolving is possible */
      if( (cons0hastwoothervars && singlevar1 != NULL) || (cons1hastwoothervars && singlevar0 != NULL) )
         continue;

      /* check if one problem variable set is a subset of the other */
      if( singlevar0 == NULL && singlevar1 == NULL )
      {
         /* both constraints are equal */
         if( !parity )
         {
            /* even parity: constraints are redundant */
            SCIPdebugMessage("xor constraints <%s> and <%s> are redundant: delete <%s>\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), SCIPconsGetName(cons1));
            SCIPdebug(SCIPprintCons(scip, cons0, NULL));
            SCIPdebug(SCIPprintCons(scip, cons1, NULL));
            SCIP_CALL( SCIPdelCons(scip, cons1) );
            (*ndelconss)++;
         }
         else
         {
            /* odd parity: constraints are contradicting */
            SCIPdebugMessage("xor constraints <%s> and <%s> are contradicting\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            SCIPdebug(SCIPprintCons(scip, cons0, NULL));
            SCIPdebug(SCIPprintCons(scip, cons1, NULL));
            *cutoff = TRUE;
         }
      }
      else if( singlevar1 == NULL )
      {
         /* cons1 is a subset of cons0 */
         if( !cons0hastwoothervars )
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;

            /* only one additional variable in cons0: fix this variable according to the parity */
            SCIPdebugMessage("xor constraints <%s> and <%s> yield sum %u == <%s>\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), parity, SCIPvarGetName(singlevar0));
            SCIPdebug(SCIPprintCons(scip, cons0, NULL));
            SCIPdebug(SCIPprintCons(scip, cons1, NULL));
            SCIP_CALL( SCIPfixVar(scip, singlevar0, parity ? 1.0 : 0.0, &infeasible, &fixed) );
            assert(infeasible || fixed);
            *cutoff = *cutoff || infeasible;
            (*nfixedvars)++;
            SCIP_CALL( SCIPdelCons(scip, cons1) );
            (*ndelconss)++;
         }
         else
         {
            int v;

            /* more than one additional variable in cons0: add cons1 to cons0, thus eliminating the equal variables */
            SCIPdebugMessage("xor constraint <%s> is superset of <%s> with parity %u\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), parity);
            SCIPdebug(SCIPprintCons(scip, cons0, NULL));
            SCIPdebug(SCIPprintCons(scip, cons1, NULL));
            for( v = 0; v < consdata1->nvars; ++v )
            {
               SCIP_CALL( addCoef(scip, cons0, conshdlrdata->eventhdlr, consdata1->vars[v]) );
            }
            SCIP_CALL( applyFixings(scip, cons0, conshdlrdata->eventhdlr, nchgcoefs) );
            assert(SCIPconsGetData(cons0) == consdata0);
            assert(consdata0->nvars >= 2); /* at least the two "other" variables should remain in the constraint */
            SCIP_CALL( consdataSort(scip, consdata0) );
         }
      }
      else if( singlevar0 == NULL )
      {
         /* cons0 is a subset of cons1 */
         if( !cons1hastwoothervars )
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;

            /* only one additional variable in cons1: fix this variable according to the parity */
            SCIPdebugMessage("xor constraints <%s> and <%s> yield sum %u == <%s>\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), parity, SCIPvarGetName(singlevar1));
            SCIPdebug(SCIPprintCons(scip, cons0, NULL));
            SCIPdebug(SCIPprintCons(scip, cons1, NULL));
            SCIP_CALL( SCIPfixVar(scip, singlevar1, parity ? 1.0 : 0.0, &infeasible, &fixed) );
            assert(infeasible || fixed);
            *cutoff = *cutoff || infeasible;
            (*nfixedvars)++;
            SCIP_CALL( SCIPdelCons(scip, cons1) );
            (*ndelconss)++;
         }
         else
         {
            int v;

            /* more than one additional variable in cons1: add cons0 to cons1, thus eliminating the equal variables */
            SCIPdebugMessage("xor constraint <%s> is subset of <%s> with parity %u\n",
               SCIPconsGetName(cons0), SCIPconsGetName(cons1), parity);
            SCIPdebug(SCIPprintCons(scip, cons0, NULL));
            SCIPdebug(SCIPprintCons(scip, cons1, NULL));
            for( v = 0; v < consdata0->nvars; ++v )
            {
               SCIP_CALL( addCoef(scip, cons1, conshdlrdata->eventhdlr, consdata0->vars[v]) );
            }
            SCIP_CALL( applyFixings(scip, cons1, conshdlrdata->eventhdlr, nchgcoefs) );
            assert(SCIPconsGetData(cons1) == consdata1);
            assert(consdata1->nvars >= 2); /* at least the two "other" variables should remain in the constraint */
            SCIP_CALL( consdataSort(scip, consdata1) );
         }
      }
      else
      {
         SCIP_Bool infeasible;
         SCIP_Bool redundant;
         SCIP_Bool aggregated;

         assert(!cons0hastwoothervars);
         assert(!cons1hastwoothervars);

         /* sum of constraints is parity == singlevar0 xor singlevar1: aggregate variables and delete cons1 */
         SCIPdebugMessage("xor constraints <%s> and <%s> yield sum %u == xor(<%s>,<%s>)\n",
            SCIPconsGetName(cons0), SCIPconsGetName(cons1), parity, SCIPvarGetName(singlevar0),
            SCIPvarGetName(singlevar1));
         if( !parity )
         {
            /* aggregate singlevar0 == singlevar1 */
            SCIP_CALL( SCIPaggregateVars(scip, singlevar1, singlevar0, 1.0, -1.0, 0.0,
                  &infeasible, &redundant, &aggregated) );
         }
         else
         {
            /* aggregate singlevar0 == 1-singlevar1 */
            SCIP_CALL( SCIPaggregateVars(scip, singlevar1, singlevar0, 1.0, 1.0, 1.0,
                  &infeasible, &redundant, &aggregated) );
         }
         assert(infeasible || redundant || SCIPdoNotAggr(scip));

         *cutoff = *cutoff || infeasible;
         if( aggregated )
            (*naggrvars)++;

         if( redundant )
         {
            SCIP_CALL( SCIPdelCons(scip, cons1) );
            (*ndelconss)++;
         }
#if 0
      /* if aggregation in the core of SCIP is not changed we do not need to call applyFixing, this would be the correct
       * way
       */
      /* remove all variables that are fixed to zero and all pairs of variables fixed to one;
       * merge multiple entries of the same or negated variables
       */
      SCIP_CALL( applyFixings(scip, cons0, conshdlrdata->eventhdlr, nchgcoefs) );
#endif

      }
   }

   return SCIP_OKAY;
}





/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyXor)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrXor(scip) );
 
   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( conshdlrdataFree(scip, &conshdlrdata) );

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitXor NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitXor NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreXor NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreXor NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolXor NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolXor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   /* release and free the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      SCIP_CALL( consdataFreeRows(scip, consdata) );
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( consdataFree(scip, consdata, conshdlrdata->eventhdlr) );

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransXor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);
   assert(sourcedata->nvars >= 1);
   assert(sourcedata->vars != NULL);

   /* create target constraint data */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->rhs, sourcedata->nvars, sourcedata->vars) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpXor)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      assert(SCIPconsIsInitial(conss[i]));
      SCIP_CALL( addRelaxation(scip, conss[i]) );
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpXor)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   int c;

   *result = SCIP_DIDNOTFIND;

   /* separate all useful constraints */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], NULL, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   } 

   /* combine constraints to get more cuts */
   /**@todo combine constraints to get further cuts */

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolXor)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   int c;

   *result = SCIP_DIDNOTFIND;

   /* separate all useful constraints */
   for( c = 0; c < nusefulconss; ++c )
   {
      SCIP_CALL( separateCons(scip, conss[c], sol, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   } 

   /* combine constraints to get more cuts */
   /**@todo combine constraints to get further cuts */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpXor)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, FALSE, &violated) );
      if( violated )
      {
         SCIP_Bool separated;

         SCIP_CALL( separateCons(scip, conss[i], NULL, &separated) );
         assert(separated); /* because the solution is integral, the separation always finds a cut */
         *result = SCIP_SEPARATED;
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsXor)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, TRUE, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckXor)
{  /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, checklprows, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;
         
         if( printreason )
         {
            int v;
            int sum = 0;
            SCIP_CONSDATA* consdata;

            consdata = SCIPconsGetData(conss[i]);
            assert( consdata != NULL );

            SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) );
            
            for( v = 0; v < consdata->nvars; ++v )
            {
               if( SCIPgetSolVal(scip, sol, consdata->vars[i]) > 0.5 )
                  sum++;
            }
            SCIPinfoMessage(scip, NULL, "violation: %d operands are set to TRUE\n", sum );
         }

         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool cutoff;
   int nfixedvars;
   int c;
   
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   nfixedvars = 0;

   /* propagate all useful constraints */
   for( c = 0; c < nusefulconss && !cutoff; ++c )
   {
      SCIP_CALL( propagateCons(scip, conss[c], conshdlrdata->eventhdlr, &cutoff, &nfixedvars) );
   }

   /* return the correct result */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolXor)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   SCIP_Bool cutoff;
   SCIP_Bool delay;
   SCIP_Bool redundant;
   SCIP_Bool aggregated;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldndelconss;
   int oldnchgcoefs;
   int firstchange;
   int c;

   assert(result != NULL);

   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;
   oldndelconss = *ndelconss;
   oldnchgcoefs = *nchgcoefs;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* process constraints */
   cutoff = FALSE;
   delay = FALSE;
   firstchange = INT_MAX;
   for( c = 0; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->propagated = FALSE;

      /* remember the first changed constraint to begin the next aggregation round with */
      if( firstchange == INT_MAX && consdata->changed )
         firstchange = c;

      /* remove all variables that are fixed to zero and all pairs of variables fixed to one;
       * merge multiple entries of the same or negated variables
       */
      SCIP_CALL( applyFixings(scip, cons, conshdlrdata->eventhdlr, nchgcoefs) );

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, cons, conshdlrdata->eventhdlr, &cutoff, nfixedvars) );

      if( !cutoff && !SCIPconsIsDeleted(cons) && !SCIPconsIsModifiable(cons) )
      {
         assert(consdata->nvars >= 2); /* otherwise, propagateCons() has deleted the constraint */

         /* if only two variables are left, both have to be equal or opposite, depending on the rhs */
         if( consdata->nvars == 2 )
         {
            SCIPdebugMessage("xor constraint <%s> has only two unfixed variables, rhs=%u\n",
               SCIPconsGetName(cons), consdata->rhs);
            
            assert(consdata->vars != NULL);
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[0]), 0.0));
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[0]), 1.0));
            assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(consdata->vars[1]), 0.0));
            assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(consdata->vars[1]), 1.0));
            
            if( !consdata->rhs )
            {
               /* aggregate variables: vars[0] - vars[1] == 0 */
               SCIPdebugMessage(" -> aggregate <%s> == <%s>\n", SCIPvarGetName(consdata->vars[0]),
                  SCIPvarGetName(consdata->vars[1]));
               SCIP_CALL( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[1], 1.0, -1.0, 0.0,
                     &cutoff, &redundant, &aggregated) );
            }
            else
            {
               /* aggregate variables: vars[0] + vars[1] == 1 */
               SCIPdebugMessage(" -> aggregate <%s> == 1 - <%s>\n", SCIPvarGetName(consdata->vars[0]),
                  SCIPvarGetName(consdata->vars[1]));
               SCIP_CALL( SCIPaggregateVars(scip, consdata->vars[0], consdata->vars[1], 1.0, 1.0, 1.0,
                     &cutoff, &redundant, &aggregated) );
            }
            assert(redundant || SCIPdoNotAggr(scip));

            if( aggregated )
            {
               assert(redundant);
               (*naggrvars)++;
            }

            if( redundant )
            {
               /* delete constraint */
               SCIP_CALL( SCIPdelCons(scip, cons) );
               (*ndelconss)++;
            }
         }
      }
   }

   /* process pairs of constraints: check them for equal operands in order to aggregate resultants;
    * only apply this expensive procedure, if the single constraint preprocessing did not find any reductions
    * (otherwise, we delay the presolving to be called again next time)
    */
   if( !cutoff )
   {
      if( *nfixedvars == oldnfixedvars && *naggrvars == oldnaggrvars )
      {
         if( firstchange < nconss && conshdlrdata->presolusehashing ) 
         {
            /* detect redundant constraints; fast version with hash table instead of pairwise comparison */
            SCIP_CALL( detectRedundantConstraints(scip, SCIPblkmem(scip), conss, nconss, &firstchange, nchgcoefs, ndelconss, &cutoff) );
         }
         if( conshdlrdata->presolpairwise )
         {
            SCIP_Longint npaircomparisons;
            int lastndelconss;
            npaircomparisons = 0;
            lastndelconss = *ndelconss;

            for( c = firstchange; c < nconss && !cutoff && !SCIPisStopped(scip); ++c )
            {
               if( SCIPconsIsActive(conss[c]) && !SCIPconsIsModifiable(conss[c]) )
               {
                  npaircomparisons += (SCIPconsGetData(conss[c])->changed) ? (SCIP_Longint) c : ((SCIP_Longint) c - (SCIP_Longint) firstchange);

                  SCIP_CALL( preprocessConstraintPairs(scip, conss, firstchange, c,
                        &cutoff, nfixedvars, naggrvars, ndelconss, nchgcoefs) );

                  if( npaircomparisons > NMINCOMPARISONS )
                  {
                     if( ((SCIP_Real) (*ndelconss - lastndelconss)) / ((SCIP_Real) npaircomparisons) < MINGAINPERNMINCOMPARISONS )
                        break;
                     lastndelconss = *ndelconss;
                     npaircomparisons = 0;
                  }
               }
            }
         }
      }
      else
         delay = TRUE;
   }

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( delay )
      *result = SCIP_DELAYED;
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars || *ndelconss > oldndelconss
      || *nchgcoefs > oldnchgcoefs )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropXor)
{  /*lint --e{715}*/
   SCIP_CALL( resolvePropagation(scip, cons, infervar, (PROPRULE)inferinfo, bdchgidx, result) );

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockXor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* external variables */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   /* internal variable */
   if( consdata->intvar != NULL )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->intvar, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveXor NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveXor NULL


/** constraint enabling notification method of constraint handler */
#define consEnableXor NULL


/** constraint disabling notification method of constraint handler */
#define consDisableXor NULL


/** variable deletion method of constraint handler */
#define consDelvarsXor NULL


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintXor)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
 
   SCIP_CALL( consdataPrint(scip, SCIPconsGetData(cons), file, FALSE) );
    
   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyXor)
{  /*lint --e{715}*/
   SCIP_VAR** sourcevars;
   SCIP_VAR** targetvars;
   const char* consname;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);

   (*valid) = TRUE;

   /* get variables and coefficients of the source constraint */
   sourcevars = SCIPgetVarsXor(sourcescip, sourcecons);
   nvars = SCIPgetNVarsXor(sourcescip, sourcecons);
   
   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   if( nvars == 0 )
   {
      SCIP_CALL( SCIPcreateConsXor(scip, cons, consname, SCIPgetRhsXor(sourcescip, sourcecons), 0, NULL,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
      return SCIP_OKAY;
   }

   /* duplicate variable array */
   SCIP_CALL( SCIPallocBufferArray(scip, &targetvars, nvars) );

   /* map variables of the source constraint to variables of the target SCIP */
   for( v = 0; v < nvars && *valid; ++v )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[v], &targetvars[v], varmap, consmap, global, valid) );
      assert(!(*valid) || targetvars[v] != NULL);
   }

   /* only create the target constraints, if all variables could be copied */
   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsXor(scip, cons, consname, SCIPgetRhsXor(sourcescip, sourcecons), nvars, targetvars,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   /* free buffer array */
   SCIPfreeBufferArray(scip, &targetvars);
   
   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseXor)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   char* endptr;
   int requiredsize;
   int varssize;
   int nvars;

   SCIPdebugMessage("parse <%s> as xor constraint\n", str);

   varssize = 100;
   nvars = 0;

   /* allocate buffer array for variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, varssize) );

   /* parse string */
   SCIP_CALL( SCIPparseVarsList(scip, str, vars, &nvars, varssize, &requiredsize, &endptr, ',', success) );

   if( *success )
   {
      SCIP_Real rhs;

      /* check if the size of the variable array was big enough */
      if( varssize < requiredsize )
      {
         /* reallocate memory */
         varssize = requiredsize;
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars, varssize) );

         /* parse string again with the correct size of the variable array */
         SCIP_CALL( SCIPparseVarsList(scip, str, vars, &nvars, varssize, &requiredsize, &endptr, ',', success) );
      }
      
      assert(*success);
      assert(varssize >= requiredsize);

      SCIPdebugMessage("successfully parsed %d variables\n", nvars);

      str = endptr;

      /* search for the equal symbol */
      while( *str != '=' )
         str++;
      
      if( SCIPstrToRealValue(str, &rhs, &endptr) )
      {
         assert(SCIPisZero(scip, rhs) || SCIPisEQ(scip, rhs, 1.0));
         /* create or constraint */
         SCIP_CALL( SCIPcreateConsXor(scip, cons, name, (rhs > 0.5 ? TRUE : FALSE), nvars, vars,
               initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
  
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, *cons, NULL) ) ); 
      }
      else 
         *success = FALSE;
   }

   /* free variable buffer */
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}



/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTEXEC(eventExecXor)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(event != NULL);

   consdata = (SCIP_CONSDATA*)eventdata;
   assert(consdata != NULL);

   consdata->propagated = FALSE;

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for xor constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrXor(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create event handler for events on variables */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL,
         NULL, NULL, NULL, NULL, NULL, NULL, eventExecXor,
         NULL) );

   /* create constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopyXor,
         consFreeXor, consInitXor, consExitXor, 
         consInitpreXor, consExitpreXor, consInitsolXor, consExitsolXor,
         consDeleteXor, consTransXor, consInitlpXor,
         consSepalpXor, consSepasolXor, consEnfolpXor, consEnfopsXor, consCheckXor, 
         consPropXor, consPresolXor, consRespropXor, consLockXor,
         consActiveXor, consDeactiveXor, 
         consEnableXor, consDisableXor,
         consDelvarsXor, consPrintXor, consCopyXor, consParseXor,
         conshdlrdata) );

   /* add xor constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/xor/presolpairwise",
         "should pairwise constraint comparison be performed in presolving?",
         &conshdlrdata->presolpairwise, TRUE, DEFAULT_PRESOLPAIRWISE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/xor/presolusehashing",
         "should hash table be used for detecting redundant constraints in advance", 
         &conshdlrdata->presolusehashing, TRUE, DEFAULT_PRESOLUSEHASHING, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a xor constraint x_0 xor ... xor x_{k-1} = rhs
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsXor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_Bool             rhs,                /**< right hand side of the constraint */
   int                   nvars,              /**< number of operator variables in the constraint */
   SCIP_VAR**            vars,               /**< array with operator variables of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the xor constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("xor constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, rhs, nvars, vars) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** gets number of variables in xor constraint */
int SCIPgetNVarsXor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an xor constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->nvars;
}

/** gets array of variables in xor constraint */
SCIP_VAR** SCIPgetVarsXor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an xor constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vars;
}

/** gets the right hand side of the xor constraint */
SCIP_Bool SCIPgetRhsXor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not an xor constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}
