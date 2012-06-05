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

/**@file   cons_soc.c
 * @brief  constraint handler for second order cone constraints \f$\sqrt{\gamma + \sum_{i=1}^{n} (\alpha_i\, (x_i + \beta_i))^2} \leq \alpha_{n+1}\, (x_{n+1}+\beta_{n+1})\f$
 * @author Stefan Vigerske
 * @author Marc Pfetsch
 *
 * @todo rhsvar == NULL is supported in some routines, but not everywhere
 * @todo merge square terms with same variables in presol/exitpre
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "scip/cons_soc.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_linear.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trysol.h"
#include "scip/intervalarith.h"
#include "nlpi/nlpi.h"
#include "nlpi/exprinterpret.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "soc"
#define CONSHDLR_DESC          "constraint handler for second order cone constraints"
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       -40 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY      -10 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING  SCIP_PROPTIMING_BEFORELP

#define QUADCONSUPGD_PRIORITY     10000 /**< priority of the constraint handler for upgrading of quadratic constraints */

#ifndef M_PI
#define M_PI           3.141592653589793238462643
#endif

/*
 * Data structures
 */

/** Eventdata for variable bound change events. */
struct SCIP_EventData
{
   SCIP_CONSDATA*        consdata;           /**< the constraint data */
   int                   varidx;             /**< the index of a variable on the left hand side which bound change is caught, or -1 for variable on right hand side */
   int                   filterpos;          /**< position of corresponding event in event filter */
};

/** constraint data for soc constraints */
struct SCIP_ConsData
{
   int                   nvars;              /**< number of variables on left hand side (n) */
   SCIP_VAR**            vars;               /**< variables on left hand side (x_i) */
   SCIP_Real*            coefs;              /**< coefficients for variables on left hand side (alpha_i) */
   SCIP_Real*            offsets;            /**< offsets for variables on left hand side (beta_i) */
   SCIP_Real             constant;           /**< constant on left hand side (gamma) */

   SCIP_VAR*             rhsvar;             /**< variable on right hand side (x_{n+1}) */
   SCIP_Real             rhscoeff;           /**< coefficient of square term on right hand side (alpha_{n+1}) */
   SCIP_Real             rhsoffset;          /**< offset for variable on right hand side (beta_{n+1}) */

   SCIP_NLROW*           nlrow;              /**< nonlinear row representation of constraint */

   SCIP_Real             lhsval;             /**< value of left hand side in current point */
   SCIP_Real             violation;          /**< violation of constraint in current point */

   SCIP_EVENTDATA*       lhsbndchgeventdatas;/**< eventdata for bound change events on left  hand side variables */
   SCIP_EVENTDATA        rhsbndchgeventdata; /**< eventdata for bound change event  on right hand side variable  */
   SCIP_Bool             ispropagated;       /**< does the domains need to be propagated? */
   SCIP_Bool             isapproxadded;      /**< has a linear outer approximation be added? */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_HEUR*            subnlpheur;         /**< a pointer to the subNLP heuristic, if available */
   SCIP_HEUR*            trysolheur;         /**< a pointer to the trysol heuristic, if available */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   int                   newsoleventfilterpos;/**< filter position of new solution event handler, if caught */
   SCIP_Bool             haveexprint;        /**< indicates whether an expression interpreter is available */
   SCIP_Bool             sepanlp;            /**< where linearization of the NLP relaxation solution added? */

   SCIP_Bool             glineur;            /**< is the Glineur outer approx preferred to Ben-Tal Nemirovski? */
   SCIP_Bool             doscaling;          /**< are constraint violations scaled? */
   SCIP_Bool             projectpoint;       /**< is the point in which a cut is generated projected onto the feasible set? */
   int                   nauxvars;           /**< number of auxiliary variables to use when creating a linear outer approx. of a SOC3 constraint */
   SCIP_Real             minefficacy;        /**< minimal efficacy of a cut to be added to LP in separation loop */
   SCIP_Bool             sparsify;           /**< whether to sparsify cuts */
   SCIP_Real             sparsifymaxloss;    /**< maximal loss in cut efficacy by sparsification */
   SCIP_Real             sparsifynzgrowth;   /**< growth rate of maximal allowed nonzeros in cuts in sparsification */
   SCIP_Bool             linfeasshift;       /**< whether to try to make solutions feasible in check by shifting the variable on the right hand side */
   char                  nlpform;            /**< formulation of SOC constraint in NLP */
   SCIP_Real             sepanlpmincont;     /**< minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation */
   SCIP_NODE*            lastenfolpnode;     /**< the node for which enforcement was called the last time (and some constraint was violated) */
   int                   nenfolprounds;      /**< counter on number of enforcement rounds for the current node */
};


/*
 * Local methods
 */

/** catch left hand side variable events */
static
SCIP_RETCODE catchLhsVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   varidx              /**< index of the variable which events to catch */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip      != NULL);
   assert(cons      != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata  != NULL);
   assert(varidx >= 0);
   assert(varidx < consdata->nvars);
   assert(consdata->lhsbndchgeventdatas != NULL);

   consdata->lhsbndchgeventdatas[varidx].consdata = consdata;
   consdata->lhsbndchgeventdatas[varidx].varidx   = varidx;
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vars[varidx], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, &consdata->lhsbndchgeventdatas[varidx], &consdata->lhsbndchgeventdatas[varidx].filterpos) );

   consdata->ispropagated = FALSE;

   return SCIP_OKAY;
}

/** catch right hand side variable events */
static
SCIP_RETCODE catchRhsVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons               /**< constraint for which to catch bound change events */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip      != NULL);
   assert(cons      != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata  != NULL);

   consdata->rhsbndchgeventdata.consdata = consdata;
   consdata->rhsbndchgeventdata.varidx   = -1;
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->rhsvar, SCIP_EVENTTYPE_UBTIGHTENED, eventhdlr, &consdata->rhsbndchgeventdata, &consdata->rhsbndchgeventdata.filterpos) );

   consdata->ispropagated = FALSE;

   return SCIP_OKAY;
}

/** catch variables events */
static
SCIP_RETCODE catchVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */      
   )
{
   SCIP_CONSDATA* consdata;
   int            i;

   assert(scip      != NULL);
   assert(cons      != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata  != NULL);
   assert(consdata->lhsbndchgeventdatas == NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->lhsbndchgeventdatas, consdata->nvars) );

   for( i = 0; i < consdata->nvars; ++i )
   {
      if( consdata->vars[i] != NULL )
      {
         SCIP_CALL( catchLhsVarEvents(scip, eventhdlr, cons, i) );
      }
   }

   if( consdata->rhsvar != NULL )
   {
      SCIP_CALL( catchRhsVarEvents(scip, eventhdlr, cons) );
   }

   return SCIP_OKAY;
}

/** drop left hand side variable events */
static
SCIP_RETCODE dropLhsVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   varidx              /**< index of the variable which events to catch */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip      != NULL);
   assert(cons      != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata  != NULL);
   assert(varidx >= 0);
   assert(varidx < consdata->nvars);
   assert(consdata->lhsbndchgeventdatas != NULL);
   assert(consdata->lhsbndchgeventdatas[varidx].varidx == varidx);

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vars[varidx], SCIP_EVENTTYPE_BOUNDTIGHTENED, eventhdlr, &consdata->lhsbndchgeventdatas[varidx], consdata->lhsbndchgeventdatas[varidx].filterpos) );

   return SCIP_OKAY;
}

/** drop right hand side variable events */
static
SCIP_RETCODE dropRhsVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons               /**< constraint for which to catch bound change events */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip      != NULL);
   assert(cons      != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata  != NULL);
   assert(consdata->rhsbndchgeventdata.varidx == -1);

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->rhsvar, SCIP_EVENTTYPE_UBTIGHTENED, eventhdlr, &consdata->rhsbndchgeventdata, consdata->rhsbndchgeventdata.filterpos) );

   return SCIP_OKAY;
}

/** drop variable events */
static
SCIP_RETCODE dropVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */      
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip      != NULL);
   assert(eventhdlr != NULL);
   assert(cons      != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata  != NULL);

   for( i = 0; i < consdata->nvars; ++i )
   {
      if( consdata->vars[i] != NULL )
      {
         SCIP_CALL( dropLhsVarEvents(scip, eventhdlr, cons, i) );
      }
   }

   SCIPfreeBlockMemoryArray(scip, &consdata->lhsbndchgeventdatas, consdata->nvars);

   if( consdata->rhsvar != NULL )
   {
      SCIP_CALL( dropRhsVarEvents(scip, eventhdlr, cons) );
   }

   return SCIP_OKAY;
}

/** process variable bound tightening event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{
   SCIP_CONSDATA* consdata;

   assert(scip      != NULL);
   assert(event     != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   consdata = eventdata->consdata;
   assert(consdata  != NULL);

   consdata->ispropagated = FALSE;
   /* @todo look at bounds on x_i to decide whether propagation makes sense */

   return SCIP_OKAY;
}

/** create a nonlinear row representation of the constraint and stores them in consdata */
static
SCIP_RETCODE createNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< SOC constraint handler */
   SCIP_CONS*            cons                /**< SOC constraint */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   char nlpform;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   nlpform = conshdlrdata->nlpform;
   if( nlpform == 'a' )
   {
      /* if the user let us choose, then we take 's' for "small" SOC constraints, but 'q' for large ones,
       * since the 's' form leads to nvars^2 elements in Hessian, while the 'q' form yields only n elements
       * however, if there is no expression interpreter, then the NLPI may have trouble, so we always use 'q' in this case
       */
      if( consdata->nvars < 100 && conshdlrdata->haveexprint )
         nlpform = 's';
      else
         nlpform = 'q';
   }

   switch( nlpform )
   {
   case 'e':
   {
      /* construct expression exp(\sqrt{\gamma + \sum_{i=1}^{n} (\alpha_i\, (x_i + \beta_i))^2} - alpha_{n+1}(x_{n+1} + beta_{n+1})) */

      if( consdata->nvars > 0 )
      {
         SCIP_EXPR* expr;
         SCIP_EXPR* exprterm;
         SCIP_EXPR* expr2;
         SCIP_EXPRTREE* exprtree;

         if( consdata->constant != 0.0 )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_CONST, consdata->constant) );  /* gamma */
         }
         else
         {
            exprterm = NULL;
         }

         for( i = 0; i < consdata->nvars; ++i )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, i) );  /* x_i */
            if( consdata->offsets[i] != 0.0 )
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, consdata->offsets[i]) );  /* beta_i */
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_PLUS, expr, expr2) );  /* x_i + beta_i */
            }
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_SQUARE, expr) );  /* (x_i + beta_i)^2 */
            if( consdata->coefs[i] != 1.0 )
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, consdata->coefs[i]) );  /* alpha_i */
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_MUL, expr, expr2) );  /* alpha_i * (x_i + beta_i)^2 */
            }
            if( exprterm != NULL )
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_PLUS, exprterm, expr) );
            }
            else
            {
               exprterm = expr;
            }
         }

         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_SQRT, exprterm) );  /* sqrt(gamma + sum_i (...)^2) */

         if( consdata->rhsvar != NULL )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, consdata->nvars) );  /* x_{n+1} */
            if( consdata->rhsoffset != 0.0 )
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, consdata->rhsoffset) );  /* beta_{n+1} */
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_PLUS, expr, expr2) );  /* x_{n+1} + beta_{n+1} */
            }
            if( consdata->rhscoeff != 1.0 )
            {
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, consdata->rhscoeff) );  /* alpha_{n+1} */
               SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_MUL, expr, expr2) );  /* alpha_{n+1} * (x_{n+1} + beta_{n+1}) */
            }
         }
         else
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_CONST, consdata->rhscoeff * consdata->rhsoffset) );
         }
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_MINUS, exprterm, expr) ); /* sqrt(gamma + sum_i (...)^2) - alpha_{n+1} * (x_{n+1} + beta_{n+1}) */

         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_EXP, exprterm) ); /* exp(sqrt(gamma + sum_i (...)^2) - alpha_{n+1} * (x_{n+1} + beta_{n+1})) */

         SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, exprterm, consdata->nvars+1, 0, NULL) );

         SCIP_CALL( SCIPexprtreeSetVars(exprtree, consdata->nvars, consdata->vars) );
         SCIP_CALL( SCIPexprtreeAddVars(exprtree, 1, &consdata->rhsvar) );

         SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons),
               0.0,
               0, NULL, NULL,
               0, NULL, 0, NULL,
               exprtree, -SCIPinfinity(scip), 1.0) );

         SCIP_CALL( SCIPexprtreeFree(&exprtree) );

         break;
      }
      /* if there are no left-hand-side variables, then we let the 's' case handle it */
   } /*lint -fallthrough */

   case 's':
   {
      /* construct expression \sqrt{\gamma + \sum_{i=1}^{n} (\alpha_i\, (x_i + \beta_i))^2} */

      SCIP_EXPR* expr;
      SCIP_EXPR* exprterm;
      SCIP_EXPR* expr2;
      SCIP_EXPRTREE* exprtree;
      SCIP_Real lincoef;

      if( consdata->constant != 0.0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_CONST, consdata->constant) );  /* gamma */
      }
      else
      {
         exprterm = NULL;
      }

      for( i = 0; i < consdata->nvars; ++i )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, i) );  /* x_i */
         if( consdata->offsets[i] != 0.0 )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, consdata->offsets[i]) );  /* beta_i */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_PLUS, expr, expr2) );  /* x_i + beta_i */
         }
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_SQUARE, expr) );  /* (x_i + beta_i)^2 */
         if( consdata->coefs[i] != 1.0 )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, consdata->coefs[i]) );  /* alpha_i */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_MUL, expr, expr2) );  /* alpha_i * (x_i + beta_i)^2 */
         }
         if( exprterm != NULL )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_PLUS, exprterm, expr) );
         }
         else
         {
            exprterm = expr;
         }
      }

      if( exprterm != NULL )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprterm, SCIP_EXPR_SQRT, exprterm) );  /* sqrt(gamma + sum_i (...)^2) */
         SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, exprterm, consdata->nvars, 0, NULL) );
         SCIP_CALL( SCIPexprtreeSetVars(exprtree, consdata->nvars, consdata->vars) );
      }
      else
      {
         assert(consdata->nvars == 0);
         assert(consdata->constant == 0.0);
         exprtree = NULL;
      }

      /* linear and constant part is -\alpha_{n+1} (x_{n+1}+\beta_{n+1}) */
      lincoef = -consdata->rhscoeff;
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons),
            -consdata->rhscoeff * consdata->rhsoffset,
            1, &consdata->rhsvar, &lincoef,
            0, NULL, 0, NULL,
            exprtree, -SCIPinfinity(scip), 0.0) );

      SCIP_CALL( SCIPexprtreeFree(&exprtree) );

      break;
   }

   case 'q':
   {
      /* construct quadratic form gamma + sum_{i=1}^{n} (alpha_i (x_i + beta_i))^2 <= (alpha_{n+1} (x_{n+1} + beta_{n+1})^2 */
      SCIP_QUADELEM sqrterm;
      SCIP_Real rhs;
      int rhsvarpos;

      /* create initial empty row with left hand side variables */
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
            0, NULL, NULL,
            consdata->nvars, consdata->vars, 0, NULL,
            NULL, -SCIPinfinity(scip), 0.0) );

      /* add gamma + sum_{i=1}^{n} (alpha_i x_i)^2 + 2 alpha_i beta_i x_i + beta_i^2 */
      rhs = -consdata->constant;
      for( i = 0; i < consdata->nvars; ++i )
      {
         sqrterm.idx1 = i;
         sqrterm.idx2 = i;
         sqrterm.coef = consdata->coefs[i] * consdata->coefs[i];
         SCIP_CALL( SCIPaddQuadElementToNlRow(scip, consdata->nlrow, sqrterm) );

         if( consdata->offsets[i] != 0.0 )
         {
            rhs -= consdata->offsets[i] * consdata->offsets[i];
            SCIP_CALL( SCIPaddLinearCoefToNlRow(scip, consdata->nlrow, consdata->vars[i], 2.0 * consdata->coefs[i] * consdata->offsets[i]) );
         }
      }

      /* add rhsvar to quadvars of nlrow, if not there yet */
      rhsvarpos = SCIPnlrowSearchQuadVar(consdata->nlrow, consdata->rhsvar);
      if( rhsvarpos == -1 )
      {
         SCIP_CALL( SCIPaddQuadVarToNlRow(scip, consdata->nlrow, consdata->rhsvar) );
         rhsvarpos = SCIPnlrowSearchQuadVar(consdata->nlrow, consdata->rhsvar);
         assert(rhsvarpos >= 0);
      }

      /* add -(alpha_{n+1} x_{n+1))^2 - 2 alpha_{n+1} beta_{n+1} x_{n+1} - beta_{n+1}^2 */
      sqrterm.idx1 = rhsvarpos;
      sqrterm.idx2 = rhsvarpos;
      sqrterm.coef = -consdata->rhscoeff * consdata->rhscoeff;
      SCIP_CALL( SCIPaddQuadElementToNlRow(scip, consdata->nlrow, sqrterm) );

      if( consdata->rhsoffset != 0.0 )
      {
         rhs += consdata->rhsoffset * consdata->rhsoffset;
         SCIP_CALL( SCIPaddLinearCoefToNlRow(scip, consdata->nlrow, consdata->rhsvar, -2.0 * consdata->rhscoeff * consdata->rhsoffset) );
      }

      SCIP_CALL( SCIPchgNlRowRhs(scip, consdata->nlrow, rhs) );

      break;
   }

   case 'd':
   {
      /* construct division form (gamma + sum_{i=1}^n (alpha_i(x_i+beta_i))^2)/(alpha_{n+1}(x_{n+1}+beta_{n+1})) <= alpha_{n+1}(x_{n+1}+beta_{n+1})
       */
      SCIP_EXPRTREE* exprtree;
      SCIP_EXPR* expr;
      SCIP_EXPR* nominator;
      SCIP_EXPR* denominator;
      SCIP_EXPR** exprs;
      SCIP_EXPRDATA_MONOMIAL** monomials;
      SCIP_Real lincoef;
      SCIP_Real one;
      SCIP_Real two;

      SCIP_CALL( SCIPallocBufferArray(scip, &exprs,     consdata->nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &monomials, consdata->nvars) );
      one = 1.0;
      two = 2.0;

      for( i = 0; i < consdata->nvars; ++i )
      {
         /* put x_i + beta_i into exprs[i] */
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &exprs[i], SCIP_EXPR_VARIDX, i) );
         if( consdata->offsets[i] != 0.0 )
         {
            SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &exprs[i], 1, &exprs[i], &one, consdata->offsets[i]) );
         }

         /* create monomial alpha_i^2 y_i^2, where y_i will be x_i + beta_i */
         SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip), &monomials[i], consdata->coefs[i] * consdata->coefs[i], 1, &i, &two) );
      }

      /* setup polynomial expression for gamma + sum_{i=1}^n alpha_i^2 (x_i+beta_i)^2 */
      SCIP_CALL( SCIPexprCreatePolynomial(SCIPblkmem(scip), &nominator, consdata->nvars, exprs, consdata->nvars, monomials, consdata->constant, FALSE) );  /*lint !e850 */

      SCIPfreeBufferArray(scip, &monomials);
      SCIPfreeBufferArray(scip, &exprs);

      /* setup alpha_{n+1}(x_{n+1}+beta_{n+1})
       * assert that this term is >= 0.0 (otherwise constraint is infeasible anyway) */
      assert(consdata->rhsvar != NULL);
      assert((consdata->rhscoeff >= 0.0 && !SCIPisNegative(scip, SCIPvarGetLbGlobal(consdata->rhsvar) + consdata->rhsoffset)) ||
         (consdata->rhscoeff <= 0.0 && !SCIPisPositive(scip, SCIPvarGetUbGlobal(consdata->rhsvar) + consdata->rhsoffset)));
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &denominator, SCIP_EXPR_VARIDX, consdata->nvars) );
      if( consdata->rhscoeff != 1.0 || consdata->rhsoffset != 0.0 )
      {
         SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &denominator, 1, &denominator, &consdata->rhscoeff, consdata->rhscoeff * consdata->rhsoffset) );
      }

      /* setup nominator/denominator */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_DIV, nominator, denominator) );

      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, 0, 0, NULL) );
      SCIP_CALL( SCIPexprtreeSetVars(exprtree, consdata->nvars, consdata->vars) );
      SCIP_CALL( SCIPexprtreeAddVars(exprtree, 1, &consdata->rhsvar) );

      /* linear and constant part is -\alpha_{n+1} (x_{n+1}+\beta_{n+1}) */
      lincoef = -consdata->rhscoeff;
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons),
            -consdata->rhscoeff * consdata->rhsoffset,
            1, &consdata->rhsvar, &lincoef,
            0, NULL, 0, NULL,
            exprtree, -SCIPinfinity(scip), 0.0) );

      SCIP_CALL( SCIPexprtreeFree(&exprtree) );

      break;
   }

   default:
      SCIPerrorMessage("unknown value for nlp formulation parameter\n");
      return SCIP_ERROR;
   }

   SCIPdebugMessage("created nonlinear row representation of SOC constraint\n");
   SCIPdebug( SCIPprintCons(scip, cons, NULL) );
   SCIPdebug( SCIPprintNlRow(scip, consdata->nlrow, NULL) );

   return SCIP_OKAY;
}

/** evaluates the left hand side of a SOC constraint */ 
static
SCIP_RETCODE evalLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to evaluate */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL if LP solution should be used */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      val;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->lhsval = consdata->constant;

   for( i = 0; i < consdata->nvars; ++i )
   {
      val = SCIPgetSolVal(scip, sol, consdata->vars[i]);

      if( SCIPisInfinity(scip, val) || SCIPisInfinity(scip, -val) )
      {
         consdata->lhsval = SCIPinfinity(scip);
         return SCIP_OKAY;
      }

      val = consdata->coefs[i] * (val + consdata->offsets[i]);
      consdata->lhsval += val * val;      
   }
   consdata->lhsval = sqrt(consdata->lhsval);

   return SCIP_OKAY;
}

/* computes the norm of the gradient of the SOC function */ 
static
SCIP_Real getGradientNorm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol                 /**< solution or NULL if LP solution should be used */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      g, h;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   g = 0.0;
   for( i = 0; i < consdata->nvars; ++i )
   {
      assert(!SCIPisInfinity(scip, ABS(SCIPgetSolVal(scip, sol, consdata->vars[i]))));  /*lint !e666*/

      h = SCIPgetSolVal(scip, sol, consdata->vars[i]) + consdata->offsets[i];
      h *= consdata->coefs[i] * consdata->coefs[i];
      g += h * h;
   }
   g /= consdata->lhsval * consdata->lhsval;
   g += consdata->rhscoeff * consdata->rhscoeff;

   return sqrt(g);
}

/** computes violation of a SOC constraint */
static
SCIP_RETCODE computeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to evaluate */
   SCIP_SOL*             sol,                /**< solution to evaluate, or NULL if LP solution should be used */
   SCIP_Bool             doscaling           /**< should the violation be scaled? */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real rhsval;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( evalLhs(scip, cons, sol) );

   if( SCIPisInfinity(scip, consdata->lhsval) )
   {
      /* infinity <= infinity is feasible
       * infinity <= finite value is not feasible and has violation infinity
       */
      if( (consdata->rhscoeff > 0.0 && SCIPisInfinity(scip,  SCIPgetSolVal(scip, sol, consdata->rhsvar))) ||
         ( consdata->rhscoeff < 0.0 && SCIPisInfinity(scip, -SCIPgetSolVal(scip, sol, consdata->rhsvar))) )
         consdata->violation = 0.0;
      else
         consdata->violation = SCIPinfinity(scip);
      return SCIP_OKAY;
   }

   rhsval = SCIPgetSolVal(scip, sol, consdata->rhsvar);
   if( SCIPisInfinity(scip,  rhsval) )
   {
      consdata->violation = consdata->rhscoeff > 0.0 ? 0.0 : SCIPinfinity(scip);
      return SCIP_OKAY;
   }
   if( SCIPisInfinity(scip, -rhsval) )
   {
      consdata->violation = consdata->rhscoeff < 0.0 ? 0.0 : SCIPinfinity(scip);
      return SCIP_OKAY;
   }

   consdata->violation = consdata->lhsval - consdata->rhscoeff * (SCIPgetSolVal(scip, sol, consdata->rhsvar) + consdata->rhsoffset);
   if( consdata->violation <= 0.0 )
   {
      /* constraint is not violated for sure */
      consdata->violation = 0.0;
      return SCIP_OKAY;
   }

   if( doscaling )
   {
      SCIP_Real norm = getGradientNorm(scip, cons, sol);
      if( norm > 1.0 )
         consdata->violation /= norm;
   }

   return SCIP_OKAY;
}

/** computes violations for a set of constraints */
static
SCIP_RETCODE computeViolations(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints to evaluate */
   int                   nconss,             /**< number of constraints to evaluate */
   SCIP_SOL*             sol,                /**< solution to evaluate, or NULL if LP solution should be used */
   SCIP_Bool             doscaling,          /**< should the violation be scaled? */
   SCIP_CONS**           maxviolcons         /**< a buffer to store pointer to maximal violated constraint, or NULL if of no interest */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      maxviol = 0.0;
   int            c;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);

   if( maxviolcons != NULL )
      *maxviolcons = NULL;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( computeViolation(scip, conss[c], sol, doscaling) );  /*lint !e613*/
      if( maxviolcons != NULL )
      {
         consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
         assert(consdata != NULL);
         if( consdata->violation > maxviol && SCIPisGT(scip, consdata->violation, SCIPfeastol(scip)) )
         {
            maxviol      = consdata->violation;
            *maxviolcons = conss[c];  /*lint !e613*/
         }
      }
   }

   return SCIP_OKAY;
}

/** generate supporting hyperplane in a given solution */
static
SCIP_RETCODE generateCutSol(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL for LP solution */
   SCIP_ROW**            row                 /**< place to store cut */
   )
{
   SCIP_CONSDATA* consdata;
   char           cutname[SCIP_MAXSTRLEN];
   SCIP_Real*     rowcoeff;
   SCIP_Real      rhs = 0.0;
   SCIP_Real      val;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(SCIPisPositive(scip, consdata->lhsval)); /* do not like to linearize in 0 */
   assert(!SCIPisInfinity(scip, consdata->lhsval));

   SCIP_CALL( SCIPallocBufferArray(scip, &rowcoeff, consdata->nvars) );

   for( i = 0; i < consdata->nvars; ++i )
   {
      val  = SCIPgetSolVal(scip, sol, consdata->vars[i]) + consdata->offsets[i];
      val *= consdata->coefs[i] * consdata->coefs[i];

      rowcoeff[i] = val / consdata->lhsval;

      val *= SCIPgetSolVal(scip, sol, consdata->vars[i]);
      rhs += val;
   }
   rhs /= consdata->lhsval;
   rhs -= consdata->lhsval - consdata->rhscoeff * consdata->rhsoffset;

   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "%s_linearization_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));

   SCIP_CALL( SCIPcreateEmptyRow(scip, row, cutname, -SCIPinfinity(scip), rhs, SCIPconsIsLocal(cons), FALSE, TRUE) );
   SCIP_CALL( SCIPaddVarsToRow(scip, *row, consdata->nvars, consdata->vars, rowcoeff) );
   SCIP_CALL( SCIPaddVarToRow(scip, *row, consdata->rhsvar, -consdata->rhscoeff) );

   SCIPfreeBufferArray(scip, &rowcoeff);

   return SCIP_OKAY;   
}

/** generate supporting hyperplane in a given point */
static
SCIP_RETCODE generateCutPoint(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            x,                  /**< point (lhs-vars) where to generate cut */
   SCIP_ROW**            row                 /**< place to store cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real*     rowcoeff;
   SCIP_Real      rhs = 0.0;
   SCIP_Real      lhsval;
   SCIP_Real      val;
   int            i;
   char           cutname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   lhsval = consdata->constant;
   for( i = 0; i < consdata->nvars; ++i )
   {
      assert(!SCIPisInfinity(scip, ABS(x[i])));

      val = consdata->coefs[i] * (x[i] + consdata->offsets[i]);
      lhsval += val * val;      
   }
   lhsval = sqrt(lhsval);

   if( SCIPisZero(scip, lhsval) )
   { /* do not like to linearize in 0 */
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &rowcoeff, consdata->nvars) );

   for( i = 0; i < consdata->nvars; ++i )
   {
      val  = x[i] + consdata->offsets[i];
      if( SCIPisZero(scip, val) )
      {
         rowcoeff[i] = 0.0;
         continue;
      }
      val *= consdata->coefs[i] * consdata->coefs[i];

      rowcoeff[i] = val / lhsval;

      val *= x[i];
      rhs += val;
   }
   rhs /= lhsval;
   rhs -= lhsval - consdata->rhscoeff * consdata->rhsoffset;

   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "%s_linearization_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));

   SCIP_CALL( SCIPcreateEmptyRow(scip, row, cutname, -SCIPinfinity(scip), rhs, SCIPconsIsLocal(cons), FALSE, TRUE) );
   SCIP_CALL( SCIPaddVarsToRow(scip, *row, consdata->nvars, consdata->vars, rowcoeff) );
   SCIP_CALL( SCIPaddVarToRow(scip, *row, consdata->rhsvar, -consdata->rhscoeff) );

   SCIPfreeBufferArray(scip, &rowcoeff);

   return SCIP_OKAY;   
}

/** generate supporting hyperplane w.r.t. solution projected on feasible set 
 * 
 * Instead of linearizing the SOC constraint in the given solution point, this function projects the point
 * first onto the feasible set of the SOC constraint (w.r.t. euclidean norm (scaled by alpha))
 * and linearizes the SOC constraint there.
 * The hope is that this produces somewhat tighter cuts.
 * 
 * The projection has only be computed for the case gamma = 0.
 * If gamma > 0, generateCut is called. 
 * 
 * Let \f$\hat x\f$ be sol or the LP solution if sol == NULL.
 * Then the projected point \f$\tilde x\f$ is given by
 * \f[
 *   \tilde x_i = \frac{\hat x_i + \lambda \beta_i}{1-\lambda},  \quad i=1,\ldots, n; \quad
 *   \tilde x_{n+1} = \frac{\hat x_{n+1} - \lambda \beta_{n+1}}{1+\lambda}
 * \f]
 * where
 * \f[
 *   \lambda = \frac{1-A}{1+A}, \qquad 
 *   A = \frac{\alpha_{n+1}(\hat x_{n+1}+\beta_{n+1})}{\sqrt{\sum_{i=1}^n (\alpha_i(\hat x_i+\beta_i))^2}}
 * \f]
 * 
 * If lambda is very close to 1, generateCut is called.
 * 
 * The generated cut is very similar to the unprojected form.
 * The only difference is in the right hand side, which is (in the case beta = 0) multiplied by 1/(1-lambda).
 * */
static
SCIP_RETCODE generateCutProjectedPoint(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL for LP solution */
   SCIP_ROW**            row                 /**< place to store cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real*     rowcoeff;
   SCIP_Real      rhs = 0.0;
   SCIP_Real      val;
   SCIP_Real      A, lambda;
   int            i;
   char           cutname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(SCIPisPositive(scip, consdata->lhsval)); /* do not like to linearize in 0 */
   assert(!SCIPisInfinity(scip, consdata->lhsval));

   if( !SCIPisZero(scip, consdata->constant) )
   {  /* have not thought about this case yet */
      SCIP_CALL( generateCutSol(scip, cons, sol, row) );
      return SCIP_OKAY;
   }

   A  = consdata->rhscoeff * (SCIPgetSolVal(scip, sol, consdata->rhsvar) + consdata->rhsoffset);
   A /= consdata->lhsval;

   lambda = (1.0 - A) / (1.0 + A);

   assert(!SCIPisNegative(scip, lambda)); /* otherwise A > 1, so constraint is not violated */

   SCIPdebugMessage("A = %g \t lambda = %g\n", A, lambda);

   if( SCIPisFeasEQ(scip, lambda, 1.0) )
   {  /* avoid numerical difficulties when dividing by (1-lambda) below */ 
      SCIP_CALL( generateCutSol(scip, cons, sol, row) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &rowcoeff, consdata->nvars) );

   for( i = 0; i < consdata->nvars; ++i )
   {
      val  = SCIPgetSolVal(scip, sol, consdata->vars[i]) + consdata->offsets[i];
      val *= consdata->coefs[i] * consdata->coefs[i];

      rowcoeff[i] = val / consdata->lhsval;

      val *= SCIPgetSolVal(scip, sol, consdata->vars[i]) + lambda * consdata->offsets[i];
      rhs += val;
   }
   rhs /= consdata->lhsval;
   rhs -= consdata->lhsval;
   rhs /= 1.0 - lambda;
   rhs -= consdata->rhscoeff * consdata->rhsoffset;

   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "%s_linearization_%d", SCIPconsGetName(cons), SCIPgetNLPs(scip));

   SCIP_CALL( SCIPcreateEmptyRow(scip, row, cutname, -SCIPinfinity(scip), rhs, SCIPconsIsLocal(cons), FALSE, TRUE) );
   SCIP_CALL( SCIPaddVarsToRow(scip, *row, consdata->nvars, consdata->vars, rowcoeff) );
   SCIP_CALL( SCIPaddVarToRow(scip, *row, consdata->rhsvar, -consdata->rhscoeff) );

   SCIPfreeBufferArray(scip, &rowcoeff);

   return SCIP_OKAY;   
}

/** generates sparsified supporting hyperplane */
static
SCIP_RETCODE generateSparseCut(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL for LP solution */
   SCIP_ROW**            row,                /**< place to store cut */
   SCIP_Real             minefficacy,        /**< minimal efficacy for a cut to be accepted */
   SCIP_Real             sparsifymaxloss,    /**< maximal loose of efficacy for a sparsified cut compared to constraint violation */
   SCIP_Real             nzgrowth            /**< growth factor for number of nonzeros */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real*     x;
   SCIP_Real*     dist;  /* distance to 0 */
   int*           ind;   /* indices */
   int            i;
   int            maxnz, nextmaxnz;
   SCIP_Real      efficacy;
   SCIP_Real      goodefficacy;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(SCIPisPositive(scip, consdata->lhsval)); /* do not like to linearize in 0 */
   assert(!SCIPisInfinity(scip, consdata->lhsval));

   if( consdata->nvars <= 3 )
   {
      SCIP_CALL( generateCutSol(scip, cons, sol, row) );
      return SCIP_OKAY;
   }

   goodefficacy = MAX((1.0-sparsifymaxloss) * consdata->violation, minefficacy);

   SCIP_CALL( SCIPallocBufferArray(scip, &x,    consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dist, consdata->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ind,  consdata->nvars) );

   SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nvars, consdata->vars, x) );
   /* distance to "-offset" * alpha_i^2 should indicate loss when moving refpoint to x[i] = -offset[i] */
   for( i = 0; i < consdata->nvars; ++i )
   {
      ind[i] = i;
      dist[i]  = ABS(x[i] + consdata->offsets[i]);
      dist[i] *= consdata->coefs[i] * consdata->coefs[i];
   }

   /* sort variables according to dist */
   SCIPsortRealInt(dist, ind, consdata->nvars);

   maxnz = 2;
   /* set all except biggest maxnz entries in x to -offset */
   for( i = 0; i < consdata->nvars - maxnz; ++i )
      x[ind[i]] = -consdata->offsets[i];

   do
   {
      /* @todo speed up a bit by computing efficacy of new cut from efficacy of old cut
       * generate row only if efficient enough */
      SCIP_CALL( generateCutPoint(scip, cons, x, row) );

      if( *row != NULL )
      {
         efficacy = -SCIPgetRowSolFeasibility(scip, *row, sol) / SCIPgetRowMaxCoef(scip, *row);

         if( SCIPisGT(scip, efficacy, goodefficacy) ||
            (maxnz >= consdata->nvars && SCIPisGT(scip, efficacy, minefficacy)) )
         {
            /* cut cuts off solution and is efficient enough */
            SCIPdebugMessage("accepted cut with %d of %d nonzeros, efficacy = %g\n", maxnz, consdata->nvars, efficacy);
            break;
         }
         SCIP_CALL( SCIPreleaseRow(scip, row) );
      }

      if( maxnz >= consdata->nvars )
      { /* cut also not efficient enough if generated in original refpoint (that's bad) */
         break;
      }

      nextmaxnz = (int)(nzgrowth * maxnz);
      if( nextmaxnz == consdata->nvars - 1)
         nextmaxnz = consdata->nvars;
      else if( nextmaxnz == maxnz )
         ++nextmaxnz;

      /* restore entries of x that are nonzero in next attempt */
      for( i = MAX(0, consdata->nvars - nextmaxnz); i < consdata->nvars - maxnz; ++i )
         x[ind[i]] = SCIPgetSolVal(scip, sol, consdata->vars[ind[i]]);

      maxnz = nextmaxnz;
   } while( TRUE );  /*lint !e506*/

   SCIPfreeBufferArray(scip, &x);
   SCIPfreeBufferArray(scip, &dist);
   SCIPfreeBufferArray(scip, &ind);

   return SCIP_OKAY;
}

/** separates a point, if possible */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of constraints that seem to be useful */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL for LP solution */
   SCIP_Bool             addweakcuts,        /**< whether also weak (only slightly violated) cuts should be added in a nonconvex constraint */
   SCIP_Bool*            success             /**< buffer to store whether the point was separated */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          minefficacy;
   int                c;
   SCIP_ROW*          row;

   assert(scip    != NULL);
   assert(conss   != NULL || nconss == 0);
   assert(nusefulconss <= nconss);
   assert(success != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *success = FALSE;

   minefficacy = addweakcuts ? SCIPfeastol(scip) : conshdlrdata->minefficacy;

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      if( SCIPisGT(scip, consdata->violation, SCIPfeastol(scip)) && !SCIPisInfinity(scip, consdata->violation) )
      {
         row = NULL;

         /* generate cut */
         if( conshdlrdata->sparsify )
         {
            SCIP_CALL( generateSparseCut(scip, conss[c], sol, &row, minefficacy, conshdlrdata->sparsifymaxloss, conshdlrdata->sparsifynzgrowth) );  /*lint !e613*/
         }  
         else if( conshdlrdata->projectpoint )
         {
            SCIP_Real efficacy;

            SCIP_CALL( generateCutProjectedPoint(scip, conss[c], sol, &row) );  /*lint !e613*/

            efficacy = -SCIPgetRowSolFeasibility(scip, row, sol) / SCIPgetRowMaxCoef(scip, row);
            if( SCIPisLE(scip, efficacy, minefficacy) )
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
         }
         else
         {
            SCIP_Real efficacy;

            SCIP_CALL( generateCutSol(scip, conss[c], sol, &row) );  /*lint !e613*/

            efficacy = -SCIPgetRowSolFeasibility(scip, row, sol) / SCIPgetRowMaxCoef(scip, row);
            if( SCIPisLE(scip, efficacy, minefficacy) )
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
         }

         if( row == NULL ) /* failed to generate (efficient enough) cut */
            continue;

         /* cut cuts off solution and efficient enough */
         SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );
         SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );  /*lint !e613*/
         *success = TRUE;
         SCIPdebugMessage("added cut with efficacy %g\n", SCIPgetCutEfficacy(scip, sol, row));

         SCIP_CALL( SCIPreleaseRow (scip, &row) );
      }

      /* enforce only useful constraints
       * others are only checked and enforced if we are still feasible or have not found a separating cut yet
       */ 
      if( c >= nusefulconss && *success )
         break;
   }

   return SCIP_OKAY;
}

/** adds linearizations cuts for convex constraints w.r.t. a given reference point to cutpool and sepastore
 * if separatedlpsol is not NULL, then cuts that separate the LP solution are added to the sepastore too
 */
static
SCIP_RETCODE addLinearizationCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< quadratic constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             ref,                /**< reference point where to linearize, or NULL for LP solution */
   SCIP_Bool*            separatedlpsol,     /**< buffer to store whether a cut that separates the current LP solution was found, or NULL if not of interest */
   SCIP_Real             minefficacy         /**< minimal efficacy of a cut when checking for separation of LP solution */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_ROW* row;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   if( separatedlpsol != NULL )
      *separatedlpsol = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613 */

      if( SCIPconsIsLocal(conss[c]) )  /*lint !e613 */
         continue;

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613 */
      assert(consdata != NULL);

      SCIP_CALL( evalLhs(scip, conss[c], ref) );  /*lint !e613 */
      if( !SCIPisPositive(scip, consdata->lhsval) || SCIPisInfinity(scip, consdata->lhsval) )
      {
         SCIPdebugMessage("skip adding linearization for <%s> since lhs is %g\n", SCIPconsGetName(conss[c]), consdata->lhsval);  /*lint !e613 */
         continue;
      }

      SCIP_CALL( generateCutSol(scip, conss[c], ref, &row) );  /*lint !e613 */

      if( row == NULL )
         continue;

      assert(!SCIProwIsLocal(row));

      /* if caller wants, then check if cut separates LP solution and add to sepastore if so */
      if( separatedlpsol != NULL )
      {
         SCIP_Real feasibility;
         SCIP_Real norm;

         feasibility = SCIPgetRowLPFeasibility(scip, row);
         norm = SCIPgetRowMaxCoef(scip, row);

         if( -feasibility / MAX(1.0, norm) >= minefficacy )
         {
            *separatedlpsol = TRUE;
            SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE) );
         }
      }

      SCIP_CALL( SCIPaddPoolCut(scip, row) );

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   return SCIP_OKAY;
}

/** processes the event that a new primal solution has been found */
static
SCIP_DECL_EVENTEXEC(processNewSolutionEvent)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS**    conss;
   int            nconss;
   SCIP_SOL*      sol;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   assert((SCIPeventGetType(event) & SCIP_EVENTTYPE_SOLFOUND) != 0);

   conshdlr = (SCIP_CONSHDLR*)eventdata;

   nconss = SCIPconshdlrGetNConss(conshdlr);

   if( nconss == 0 )
      return SCIP_OKAY;

   sol = SCIPeventGetSol(event);
   assert(sol != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* we are only interested in solution coming from some heuristic other than trysol, but not from the tree
    * the reason for ignoring trysol solutions is that they may come from an NLP solve in sepalp, where we already added linearizations,
    * or are from the tree, but postprocessed via proposeFeasibleSolution
    */
   if( SCIPsolGetHeur(sol) == NULL || SCIPsolGetHeur(sol) == conshdlrdata->trysolheur )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);
   assert(conss != NULL);

   SCIPdebugMessage("caught new sol event %x from heur <%s>; have %d conss\n", SCIPeventGetType(event), SCIPheurGetName(SCIPsolGetHeur(sol)), nconss);

   SCIP_CALL( addLinearizationCuts(scip, conshdlr, conss, nconss, sol, NULL, 0.0) );

   return SCIP_OKAY;
}

/** removes fixed variables, replace aggregated and negated variables
 *
 * repeats replacements until no further change is found;
 * takes care of capture/release and locks, but not of variable events (assumes that var events are not caught yet) 
 */
static
SCIP_RETCODE presolveRemoveFixedVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for signpower constraints */
   SCIP_CONS*            cons,               /**< constraint */
   int*                  ndelconss,          /**< counter for number of deleted constraints */
   int*                  nupgdconss,         /**< counter for number of upgraded constraints */
   int*                  nchgbds,            /**< counter for number of bound changes */
   int*                  nfixedvars,         /**< counter for number of fixed variables */
   SCIP_Bool*            iscutoff,           /**< to store whether constraint cannot be satisfied */
   SCIP_Bool*            isdeleted           /**< to store whether constraint is redundant and can be deleted */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool      havechange;
   SCIP_Bool      haveremovedvar;
   int            i;
   SCIP_VAR*      x;
   SCIP_Real      coef;
   SCIP_Real      offset;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(iscutoff != NULL);
   assert(isdeleted != NULL);

   *iscutoff  = FALSE;
   *isdeleted = FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMessage("remove fixed variables from constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );

   havechange     = FALSE;
   haveremovedvar = FALSE;

   /* process variables on left hand side */
   for( i = 0; i < consdata->nvars; ++i )
   {
      x = consdata->vars[i];
      assert(x != NULL);
      assert(SCIPvarGetStatus(x) != SCIP_VARSTATUS_ORIGINAL);

      if( SCIPvarIsActive(x) || SCIPvarGetStatus(x) == SCIP_VARSTATUS_MULTAGGR )
         continue;

      havechange = TRUE;

      /* drop variable event and unlock and release variable */
      SCIP_CALL( dropLhsVarEvents(scip, conshdlrdata->eventhdlr, cons, i) );
      SCIP_CALL( SCIPunlockVarCons(scip, x, cons, TRUE, TRUE) );
      SCIP_CALL( SCIPreleaseVar(scip, &consdata->vars[i]) );

      coef = 1.0;
      offset = consdata->offsets[i];
      SCIP_CALL( SCIPvarGetProbvarSum(&x, &coef, &offset) );

      SCIPdebugMessage("  lhs term at position %d is replaced by %g * <%s> + %g\n",
         i, coef, SCIPvarGetName(x), offset);

      /* if variable has been fixed, add (alpha*offset)^2 to gamma and continue */
      if( coef == 0.0 || x == NULL )
      {
         consdata->constant  += consdata->coefs[i] * consdata->coefs[i] * offset * offset;
         consdata->offsets[i] = 0.0;
         haveremovedvar = TRUE;
         continue;
      }

      assert(SCIPvarIsActive(x) || SCIPvarGetStatus(x) == SCIP_VARSTATUS_MULTAGGR);

      /* replace coefs[i] * (vars[i] + offsets[i]) by coefs[i]*coef * (x + offsets[i]/coef) */
      consdata->offsets[i] = offset;
      if( coef != 1.0 )
      {
         consdata->coefs[i]    = REALABS(coef * consdata->coefs[i]);
         consdata->offsets[i] /= coef;
      }
      consdata->vars[i] = x;

      /* capture and lock new variable, catch variable events */
      SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[i]) );
      SCIP_CALL( SCIPlockVarCons(scip, consdata->vars[i], cons, TRUE, TRUE) );
      SCIP_CALL( catchLhsVarEvents(scip, conshdlrdata->eventhdlr, cons, i) );
   }

   /* process variable on right hand side */
   x = consdata->rhsvar;
   assert(x != NULL);
   if( !SCIPvarIsActive(x) && SCIPvarGetStatus(x) != SCIP_VARSTATUS_MULTAGGR )
   {
      havechange = TRUE;

      /* drop variable event and unlock and release variable */
      SCIP_CALL( dropRhsVarEvents(scip, conshdlrdata->eventhdlr, cons) );
      SCIP_CALL( SCIPreleaseVar(scip, &consdata->rhsvar) );
      SCIP_CALL( SCIPunlockVarCons(scip, x, cons, consdata->rhscoeff > 0.0, consdata->rhscoeff < 0.0) );

      coef = 1.0;
      offset = 0.0;
      SCIP_CALL( SCIPvarGetProbvarSum(&x, &coef, &offset) );

      SCIPdebugMessage("  rhs variable is replaced by %g * <%s> + %g\n", coef, SCIPvarGetName(x), offset);

      if( coef == 0.0 || x == NULL )
      {
         /* if variable has been fixed, add offset to rhsoffset */
         consdata->rhsoffset += offset;
      }
      else
      {
         /* replace rhscoef * (rhsvar + rhsoffset) by rhscoef*coef * (x + offset/coef + rhsoffset/coef) */
         assert(SCIPvarIsActive(x) || SCIPvarGetStatus(x) == SCIP_VARSTATUS_MULTAGGR);

         consdata->rhsoffset = (consdata->rhsoffset + offset) / coef;
         consdata->rhscoeff *= coef;
         consdata->rhsvar = x;

         /* capture and lock new variable, catch variable events */
         SCIP_CALL( SCIPcaptureVar(scip, consdata->rhsvar) );
         SCIP_CALL( SCIPlockVarCons(scip, consdata->rhsvar, cons, consdata->rhscoeff > 0.0, consdata->rhscoeff < 0.0) );
         SCIP_CALL( catchRhsVarEvents(scip, conshdlrdata->eventhdlr, cons) );
      }
   }

   if( !havechange )
      return SCIP_OKAY;

   /* free nonlinear row representation */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* if a variable has been removed, close gaps in vars array */
   if( haveremovedvar )
   {
      int oldnvars;

      /* due to the realloc of the block memory below and the way we store the eventdata in consdata, we best drop all events here and catch them again below */
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );

      oldnvars = consdata->nvars;
      for( i = 0; i < consdata->nvars; ++i )
      {
         /* forget about empty places at end of vars array */
         while( consdata->nvars && consdata->vars[consdata->nvars-1] == NULL )
            --consdata->nvars;

         /* all variables at index >= i have been removed */
         if( i == consdata->nvars )
            break;

         if( consdata->vars[i] != NULL )
            continue;

         /* move variable from position nvars-1 to position i */

         assert(consdata->nvars >= 1);
         assert(consdata->vars[consdata->nvars-1] != NULL);

         consdata->vars[i]    = consdata->vars[consdata->nvars-1];
         consdata->offsets[i] = consdata->offsets[consdata->nvars-1];
         consdata->coefs[i]   = consdata->coefs[consdata->nvars-1];

         --consdata->nvars;
      }

      assert(consdata->nvars < oldnvars);

      /* shrink arrays in consdata */
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars,    oldnvars, consdata->nvars) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->offsets, oldnvars, consdata->nvars) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->coefs,   oldnvars, consdata->nvars) );

      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   SCIPdebugMessage("\t-> ");
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );

   if( consdata->nvars == 0 )
   { /* all variables on left hand size have been removed, remaining constraint is sqrt(gamma) <= ... */
      assert(!SCIPisNegative(scip, consdata->constant));
      if( consdata->rhsvar == NULL )
      { /* also rhsvar has been removed, remaining constraint is sqrt(gamma) <= rhscoeff * rhsoffset */
         if( SCIPisFeasLE(scip, sqrt(consdata->constant), consdata->rhscoeff*consdata->rhsoffset) )
         {
            SCIPdebugMessage("remove redundant constraint <%s> after fixing all variables\n", SCIPconsGetName(cons));
         }
         else
         {
            SCIPdebugMessage("found problem infeasible after fixing all variables in <%s>\n", SCIPconsGetName(cons));
            *iscutoff = TRUE;
         }
         ++*ndelconss;
      }
      else if( !SCIPvarIsActive(consdata->rhsvar) )
      { /* remaining constraint is sqrt(gamma) - rhscoeff * rhsoffset <= rhscoeff * rhsvar, and rhsvar is probably multi-aggregated */
         SCIP_CONS* lincons;

         SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons), 1, &consdata->rhsvar, &consdata->rhscoeff,
               sqrt(consdata->constant) - consdata->rhscoeff * consdata->rhsoffset, SCIPinfinity(scip),
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
               SCIPconsIsStickingAtNode(cons)) );
         SCIP_CALL( SCIPaddCons(scip, lincons) );
         SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
         ++*nupgdconss;
      }
      else if( consdata->rhscoeff > 0.0 )
      { /* remaining constraint is sqrt(gamma) / rhscoeff - rhsoffset <= rhsvar */
         SCIP_Bool tightened;
         SCIP_CALL( SCIPtightenVarLb(scip, consdata->rhsvar, sqrt(consdata->constant) / consdata->rhscoeff - consdata->rhsoffset, TRUE, iscutoff, &tightened) );
         if( *iscutoff )
         {
            SCIPdebugMessage("found problem infeasible after fixing all lhs variables in <%s> and tightening lower bound of rhs var\n", SCIPconsGetName(cons));
         }
         else if( tightened )
         {
            SCIPdebugMessage("remove redundant constraint <%s> after fixing all lhs variables and tightening lower bound of rhs var\n", SCIPconsGetName(cons));
            ++*nchgbds;
         }
         else
         {
            SCIPdebugMessage("remove redundant constraint <%s> after fixing all lhs variables\n", SCIPconsGetName(cons));
         }
         ++*ndelconss;
      }
      else
      { /* remaining constraint is sqrt(gamma) / rhscoeff - rhsoffset >= rhsvar */
         SCIP_Bool tightened;
         SCIP_CALL( SCIPtightenVarUb(scip, consdata->rhsvar, sqrt(consdata->constant) / consdata->rhscoeff - consdata->rhsoffset, TRUE, iscutoff, &tightened) );
         if( *iscutoff )
         {
            SCIPdebugMessage("found problem infeasible after fixing all lhs variables in <%s> and tightening upper bound of rhs var\n", SCIPconsGetName(cons));
         }
         else if( tightened )
         {
            SCIPdebugMessage("remove redundant constraint <%s> after fixing all lhs variables and tightening upper bound of rhs var\n", SCIPconsGetName(cons));
            ++*nchgbds;
         }
         else
         {
            SCIPdebugMessage("remove redundant constraint <%s> after fixing all lhs variables\n", SCIPconsGetName(cons));
         }
         ++*ndelconss;
      }
      SCIP_CALL( SCIPdelCons(scip, cons) );
      *isdeleted = TRUE;
      return SCIP_OKAY;
   }

   if( consdata->rhsvar == NULL )
   { /* constraint becomes sum_i (alpha_i*(x_i+beta_i))^2 <= (rhscoeff*rhsoffset)^2 - gamma */
      if( consdata->nvars > 1 )
      { /* upgrade to quadratic constraint */
         SCIP_CONS* quadcons;
         SCIP_QUADVARTERM* quadvarterms;
         SCIP_Real  rhs;

         SCIP_CALL( SCIPallocBufferArray(scip, &quadvarterms, consdata->nvars) );
         BMSclearMemoryArray(quadvarterms, consdata->nvars);
         rhs = consdata->rhscoeff * consdata->rhsoffset;
         rhs = rhs * rhs - consdata->constant;

         for( i = 0; i < consdata->nvars; ++i )
         {
            quadvarterms[i].var = consdata->vars[i];
            quadvarterms[i].sqrcoef = consdata->coefs[i] * consdata->coefs[i];
            if( consdata->offsets[i] != 0.0 )
            {
               quadvarterms[i].lincoef = 2 * consdata->offsets[i] * quadvarterms[i].sqrcoef;
               rhs -= quadvarterms[i].sqrcoef * consdata->offsets[i]*consdata->offsets[i];
            }
         }

         assert(!SCIPconsIsStickingAtNode(cons));
         SCIP_CALL( SCIPcreateConsQuadratic2(scip, &quadcons, SCIPconsGetName(cons), 0, NULL, NULL,
               consdata->nvars, quadvarterms, 0, NULL, -SCIPinfinity(scip), rhs,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
         SCIP_CALL( SCIPaddCons(scip, quadcons) );
         SCIPdebugMessage("upgraded <%s> to quadratic constraint: ", SCIPconsGetName(cons));
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, quadcons, NULL) ) );

         SCIP_CALL( SCIPreleaseCons(scip, &quadcons) );

         SCIPfreeBufferArray(scip, &quadvarterms);

         ++*nupgdconss;
      }
      else if( !SCIPvarIsActive(consdata->vars[0]) )
      { /* constraint is |alpha*(x+beta)| <= sqrt((rhscoeff*rhsoffset)^2 - gamma), but x is probably multaggr. -> turn into ranged linear constraint */
         SCIP_CONS* lincons;

         /* create constraint alpha*x <=  sqrt((rhscoeff*rhsoffset)^2 - gamma) - alpha*beta
          *                   alpha*x >= -sqrt((rhscoeff*rhsoffset)^2 - gamma) - alpha*beta */
         SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons), 1, &consdata->vars[0], &consdata->coefs[0],
               -sqrt(consdata->rhscoeff * consdata->rhscoeff * consdata->rhsoffset * consdata->rhsoffset - consdata->constant) - consdata->coefs[0] * consdata->offsets[0],
               +sqrt(consdata->rhscoeff * consdata->rhscoeff * consdata->rhsoffset * consdata->rhsoffset - consdata->constant) - consdata->coefs[0] * consdata->offsets[0],
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
               SCIPconsIsStickingAtNode(cons)) );
         SCIP_CALL( SCIPaddCons(scip, lincons) );
         SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

         ++*nupgdconss;
      }
      else
      { /* constraint is |alpha*(x+beta)| <= sqrt((rhscoeff*rhsoffset)^2 - gamma) -> propagate bounds */
         SCIP_Bool tightened;
         SCIP_Real rhs;
         assert(consdata->nvars == 1); /* case == 0 handled before */
         rhs = consdata->rhscoeff * consdata->rhsoffset;
         rhs = rhs * rhs;
         if( SCIPisNegative(scip, rhs - consdata->constant) )
         { /* take this as infeasible */
            SCIPdebugMessage("found problem infeasible after fixing rhs and all except one lhs variables in <%s>\n", SCIPconsGetName(cons));
            *iscutoff = TRUE;
         }
         else
         {
            rhs -= consdata->constant;
            rhs  = rhs < 0.0 ? 0.0 : sqrt(rhs);

            if( SCIPisZero(scip, rhs) )
            { /* constraint is x = -beta */
               SCIP_CALL( SCIPfixVar(scip, consdata->vars[0], -consdata->offsets[0], iscutoff, &tightened) );
               if( *iscutoff )
               {
                  SCIPdebugMessage("found problem infeasible after fixing rhs and all except one lhs variables and fixing remaining lhs var in <%s>\n", SCIPconsGetName(cons));
               }
               else if( tightened )
               {
                  SCIPdebugMessage("remove redundant constraint <%s> after fixing rhs and all except one lhs variables and fixing remaining lhs var\n", SCIPconsGetName(cons));
                  ++*nfixedvars;
               }
               else
               {
                  SCIPdebugMessage("remove redundant constraint <%s> after fixing rhs and all except one lhs variables and fixing remaining lhs var\n", SCIPconsGetName(cons));
               }
            }
            else
            { /* constraint is -rhs/|alpha| - beta <= x <= rhs/|alpha| - beta */
               rhs /= ABS(consdata->coefs[0]);
               SCIP_CALL( SCIPtightenVarLb(scip, consdata->vars[0], -rhs - consdata->offsets[0], TRUE, iscutoff, &tightened) );
               if( *iscutoff )
               {
                  SCIPdebugMessage("found problem infeasible after fixing rhs and all except one lhs variables and tightening lower bound of remaining lhs var in <%s>\n", SCIPconsGetName(cons));
               }
               else
               {
                  if( tightened )
                     ++*nchgbds;
                  SCIP_CALL( SCIPtightenVarUb(scip, consdata->vars[0], rhs - consdata->offsets[0], TRUE, iscutoff, &tightened) );
                  if( *iscutoff )
                  {
                     SCIPdebugMessage("found problem infeasible after fixing rhs and all except one lhs variables and tightening upper bound of remaining lhs var in <%s>\n", SCIPconsGetName(cons));
                  }
                  else if( tightened )
                     ++*nchgbds;
               }
               if( !*iscutoff )
               {
                  SCIPdebugMessage("remove redundant constraint <%s> after fixing rhs and all except one lhs variables and tightening bounds on remaining lhs var\n", SCIPconsGetName(cons));
               }
            }
         }
         ++*ndelconss;
      }
      *isdeleted = TRUE;
      SCIP_CALL( SCIPdelCons(scip, cons) );
      return SCIP_OKAY;
   }

   if( consdata->nvars == 1 && SCIPisZero(scip, consdata->constant) )
   { /* one variable on lhs left and no constant, constraint becomes |alpha*(x+beta)| <= rhscoef*(rhsvar+rhsoffset) -> upgrade to two linear constraints */
      SCIP_CONS* lincons;
      SCIP_VAR*  vars[2];
      SCIP_Real  coefs[2];
      SCIP_Real  rhs;
      assert(consdata->rhsvar != NULL); /* case == NULL has been handled before */

      vars[0] = consdata->vars[0];
      vars[1] = consdata->rhsvar;
      coefs[0] = consdata->coefs[0];
      coefs[1] = -consdata->rhscoeff;
      rhs = consdata->rhscoeff * consdata->rhsoffset - coefs[0] * consdata->offsets[0];

      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons), 2, vars, coefs, -SCIPinfinity(scip), rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

      coefs[0] = -coefs[0];
      rhs = consdata->rhscoeff * consdata->rhsoffset - coefs[0] * consdata->offsets[0];

      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons), 2, vars, coefs, -SCIPinfinity(scip), rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

      SCIPdebugMessage("upgraded <%s> to two linear constraint\n", SCIPconsGetName(cons));

      ++*nupgdconss;
      SCIP_CALL( SCIPdelCons(scip, cons) );
      *isdeleted = TRUE;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}


/** adds the linear outer-approximation of Glineur et.al. for a SOC constraint of dimension 3
 * 
 * Input is the data for a constraint \f$\sqrt{(\alpha_1(x_1+offset1))^2 + (\alpha_2(x_2+offset2))^2) \leq \alpha_3(x_3+offset3)}\f$.
 * Here constant >= 0.0, alpha3 > 0.0, and the lower bound of x3 >= -offset3.
 * Also x2 = NULL is allowed, in which case the second term is assumed to be constant, and offset2 != 0 is needed.
 */
static
SCIP_RETCODE presolveCreateGlineurApproxDim3(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< original constraint */
   SCIP_VAR*             x1,                 /**< variable x1 */
   SCIP_VAR*             x2,                 /**< variable x2 */
   SCIP_VAR*             x3,                 /**< variable x3 */
   SCIP_Real             alpha1,             /**< coefficient of x1 */
   SCIP_Real             alpha2,             /**< coefficient of x2 */
   SCIP_Real             alpha3,             /**< coefficient of x3 */
   SCIP_Real             offset1,            /**< offset of x1 */
   SCIP_Real             offset2,            /**< offset of x2 */
   SCIP_Real             offset3,            /**< offset of x3 */
   int                   N,                  /**< size of linear approximation, need to be >= 1 */
   const char*           basename,           /**< string to use for building variable and constraint names */
   int*                  naddconss           /**< buffer where to add the number of added constraints */
   )
{
   SCIP_CONS*     lincons;
   SCIP_VAR*      vars[3];
   SCIP_Real      vals[3];
   char           varname[255];
   char           linname[255];
   int            i;
   SCIP_VAR**     avars;
   SCIP_VAR**     bvars;
   SCIP_Real      val;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(x1   != NULL);
   assert(x2   != NULL || !SCIPisZero(scip, offset2));
   assert(x3   != NULL);
   assert(SCIPisPositive(scip, alpha3));
   assert(SCIPisGE(scip, SCIPconsIsLocal(cons) ? SCIPvarGetLbLocal(x3) : SCIPvarGetLbGlobal(x3), -offset3));
   assert(basename != NULL);
   assert(N >= 1);
   assert(naddconss != NULL);

   SCIPdebugMessage("Creating linear Glineur outer-approximation for <%s>.\n", basename);
   SCIPdebugMessage("sqr(%g(%s+%g)) + sqr(%g(%s+%g)) <= sqr(%g(%s+%g)).\n", 
      alpha1, SCIPvarGetName(x1), offset1, alpha2, x2 ? SCIPvarGetName(x2) : "0", offset2, alpha3, SCIPvarGetName(x3), offset3
      );

   SCIP_CALL( SCIPallocBufferArray(scip, &avars, N+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bvars, N+1) );

   /* create additional variables; we do not use avars[0] and bvars[0] */
   for( i = 1; i <= N; ++i )
   {
      (void) SCIPsnprintf(varname, 255, "soc#%s_a%d", basename, i);
      SCIP_CALL( SCIPcreateVar(scip, &avars[i], varname, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, 
            SCIP_VARTYPE_CONTINUOUS, SCIPconsIsInitial(cons), FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, avars[i]) );

      (void) SCIPsnprintf(varname, 255, "soc#%s_b%d", basename, i);
      SCIP_CALL( SCIPcreateVar(scip, &bvars[i], varname, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, SCIPconsIsInitial(cons), FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, bvars[i]) );
   }

   /* create linear constraints for the first case
    * cos(pi) = -1, sin(pi) = 0
    * -> a_1  = - alpha1 (x1 + offset1)    ->  -alpha1*x1 - a_1  =  alpha1*offset1 
    * -> b_1 >= | alpha2 (x2 + offset2) |  ->   alpha2*x2 - b_1 <= -alpha2*offset2
    *                                           alpha2*x2 + b_1 >= -alpha2*offset2
    */

   vars[0] = x1;
   vals[0] = -alpha1;
   vars[1] = avars[1];
   vals[1] = -1.0;

   (void) SCIPsnprintf(linname, 255, "soc#%s#a%d", basename, 0);
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, alpha1*offset1, alpha1*offset1,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
         SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, lincons, NULL) ) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   ++*naddconss;

   if( x2 != NULL )
   {
      vars[0] = x2;
      vals[0] = alpha2;
      vars[1] = bvars[1];
      vals[1] = -1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#b%d", basename, 0);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, -SCIPinfinity(scip), -alpha2*offset2,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, lincons, NULL) ) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++*naddconss;

      vars[0] = x2;
      vals[0] = alpha2;
      vars[1] = bvars[1];
      vals[1] = 1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#B%d", basename, 0);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, -alpha2*offset2, SCIPinfinity(scip),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, lincons, NULL) ) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++*naddconss;
   }
   else
   { /* x2 == NULL ->  b_1 >= |alpha2*offset2| */
      SCIP_Bool infeas;
      SCIP_Bool tightened;
      SCIP_CALL( SCIPtightenVarLb(scip, bvars[1], ABS(alpha2 * offset2), TRUE, &infeas, &tightened) );
      if( infeas == TRUE )
      {
         SCIPwarningMessage("creating glineur outer approximation of SOC3 constraint found problem infeasible.\n");
      }
   }

   /* create intermediate linear constraints */
   val = M_PI;
   for( i = 1; i < N; ++i )
   {
      val /= 2.0;

      vars[0] = avars[i];
      vals[0] = cos(val);
      vars[1] = bvars[i];
      vals[1] = sin(val);
      vars[2] = avars[i+1];
      vals[2] = -1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#a%d", basename, i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, 0.0, 0.0,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, lincons, NULL) ) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++*naddconss;

      vars[0] = avars[i];
      vals[0] = -sin(val);
      vars[1] = bvars[i];
      vals[1] = cos(val);
      vars[2] = bvars[i+1];
      vals[2] = -1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#b%d", basename, i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, -SCIPinfinity(scip), 0.0,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, lincons, NULL) ) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++*naddconss;

      vars[0] = avars[i];
      vals[0] = -sin(val);
      vars[1] = bvars[i];
      vals[1] = cos(val);
      vars[2] = bvars[i+1];
      vals[2] = 1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#B%d", basename, i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, 0.0, SCIPinfinity(scip),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, lincons, NULL) ) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++*naddconss;
   }

   /* create last linear constraint */
   val /= 2.0;
   vars[0] = avars[N];
   vals[0] = -cos(val);
   vars[1] = bvars[N];
   vals[1] = -sin(val);
   vars[2] = x3;
   vals[2] = alpha3;

   (void) SCIPsnprintf(linname, 255, "soc#%s#a%d", basename, N);
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, -alpha3*offset3, -alpha3*offset3,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
         SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, lincons, NULL) ) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   ++*naddconss;

   for( i = 1; i <= N; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &avars[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &bvars[i]) );
   }
   SCIPfreeBufferArray(scip, &avars);
   SCIPfreeBufferArray(scip, &bvars);

   return SCIP_OKAY;
}

/** adds the linear outer-approximation of Ben-Tal and Nemirovski for a SOC constraint of dimension 3
 * 
 * Input is the data for a constraint \f$\sqrt{constant + (\alpha_1(x_1+offset1))^2 + (\alpha_2(x_2+offset2))^2) \leq \alpha_3(x_3+offset3)}\f$.
 * Here constant >= 0.0, alpha3 > 0.0, and the lower bound of x3 >= -offset3.
 * Also x2 = NULL is allowed, in which case the second term is assumed to be constant, and offset2 != 0 is needed.
 * */
static
SCIP_RETCODE presolveCreateBenTalNemirovskiApproxDim3(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< original constraint */
   SCIP_VAR*             x1,                 /**< variable x1 */
   SCIP_VAR*             x2,                 /**< variable x2 */
   SCIP_VAR*             x3,                 /**< variable x3 */
   SCIP_Real             alpha1,             /**< coefficient of x1 */
   SCIP_Real             alpha2,             /**< coefficient of x2 */
   SCIP_Real             alpha3,             /**< coefficient of x3 */
   SCIP_Real             offset1,            /**< offset of x1 */
   SCIP_Real             offset2,            /**< offset of x2 */
   SCIP_Real             offset3,            /**< offset of x3 */
   int                   N,                  /**< size of linear approximation, need to be >= 1 */
   const char*           basename,           /**< string to use for building variable and constraint names */
   int*                  naddconss           /**< buffer where to add the number of added constraints */
   )
{
   SCIP_CONS*     lincons;
   SCIP_VAR*      vars[3];
   SCIP_Real      vals[3];
   char           varname[255];
   char           linname[255];
   int            i;
   SCIP_VAR**     avars;
   SCIP_VAR**     bvars;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(x1   != NULL);
   assert(x2   != NULL || !SCIPisZero(scip, offset2));
   assert(x3   != NULL);
   assert(SCIPisPositive(scip, alpha3));
   assert(SCIPisGE(scip, SCIPconsIsLocal(cons) ? SCIPvarGetLbLocal(x3) : SCIPvarGetLbGlobal(x3), -offset3));
   assert(basename != NULL);
   assert(N >= 1);
   assert(naddconss != NULL);

   SCIPdebugMessage("Creating linear Ben-Tal Nemirovski outer-approximation for <%s>.\n", basename);

   SCIP_CALL( SCIPallocBufferArray(scip, &avars, N+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bvars, N+1) );

   /* create additional variables */
   for( i = 0; i <= N; ++i )
   {
      (void) SCIPsnprintf(varname, 255, "soc#%s_a%d", basename, i);
      SCIP_CALL( SCIPcreateVar(scip, &avars[i], varname, 0.0, SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, SCIPconsIsLocal(cons), TRUE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, avars[i]) );

      (void) SCIPsnprintf(varname, 255, "soc#%s_b%d", basename, i);
      SCIP_CALL( SCIPcreateVar(scip, &bvars[i], varname, 0.0, SCIPinfinity(scip), 0.0, 
            SCIP_VARTYPE_CONTINUOUS, SCIPconsIsLocal(cons), TRUE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, bvars[i]) );
   }

   /* create first linear constraints - split into two because of the absolute value */
   vars[0] = avars[0];
   vals[0] = 1.0;
   vars[1] = x1;
   vals[1] = -alpha1;

   (void) SCIPsnprintf(linname, 255, "soc#%s#a%d", basename, 0);
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, alpha1 * offset1, SCIPinfinity(scip),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
         SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   ++*naddconss;

   vars[0] = avars[0];
   vals[0] = 1.0;
   vars[1] = x1;
   vals[1] = alpha1;

   (void) SCIPsnprintf(linname, 255, "soc#%s#A%d", basename, 0);
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, -alpha1 * offset1, SCIPinfinity(scip),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
         SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   ++*naddconss;

   if( x2 != NULL )
   {
      vars[0] = bvars[0];
      vals[0] = 1.0;
      vars[1] = x2;
      vals[1] = -alpha2;

      (void) SCIPsnprintf(linname, 255, "soc#%s#b%d", basename, 0);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, alpha2 * offset2, SCIPinfinity(scip),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++*naddconss;

      vars[0] = bvars[0];
      vals[0] = 1.0;
      vars[1] = x2;
      vals[1] = alpha2;

      (void) SCIPsnprintf(linname, 255, "soc#%s#B%d", basename, 0);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, -alpha2 * offset2, SCIPinfinity(scip),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++*naddconss;
   }
   else
   { /* second summand is just a constant */
      if( SCIPconsIsLocal(cons) )
      {
         SCIP_CALL( SCIPchgVarLbNode(scip, NULL, bvars[0], ABS(alpha2 * offset2)) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarLbGlobal(scip, bvars[0], ABS(alpha2 * offset2)) );
      }
   }

   /* create intermediate linear constraints */
   for( i = 1; i <= N; ++i )
   {
      SCIP_Real val;

      val = M_PI / pow(2.0, (double) (i+1));

      vars[0] = avars[i-1];
      vals[0] = cos(val);
      vars[1] = bvars[i-1];
      vals[1] = sin(val);
      vars[2] = avars[i];
      vals[2] = -1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#a%d", basename, i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, 0.0, 0.0,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++*naddconss;

      vars[0] = avars[i-1];
      vals[0] = sin(val);
      vars[1] = bvars[i-1];
      vals[1] = -cos(val);
      vars[2] = bvars[i];
      vals[2] = 1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#b%d", basename, i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, 0.0, SCIPinfinity(scip),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++*naddconss;

      vars[0] = avars[i-1];
      vals[0] = -sin(val);
      vars[1] = bvars[i-1];
      vals[1] = cos(val);
      vars[2] = bvars[i];
      vals[2] = 1.0;

      (void) SCIPsnprintf(linname, 255, "soc#%s#B%d", basename, i);
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 3, vars, vals, 0.0, SCIPinfinity(scip),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
            SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
            SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
            TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCons(scip, lincons) );
      SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
      ++*naddconss;
   }

   /* create last linear constraints */
   vars[0] = x3;
   vals[0] = alpha3;
   vars[1] = avars[N];
   vals[1] = -1.0;

   (void) SCIPsnprintf(linname, 255, "soc#%s#a%d", basename, N);
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, -alpha3 * offset3, SCIPinfinity(scip),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
         SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   ++*naddconss;

   vars[0] = avars[N];
   vals[0] = tan( M_PI / pow(2.0, (double) (N+1)) );
   vars[1] = bvars[N];
   vals[1] = -1.0;

   (void) SCIPsnprintf(linname, 255, "soc#%s#b%d", basename, i);
   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, linname, 2, vars, vals, 0.0, SCIPinfinity(scip),
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons),
         SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
         SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons),
         TRUE /* removable */, SCIPconsIsStickingAtNode(cons)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
   ++*naddconss;

   for( i = 0; i <= N; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &avars[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &bvars[i]) );
   }
   SCIPfreeBufferArray(scip, &avars);
   SCIPfreeBufferArray(scip, &bvars);

   return SCIP_OKAY;
}

/** adds a linear outer approx for a three dimensional SOC constraint
 * 
 * chooses between Ben-Tan/Nemirovski and Glineur and calls the corresponding function
 */
static
SCIP_RETCODE presolveCreateOuterApproxDim3(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< original constraint */
   SCIP_VAR*             x1,                 /**< variable x1 */
   SCIP_VAR*             x2,                 /**< variable x2 */
   SCIP_VAR*             x3,                 /**< variable x3 */
   SCIP_Real             alpha1,             /**< coefficient of x1 */
   SCIP_Real             alpha2,             /**< coefficient of x2 */
   SCIP_Real             alpha3,             /**< coefficient of x3 */
   SCIP_Real             offset1,            /**< offset of x1 */
   SCIP_Real             offset2,            /**< offset of x2 */
   SCIP_Real             offset3,            /**< offset of x3 */
   int                   N,                  /**< size of linear approximation, need to be >= 1 */
   SCIP_Bool             glineur,            /**< whether to prefer Glineur to Ben-Tal Nemirovski */
   const char*           basename,           /**< string to use for building variable and constraint names */
   int*                  naddconss           /**< buffer where to add the number of added constraints */
   )
{
   if( glineur )
   {
      SCIP_CALL( presolveCreateGlineurApproxDim3(scip, cons, x1, x2, x3, alpha1, alpha2, alpha3, offset1, offset2, offset3, N, basename, naddconss) );
   }
   else
   {
      SCIP_CALL( presolveCreateBenTalNemirovskiApproxDim3(scip, cons, x1, x2, x3, alpha1, alpha2, alpha3, offset1, offset2, offset3, N, basename, naddconss) );
   }

   return SCIP_OKAY;
}

/** adds linear outer approximation of Ben-Tal and Nemirovski for a constraint \f$\gamma + \sum_{i=1}^n (\alpha_i (x_i + \beta_i))^2 <= (\alpha_{n+1} (x_{n+1} + \beta_{n+1}))^2\f$ to the LP
 * 
 * if n>2, calls same function recursively;
 * if n=2, calls presolveCreateBenTalNemirovskiApproxDim3
 */
static
SCIP_RETCODE presolveCreateOuterApprox(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nlhsvars,           /**< number of variables on left hand side (n) */
   SCIP_VAR**            lhsvars,            /**< variables on left hand side (x_i) */
   SCIP_Real*            lhscoefs,           /**< coefficients of variables on left hand side (alpha_i) */
   SCIP_Real*            lhsoffsets,         /**< offsets of variable on left hand side (beta_i) */
   SCIP_VAR*             rhsvar,             /**< variable on right hand side (y) */
   SCIP_Real             rhscoeff,           /**< coefficient of variable on right hand side (alpha_{n+1}) */
   SCIP_Real             rhsoffset,          /**< offset of variable on right hand side (beta_{n+1}) */
   SCIP_Real             constant,           /**< constant term (gamma) */
   const char*           basename,           /**< prefix for variable and constraint name */
   SCIP_CONS*            origcons,           /**< original constraint for which this SOC3 set is added */
   int                   soc3_nr_auxvars,    /**< number of auxiliary variables to use for a SOC3 constraint, or 0 if automatic */
   SCIP_Bool             glineur,            /**< whether Glineur should be preferred to Ben-Tal Nemirovski */
   int*                  naddconss           /**< buffer where to add the number of added constraints */
   )
{
   char       name[255];
   SCIP_VAR*  auxvar1;
   SCIP_VAR*  auxvar2;

   assert(scip     != NULL);
   assert(lhsvars  != NULL);
   assert(nlhsvars >= 2);
   assert(lhscoefs != NULL);
   assert(lhsoffsets != NULL);
   assert(rhsvar   != NULL);
   assert(basename != NULL);
   assert(!SCIPisNegative(scip, constant));
   assert(naddconss != NULL);

   if( nlhsvars == 1 )
   { /* end of recursion */
      assert(SCIPisPositive(scip, constant));
      SCIP_CALL( presolveCreateOuterApproxDim3(scip, origcons,
            lhsvars[0],    NULL,           rhsvar,
            lhscoefs[0],   1.0,            rhscoeff,
            lhsoffsets[0], sqrt(constant), rhsoffset,
            soc3_nr_auxvars, glineur, basename, naddconss) );

      return SCIP_OKAY;
   }

   if( nlhsvars == 2 && SCIPisZero(scip, constant) )
   { /* end of recursion */
      assert(lhsvars[0] != NULL);
      assert(lhsvars[1] != NULL);
      assert(rhsvar     != NULL);
      SCIP_CALL( presolveCreateOuterApproxDim3(scip, origcons,
            lhsvars[0],    lhsvars[1],    rhsvar,
            lhscoefs[0],   lhscoefs[1],   rhscoeff,
            lhsoffsets[0], lhsoffsets[1], rhsoffset,
            soc3_nr_auxvars, glineur, basename, naddconss) );

      return SCIP_OKAY;
   }

   if( nlhsvars == 3 || (nlhsvars == 2 && !SCIPisZero(scip, constant)) )
   { 
      /* a bit special case too */
      /* for first two variables on lhs, create a new aux.var and a new SOC3 */
      (void) SCIPsnprintf(name, 255, "%s#z1", basename);
      SCIP_CALL( SCIPcreateVar(scip, &auxvar1, name, 0.0, SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar1) );

      /* constraint alpha_0 (x_0+beta0)^2 + alpha_1 (x_1+beta1)^2 <= auxvar^2 */
      SCIP_CALL( presolveCreateOuterApproxDim3(scip, origcons,
            lhsvars[0],    lhsvars[1],    auxvar1,
            lhscoefs[0],   lhscoefs[1],   1.0,
            lhsoffsets[0], lhsoffsets[1], 0.0,
            soc3_nr_auxvars, glineur, name, naddconss) );

      (void) SCIPsnprintf(name, 255, "%s_soc3", basename);
      if( nlhsvars == 3 )
      { /* create new constraint alpha_2 (x_2+beta2)^2 + auxvar^2 <= (rhscoeff * (rhsvar+rhsoffset))^2 */
         SCIP_CALL( presolveCreateOuterApproxDim3(scip, origcons,
               lhsvars[2],    auxvar1, rhsvar,
               lhscoefs[2],   1.0,     rhscoeff,
               lhsoffsets[2], 0.0,     rhsoffset,
               soc3_nr_auxvars, glineur, name, naddconss) );
      }
      else
      { /* create new constraint auxvar^2 + sqrt(constant)^2 <= (rhscoeff * (rhsvar+rhsoffset))^2 */
         SCIP_CALL( presolveCreateOuterApproxDim3(scip, origcons,
               auxvar1, NULL,           rhsvar,
               1.0,     1.0,            rhscoeff,
               0.0,     sqrt(constant), rhsoffset,
               soc3_nr_auxvars, glineur, name, naddconss) );
      }

      SCIP_CALL( SCIPreleaseVar(scip, &auxvar1) );

      return SCIP_OKAY;
   }

   /* nlhsvars >= 4 */

   (void) SCIPsnprintf(name, 255, "%s#z1", basename);
   SCIP_CALL( SCIPcreateVar(scip, &auxvar1, name, 0.0, SCIPinfinity(scip), 0.0,
         SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, auxvar1) );

   /* approx for left half of lhs */
   SCIP_CALL( presolveCreateOuterApprox(scip,
         nlhsvars/2, lhsvars, lhscoefs, lhsoffsets,
         auxvar1, 1.0, 0.0,
         constant, name, origcons, soc3_nr_auxvars, glineur, naddconss) );

   (void) SCIPsnprintf(name, 255, "%s#z2", basename);
   SCIP_CALL( SCIPcreateVar(scip, &auxvar2, name, 0., SCIPinfinity(scip), 0.0, 
         SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, auxvar2) );

   /* approx for right half of lhs */
   SCIP_CALL( presolveCreateOuterApprox(scip,
         nlhsvars-nlhsvars/2, &lhsvars[nlhsvars/2], &lhscoefs[nlhsvars/2], &lhsoffsets[nlhsvars/2],
         auxvar2, 1.0, 0.0,
         0.0, name, origcons, soc3_nr_auxvars, glineur, naddconss) );

   /* SOC constraint binding both auxvar's */
   (void)SCIPsnprintf(name, 255, "%s_soc3", basename);
   SCIP_CALL( presolveCreateOuterApproxDim3(scip, origcons,
         auxvar1, auxvar2, rhsvar,
         1.0,     1.0,     rhscoeff,
         0.0,     0.0,     rhsoffset,
         soc3_nr_auxvars, glineur, name, naddconss) );

   SCIP_CALL( SCIPreleaseVar(scip, &auxvar1) );
   SCIP_CALL( SCIPreleaseVar(scip, &auxvar2) );

   return SCIP_OKAY;
}

/** propagates variable bounds */
static
SCIP_RETCODE propagateBounds(
   SCIP*           scip,      /**< SCIP data structure */
   SCIP_CONS*      cons,      /**< constraint */
   SCIP_RESULT*    result,    /**< buffer to store result of propagation */
   int*            nchgbds    /**< buffer where to add number of tightened bounds */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_INTERVAL  lhsrange;
   SCIP_INTERVAL* lhsranges;
   SCIP_INTERVAL  rhsrange;
   SCIP_INTERVAL  a, b, c;
   SCIP_ROUNDMODE roundmode;
   SCIP_Bool      infeas, tightened;
   int            i;
   SCIP_Real      lb, ub;

   assert(scip   != NULL);
   assert(cons   != NULL);
   assert(result != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->ispropagated )
   {
      SCIPdebugMessage("skip propagation for constraint %s\n", SCIPconsGetName(cons));
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
   else
   {
      SCIPdebugMessage("try propagation for constraint %s\n", SCIPconsGetName(cons));
   }

   *result = SCIP_DIDNOTFIND;
   consdata->ispropagated = TRUE;

   /* @todo do something clever to decide whether propagation should be tried */

   SCIPintervalSetBounds(&lhsrange, consdata->constant - SCIPepsilon(scip), consdata->constant + SCIPepsilon(scip));

   SCIP_CALL( SCIPallocBufferArray(scip, &lhsranges, consdata->nvars) );
   for( i = 0; i < consdata->nvars; ++i )
   {
      lb = SCIPcomputeVarLbLocal(scip, consdata->vars[i]) - SCIPepsilon(scip);
      ub = SCIPcomputeVarUbLocal(scip, consdata->vars[i]) + SCIPepsilon(scip);
      SCIPintervalSetBounds(&lhsranges[i], MIN(lb, ub), MAX(lb, ub));
      if( consdata->offsets[i] != 0.0 )
         SCIPintervalAddScalar(SCIPinfinity(scip), &lhsranges[i], lhsranges[i], consdata->offsets[i]);
      if( consdata->coefs[i]   != 1.0 )
         SCIPintervalMulScalar(SCIPinfinity(scip), &lhsranges[i], lhsranges[i], consdata->coefs[i]);
      SCIPintervalSquare(SCIPinfinity(scip), &lhsranges[i], lhsranges[i]);

      SCIPintervalAdd(SCIPinfinity(scip), &lhsrange, lhsrange, lhsranges[i]);
   }

   if( SCIPvarGetStatus(consdata->rhsvar) != SCIP_VARSTATUS_MULTAGGR )
   {
      SCIPintervalSquareRoot(SCIPinfinity(scip), &a, lhsrange);
      if( consdata->rhscoeff != 1.0 )
         SCIPintervalDivScalar(SCIPinfinity(scip), &a, a, consdata->rhscoeff);
      if( consdata->rhsoffset != 0.0 )
         SCIPintervalSubScalar(SCIPinfinity(scip), &a, a, consdata->rhsoffset);
      SCIP_CALL( SCIPtightenVarLb(scip, consdata->rhsvar, SCIPintervalGetInf(a), FALSE, &infeas, &tightened) );
      if( infeas )
      {
         SCIPdebugMessage("propagation found constraint <%s> infeasible\n", SCIPconsGetName(cons));
         *result = SCIP_CUTOFF;
      }
      else if( tightened )
      {
         SCIPdebugMessage("propagation tightened bounds of rhs variable <%s> in constraint <%s>\n", SCIPvarGetName(consdata->rhsvar), SCIPconsGetName(cons));
         *result = SCIP_REDUCEDDOM;
         ++*nchgbds;
      }
   }

   if( *result != SCIP_CUTOFF )
   {
      lb = SCIPcomputeVarLbLocal(scip, consdata->rhsvar) - SCIPepsilon(scip);
      ub = SCIPcomputeVarUbLocal(scip, consdata->rhsvar) + SCIPepsilon(scip);
      SCIPintervalSetBounds(&rhsrange, MIN(lb, ub), MAX(lb, ub));
      if( consdata->rhsoffset != 0.0 )
         SCIPintervalAddScalar(SCIPinfinity(scip), &rhsrange, rhsrange, consdata->rhsoffset);
      if( consdata->rhscoeff  != 1.0 )
         SCIPintervalMulScalar(SCIPinfinity(scip), &rhsrange, rhsrange, consdata->rhscoeff);
      SCIPintervalSquare(SCIPinfinity(scip), &rhsrange, rhsrange);
      /* rhsrange = sqr(rhscoeff * (rhsvar + rhsoffset) ) */

      if( lhsrange.inf > rhsrange.sup )
      {
         SCIPdebugMessage("propagation found constraint <%s> infeasible: lhs = [%.15g,%.15g] > rhs = [%.15g,%.15g]\n",
            SCIPconsGetName(cons), lhsrange.inf, lhsrange.sup, rhsrange.inf, rhsrange.sup);
         *result = SCIP_CUTOFF;
      }
   }

   if( *result != SCIP_CUTOFF )
   {
      SCIPintervalSub(SCIPinfinity(scip), &b, rhsrange, lhsrange);  /*lint !e644 */
      for( i = 0; i < consdata->nvars; ++i )
      {
         if( SCIPvarGetStatus(consdata->vars[i]) == SCIP_VARSTATUS_MULTAGGR )
            continue;

         roundmode = SCIPintervalGetRoundingMode();
         if( !SCIPisInfinity(scip, b.sup) )
         {
            SCIPintervalSetRoundingModeUpwards();
            a.sup = b.sup + lhsranges[i].inf;
         }
         else
         {
            a.sup = SCIPinfinity(scip);
         }
         if( !SCIPisInfinity(scip, -b.inf) )
         {
            SCIPintervalSetRoundingModeDownwards();
            a.inf = b.inf + lhsranges[i].sup;
         }
         else
         {
            a.inf = -SCIPinfinity(scip);
         }
         SCIPintervalSetRoundingMode(roundmode);
         SCIPintervalSquareRoot(SCIPinfinity(scip), &a, a);

         assert(consdata->coefs[i] >= 0.0); /* should be ensured in create and presolveRemoveFixed */

         c = a;
         if( consdata->coefs[i]   != 1.0 )
            SCIPintervalDivScalar(SCIPinfinity(scip), &c, c, consdata->coefs[i]);
         if( consdata->offsets[i] != 0.0 )
            SCIPintervalSubScalar(SCIPinfinity(scip), &c, c, consdata->offsets[i]);

         SCIP_CALL( SCIPtightenVarUb(scip, consdata->vars[i], SCIPintervalGetSup(c), FALSE, &infeas, &tightened) );
         if( infeas )
         {
            SCIPdebugMessage("propagation found constraint <%s> infeasible\n", SCIPconsGetName(cons));
            *result = SCIP_CUTOFF;
            break;
         }
         else if( tightened )
         {
            SCIPdebugMessage("propagation tightened bounds of lhs variable <%s> in constraint <%s>\n", SCIPvarGetName(consdata->vars[i]), SCIPconsGetName(cons));
            *result = SCIP_REDUCEDDOM;
            ++*nchgbds;
         }

         c = a;
         SCIPintervalDivScalar(SCIPinfinity(scip), &c, c, -consdata->coefs[i]);
         if( consdata->offsets[i] != 0.0 )
            SCIPintervalSubScalar(SCIPinfinity(scip), &c, c, consdata->offsets[i]);

         SCIP_CALL( SCIPtightenVarLb(scip, consdata->vars[i], SCIPintervalGetInf(c), FALSE, &infeas, &tightened) );
         if( infeas )
         {
            SCIPdebugMessage("propagation found constraint <%s> infeasible\n", SCIPconsGetName(cons));
            *result = SCIP_CUTOFF;
            break;
         }
         else if( tightened )
         {
            SCIPdebugMessage("propagation tightened bounds of lhs variable <%s> in constraint <%s>\n", SCIPvarGetName(consdata->vars[i]), SCIPconsGetName(cons));
            *result = SCIP_REDUCEDDOM;
            ++*nchgbds;
         }
      }
   }

   SCIPfreeBufferArray(scip, &lhsranges);

   if( *result != SCIP_DIDNOTFIND )
   {
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
   }

   return SCIP_OKAY;
}

/** tries to adjust a solution such that it satisfies a given constraint by increasing the value for the constraints right hand side variable */
static
SCIP_RETCODE polishSolution(
   SCIP*           scip,      /**< SCIP data structure */
   SCIP_CONS*      cons,      /**< constraint */
   SCIP_SOL*       sol,       /**< solution to polish */
   SCIP_Bool*      success    /**< buffer to store whether polishing was successful */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real rhsval;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(sol  != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, consdata->rhscoeff));

   /* compute minimal rhs variable value so that constraint is satisfied */
   if( !SCIPisInfinity(scip, consdata->lhsval) )
      rhsval = consdata->lhsval / consdata->rhscoeff - consdata->rhsoffset;
   else
      rhsval = consdata->rhscoeff > 0.0 ? SCIPinfinity(scip) : -SCIPinfinity(scip);

   if( consdata->rhscoeff > 0.0 )
   {
      assert(SCIPvarMayRoundUp(consdata->rhsvar));

      /* round rhsval up, if variable is integral */
      if( SCIPvarIsIntegral(consdata->rhsvar) && !SCIPisInfinity(scip, rhsval) )
         rhsval = SCIPceil(scip, rhsval);

      /* if new value is above upper bound, we are lost */
      if( SCIPisGT(scip, rhsval, SCIPvarGetUbGlobal(consdata->rhsvar)) )
      {
         *success = FALSE;
      }
      else
      {
         /* if new value is larger then current one, increase to new value */
         if( rhsval > SCIPgetSolVal(scip, sol, consdata->rhsvar) )
         {
            SCIPdebugMessage("increase <%s> to %g\n", SCIPvarGetName(consdata->rhsvar), rhsval);
            SCIP_CALL( SCIPsetSolVal(scip, sol, consdata->rhsvar, rhsval) );
         }

         *success = TRUE;
      }
   }
   else
   {
      assert(SCIPvarMayRoundDown(consdata->rhsvar));

      /* round rhsval down, if variable is integral */
      if( SCIPvarIsIntegral(consdata->rhsvar) )
         rhsval = SCIPfloor(scip, rhsval);

      /* if new value is below lower bound, we are lost */
      if( SCIPisLT(scip, rhsval, SCIPvarGetLbGlobal(consdata->rhsvar)) )
      {
         *success = FALSE;
      }
      else
      {
         /* if new value is below current one, decrease to new value */
         if( rhsval < SCIPgetSolVal(scip, sol, consdata->rhsvar) )
         {
            SCIPdebugMessage("decrease <%s> to %g\n", SCIPvarGetName(consdata->rhsvar), rhsval);
            SCIP_CALL( SCIPsetSolVal(scip, sol, consdata->rhsvar, rhsval) );
         }

         *success = TRUE;
      }
   }

   SCIPdebugMessage("polishing solution for constraint <%s> was %ssuccessful\n", SCIPconsGetName(cons), *success ? "" : "not ");

   return SCIP_OKAY;
}

/*
 * Quadratic constraint upgrading
 */


#ifdef QUADCONSUPGD_PRIORITY
/** tries to upgrade a quadratic constraint to a SOC constraint
 * @todo more general quadratic constraints then sums of squares might allow an upgrade to a SOC
 */
static
SCIP_DECL_QUADCONSUPGD(upgradeConsQuadratic)
{
   int            nquadvars;
   SCIP_QUADVARTERM* term;
   SCIP_VAR**     lhsvars;
   SCIP_Real*     lhscoefs;
   SCIP_Real*     lhsoffsets;
   SCIP_Real      lhsconstant;
   int            lhscount;
   SCIP_VAR*      rhsvar; 
   SCIP_Real      rhscoef;
   SCIP_Real      rhsoffset;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nupgdconss != NULL);
   assert(upgdconss  != NULL);

   *nupgdconss = 0;

   SCIPdebugMessage("upgradeConsQuadratic called for constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );

   /* currently do not support linear parts in upgrading of SOC constraints */
   if( SCIPgetNLinearVarsQuadratic(scip, cons) )
      return SCIP_OKAY;

   /* currently do not support bilinear parts in upgrading of SOC constraints */
   if( SCIPgetNBilinTermsQuadratic(scip, cons) )
      return SCIP_OKAY;

   nquadvars = SCIPgetNQuadVarTermsQuadratic(scip, cons);

   /* currently, a proper SOC constraint needs at least 3 variables */
   if( nquadvars < 3 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &lhsvars,    nquadvars-1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhscoefs,   nquadvars-1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lhsoffsets, nquadvars-1) );

   lhsconstant = 0.0;
   lhscount = 0;
   rhsvar = NULL; 
   rhscoef = 0.0;
   rhsoffset = 0.0;

   if( !SCIPisInfinity(scip, SCIPgetRhsQuadratic(scip, cons)) )
   { /* try whether constraint on right hand side is SOC */
      lhsconstant = -SCIPgetRhsQuadratic(scip, cons);

      for( i = 0; i < nquadvars; ++i )
      {
         term = &SCIPgetQuadVarTermsQuadratic(scip, cons)[i];

         /* if there is a linear variable that is still considered as quadratic (constraint probably not presolved yet), then give up */
         if( term->sqrcoef == 0.0 )
            goto cleanup;

         if( term->sqrcoef > 0.0 )
         {
            if( lhscount >= nquadvars - 1 )
            { /* too many variables on lhs, i.e., all variables seem to have positive coefficient */
               rhsvar = NULL;
               break;
            }

            lhsvars[lhscount]  = term->var;
            lhscoefs[lhscount] = sqrt(term->sqrcoef);

            if( term->lincoef != 0.0 )
            {
               lhsoffsets[lhscount] = term->lincoef / (2 * term->sqrcoef);
               lhsconstant -= term->lincoef * term->lincoef / (4 * term->sqrcoef);
            }
            else
            {
               lhsoffsets[lhscount] = 0.0;
            }

            ++lhscount;
         }
         else if( rhsvar != NULL || SCIPisLT(scip, SCIPcomputeVarLbGlobal(scip, term->var), term->lincoef / (2*term->sqrcoef)) )
         { /* second variable with negative coefficient -> cannot be SOC */
            /* or lower bound of variable does not ensure positivity of right hand side */
            rhsvar = NULL;
            break;
         }
         else
         {
            rhsvar       = term->var;
            rhscoef      = sqrt(-term->sqrcoef);
            rhsoffset    = term->lincoef / (2 * term->sqrcoef);
            lhsconstant -= term->lincoef * term->lincoef / (4 * term->sqrcoef);
         }
      }
   }

   if( rhsvar != NULL && lhscount >= 2 && !SCIPisNegative(scip, lhsconstant) )
   { /* found SOC constraint, so upgrade to SOC constraint(s) (below) and relax right hand side */
      SCIPdebugMessage("found right hand side of constraint <%s> to be SOC\n", SCIPconsGetName(cons));

      /* check if upgdconss is long enough to store upgrade constraints:
       * we need two if we will have a quadratic constraint for the left hand side left */
      *nupgdconss = SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, cons)) ? 1 : 2;
      if( *nupgdconss > upgdconsssize )
      {
         /* signal that we need more memory and return */
         *nupgdconss = -*nupgdconss;
         goto cleanup;
      }

      SCIP_CALL( SCIPcreateConsSOC(scip, &upgdconss[0], SCIPconsGetName(cons),
            lhscount, lhsvars, lhscoefs, lhsoffsets, MAX(lhsconstant, 0.0),
            rhsvar, rhscoef, rhsoffset,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, upgdconss[0], NULL) ) );

      /* create constraint that is equal to cons except that rhs is now infinity */
      if( !SCIPisInfinity(scip, -SCIPgetLhsQuadratic(scip, cons)) )
      {
         SCIP_CALL( SCIPcreateConsQuadratic2(scip, &upgdconss[1], SCIPconsGetName(cons),
               SCIPgetNLinearVarsQuadratic(scip, cons), SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
               SCIPgetNQuadVarTermsQuadratic(scip, cons), SCIPgetQuadVarTermsQuadratic(scip, cons),
               SCIPgetNBilinTermsQuadratic(scip, cons), SCIPgetBilinTermsQuadratic(scip, cons),
               SCIPgetLhsQuadratic(scip, cons), SCIPinfinity(scip),
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
      }
   }
   else if( !SCIPisInfinity(scip, - SCIPgetLhsQuadratic(scip, cons)) )
   { /* if the first failed, try if constraint on left hand side is SOC (using negated coefficients) */
      lhscount = 0;
      rhsvar = NULL;

      lhsconstant = SCIPgetLhsQuadratic(scip, cons);

      for( i = 0; i < nquadvars; ++i )
      {
         term = &SCIPgetQuadVarTermsQuadratic(scip, cons)[i];

         /* if there is a linear variable that is still considered as quadratic (constraint probably not presolved yet), then give up */
         if( term->sqrcoef == 0.0 )
            goto cleanup;

         if( term->sqrcoef < 0.0 )
         {
            if( lhscount >= nquadvars )
            { /* too many variables on lhs, i.e., all variables seem to have negative coefficient */
               rhsvar = NULL;
               break;
            }

            lhsvars[lhscount]  = term->var;
            lhscoefs[lhscount] = sqrt(-term->sqrcoef);

            if( term->lincoef != 0.0 )
            {
               lhsoffsets[lhscount] = term->lincoef / (2 * term->sqrcoef);
               lhsconstant += term->lincoef * term->lincoef / (4 * term->sqrcoef);
            }
            else
            {
               lhsoffsets[lhscount] = 0.0;
            }

            ++lhscount;
         }
         else if( rhsvar || SCIPisLT(scip, SCIPcomputeVarLbGlobal(scip, term->var), -term->lincoef / (2*term->sqrcoef)) )
         { /* second variable with positive coefficient -> cannot be SOC */
            /* or lower bound of variable does not ensure positivity of right hand side */
            rhsvar = NULL;
            break;
         }
         else
         {
            rhsvar       = term->var;
            rhscoef      = sqrt(term->sqrcoef);
            rhsoffset    = term->lincoef / (2 * term->sqrcoef);
            lhsconstant += term->lincoef * term->lincoef / (4 * term->sqrcoef);
         }
      }

      if( rhsvar && lhscount >= 2 && !SCIPisNegative(scip, lhsconstant) )
      { /* found SOC constraint, so upgrade to SOC constraint(s) (below) and relax left hand side */
         SCIPdebugMessage("found left hand side of constraint <%s> to be SOC\n", SCIPconsGetName(cons));

         /* check if upgdconss is long enough to store upgrade constraints:
          * we need two if we will have a quadratic constraint for the right hand side left */
         *nupgdconss = SCIPisInfinity(scip, SCIPgetRhsQuadratic(scip, cons)) ? 1 : 2;
         if( *nupgdconss > upgdconsssize )
         {
            /* signal that we need more memory and return */
            *nupgdconss = -*nupgdconss;
            goto cleanup;
         }

         SCIP_CALL( SCIPcreateConsSOC(scip, &upgdconss[0], SCIPconsGetName(cons),
               lhscount, lhsvars, lhscoefs, lhsoffsets, MAX(lhsconstant, 0.0),
               rhsvar, rhscoef, rhsoffset,
               SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
               SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
               SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, upgdconss[0], NULL) ) );

         /* create constraint that is equal to cons except that lhs is now -infinity */
         if( !SCIPisInfinity(scip, SCIPgetRhsQuadratic(scip, cons)) )
         {
            SCIP_CALL( SCIPcreateConsQuadratic2(scip, &upgdconss[1], SCIPconsGetName(cons),
                  SCIPgetNLinearVarsQuadratic(scip, cons), SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
                  SCIPgetNQuadVarTermsQuadratic(scip, cons), SCIPgetQuadVarTermsQuadratic(scip, cons),
                  SCIPgetNBilinTermsQuadratic(scip, cons), SCIPgetBilinTermsQuadratic(scip, cons),
                  -SCIPinfinity(scip), SCIPgetRhsQuadratic(scip, cons),
                  SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
                  SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
                  SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons)) );
         }
      }
   }

 cleanup:
   SCIPfreeBufferArray(scip, &lhsvars);
   SCIPfreeBufferArray(scip, &lhscoefs);
   SCIPfreeBufferArray(scip, &lhsoffsets);

   return SCIP_OKAY;
} /*lint !e715*/
#endif

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySOC)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrSOC(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeSOC)
{  
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitSOC)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur  = SCIPfindHeur(scip, "subnlp");
   conshdlrdata->trysolheur  = SCIPfindHeur(scip, "trysol");
   conshdlrdata->haveexprint = (strcmp(SCIPexprintGetName(), "NONE") != 0);

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitSOC)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur  = NULL;
   conshdlrdata->trysolheur  = NULL;
   conshdlrdata->haveexprint = FALSE;

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreSOC)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of soc constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreSOC NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreSOC)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool iscutoff;
   SCIP_Bool isdeleted;
   int dummy;
   int c;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   if( nconss == 0 )
      return SCIP_OKAY;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( presolveRemoveFixedVariables(scip, conshdlr, conss[c], &dummy, &dummy, &dummy, &dummy, &iscutoff, &isdeleted) ); /*lint !e613*/
      if( iscutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* if conss[c] has been deleted, skip the rest */
      if( isdeleted )
         continue;

      /* tell SCIP that we have something nonlinear */
      SCIPmarkNonlinearitiesPresent(scip);
      if( !SCIPhasContinuousNonlinearitiesPresent(scip) )
      {
         consdata = SCIPconsGetData(conss[c]); /*lint !e613*/
         assert(consdata != NULL);

         for( i = 0; i < consdata->nvars; ++i )
            if( SCIPvarGetType(consdata->vars[i]) >= SCIP_VARTYPE_CONTINUOUS )
            {
               SCIPmarkContinuousNonlinearitiesPresent(scip);
               break;
            }
      }
   }

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolSOC)
{  
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr);

   /* add nlrow representation to NLP, if NLP has been enabled */
   if( SCIPisNLPConstructed(scip) )
   {
      for( c = 0; c < nconss; ++c )
      {
         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);

         if( consdata->nlrow == NULL )
         {
            SCIP_CALL( createNlRow(scip, conshdlr, conss[c]) );
            assert(consdata->nlrow != NULL);
         }
         SCIP_CALL( SCIPaddNlRow(scip, consdata->nlrow) );
      }
   }

   conshdlrdata->newsoleventfilterpos = -1;
   if( nconss != 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, &conshdlrdata->newsoleventfilterpos) );
   }

   /* reset flags and counters */
   conshdlrdata->sepanlp = FALSE;
   conshdlrdata->lastenfolpnode = NULL;
   conshdlrdata->nenfolprounds = 0;

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolSOC)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip     != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr);

   if( conshdlrdata->newsoleventfilterpos >= 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, conshdlrdata->newsoleventfilterpos) );
      conshdlrdata->newsoleventfilterpos = -1;
   }

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* free nonlinear row representation */
      if( consdata->nlrow != NULL )
      {
         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }
   }

   return SCIP_OKAY;
} /*lint !e715*/


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSOC)
{
   int i;

   assert(scip      != NULL);
   assert(conshdlr  != NULL);
   assert(cons      != NULL);
   assert(consdata  != NULL);
   assert(*consdata != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert((*consdata)->nlrow == NULL); /* should have been freed in exitsol */

   SCIPdebugMessage("Deleting SOC constraint <%s>.\n", SCIPconsGetName(cons) );

   if( SCIPconsIsTransformed(cons) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   for( i = 0; i < (*consdata)->nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->vars[i]) );
   }

   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars,    (*consdata)->nvars);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->coefs,   (*consdata)->nvars);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->offsets, (*consdata)->nvars);
   assert((*consdata)->lhsbndchgeventdatas == NULL);

   if( (*consdata)->rhsvar != NULL )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->rhsvar) );
   }

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransSOC)
{  
   SCIP_CONSDATA*     consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     sourcedata;
   char               s[SCIP_MAXSTRLEN];
   int                i;

   assert(scip       != NULL);
   assert(conshdlr   != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMessage("Transforming SOC constraint: <%s>.\n", SCIPconsGetName(sourcecons) );

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata       != NULL);
   assert(sourcedata->vars != NULL);
   assert(sourcedata->coefs != NULL);
   assert(sourcedata->offsets != NULL);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->nvars = sourcedata->nvars;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, consdata->nvars) );
   SCIP_CALL( SCIPgetTransformedVars(scip, consdata->nvars, sourcedata->vars, consdata->vars) );
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[i]) );
   }

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->coefs,   sourcedata->coefs,   consdata->nvars) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->offsets, sourcedata->offsets, consdata->nvars) );
   consdata->constant = sourcedata->constant;

   SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->rhsvar, &consdata->rhsvar) );
   consdata->rhsoffset = sourcedata->rhsoffset;
   consdata->rhscoeff  = sourcedata->rhscoeff;

   SCIP_CALL( SCIPcaptureVar(scip, consdata->rhsvar) );

   consdata->nlrow = NULL;
   consdata->lhsbndchgeventdatas = NULL;
   consdata->ispropagated  = FALSE;
   consdata->isapproxadded = FALSE;

   /* create transformed constraint with the same flags */
   (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));
   SCIP_CALL( SCIPcreateCons(scip, targetcons, s, conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, *targetcons) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpSOC)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of soc constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpSOC NULL
#endif


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpSOC)
{  
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          sepasuccess;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, conshdlrdata->doscaling, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   /* at root, check if we want to solve the NLP relaxation and use its solutions as reference point
    * if there is something convex, then linearizing in the solution of the NLP relaxation can be very useful
    */
   if( SCIPgetDepth(scip) == 0 && !conshdlrdata->sepanlp &&
      (SCIPgetNContVars(scip) >= conshdlrdata->sepanlpmincont * SCIPgetNVars(scip) || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY) &&
      SCIPisNLPConstructed(scip) && SCIPgetNNlpis(scip) > 0 )
   {
      SCIP_NLPSOLSTAT solstat;
      SCIP_Bool       solvednlp;

      solstat = SCIPgetNLPSolstat(scip);
      solvednlp = FALSE;
      if( solstat == SCIP_NLPSOLSTAT_UNKNOWN )
      {
         /* NLP is not solved yet, so we do it now and update solstat */

         /* ensure linear conss are in NLP */
         if( conshdlrdata->subnlpheur != NULL )
         {
            SCIP_CALL( SCIPaddLinearConsToNlpHeurSubNlp(scip, conshdlrdata->subnlpheur, TRUE, TRUE) );
         }

         /* set LP solution as starting values, if available */
         if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
         {
            SCIP_CALL( SCIPsetNLPInitialGuessSol(scip, NULL) );
         }

         /* SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_VERBLEVEL, 1) ); */
         SCIP_CALL( SCIPsolveNLP(scip) );

         solstat = SCIPgetNLPSolstat(scip);
         SCIPdebugMessage("solved NLP relax, solution status: %d\n", solstat);

         solvednlp = TRUE;
      }

      conshdlrdata->sepanlp = TRUE;

      if( solstat == SCIP_NLPSOLSTAT_GLOBINFEASIBLE )
      {
         SCIPdebugMessage("NLP relaxation is globally infeasible, thus can cutoff node\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if( solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         /* if we have feasible NLP solution, generate linearization cuts there */
         SCIP_Bool lpsolseparated;
         SCIP_SOL* nlpsol;

         SCIP_CALL( SCIPcreateNLPSol(scip, &nlpsol, NULL) );
         assert(nlpsol != NULL);

         /* if we solved the NLP and solution is integral, then pass it to trysol heuristic */
         if( solvednlp && conshdlrdata->trysolheur != NULL )
         {
            int nfracvars;

            nfracvars = 0;
            if( SCIPgetNBinVars(scip) > 0 || SCIPgetNIntVars(scip) > 0 )
            {
               SCIP_CALL( SCIPgetNLPFracVars(scip, NULL, NULL, NULL, &nfracvars, NULL) );
            }

            if( nfracvars == 0 )
            {
               SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, nlpsol) );
            }
         }

         SCIP_CALL( addLinearizationCuts(scip, conshdlr, conss, nconss, nlpsol, &lpsolseparated, conshdlrdata->minefficacy) );

         SCIP_CALL( SCIPfreeSol(scip, &nlpsol) );

         /* if a cut that separated the LP solution was added, then return, otherwise continue with usual separation in LP solution */
         if( lpsolseparated )
         {
            SCIPdebugMessage("linearization cuts separate LP solution\n");

            *result = SCIP_SEPARATED;

            return SCIP_OKAY;
         }
      }
   }
   /* if we do not want to try solving the NLP, or have no NLP, or have no NLP solver, or solving the NLP failed,
    * or separating with NLP solution as reference point failed, then try (again) with LP solution as reference point
    */

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, FALSE, &sepasuccess) );
   if( sepasuccess )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolSOC)
{  
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          sepasuccess;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);
   assert(sol      != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, sol, conshdlrdata->doscaling, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, FALSE, &sepasuccess) );
   if( sepasuccess )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpSOC)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_CONS*         maxviolcons;
   SCIP_Bool          success;
   int                nbndchg;
   int                c;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result   != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, conshdlrdata->doscaling, &maxviolcons) );

   if( maxviolcons == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   /* if we are above the 100'th enforcement round for this node, something is strange
    * (maybe the LP does not think that the cuts we add are violated, or we do ECP on a high-dimensional convex function)
    * in this case, check if some limit is hit or SCIP should stop for some other reason and terminate enforcement by creating a dummy node
    * (in optimized more, returning SCIP_INFEASIBLE in *result would be sufficient, but in debug mode this would give an assert in scip.c)
    * the reason to wait for 100 rounds is to avoid calls to SCIPisStopped in normal runs, which may be expensive
    * we only increment nenfolprounds until 101 to avoid an overflow
    */
   if( conshdlrdata->lastenfolpnode == SCIPgetCurrentNode(scip) )
   {
      if( conshdlrdata->nenfolprounds > 100 )
      {
         if( SCIPisStopped(scip) )
         {
            SCIP_NODE* child;

            SCIP_CALL( SCIPcreateChild(scip, &child, 1.0, SCIPnodeGetEstimate(SCIPgetCurrentNode(scip))) );
            *result = SCIP_BRANCHED;

            return SCIP_OKAY;
         }
      }
      else
         ++conshdlrdata->nenfolprounds;
   }
   else
   {
      conshdlrdata->lastenfolpnode = SCIPgetCurrentNode(scip);
      conshdlrdata->nenfolprounds = 0;
   }

   /* try separation, this should usually work */
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, TRUE, &success) );
   if( success )
   {
      SCIPdebugMessage("enforced by separation\n");
      *result = SCIP_SEPARATED;
      return SCIP_OKAY;
   }

   /* try propagation */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      if( !SCIPisGT(scip, consdata->violation, SCIPfeastol(scip)) )
         continue;

      nbndchg = 0;
      SCIP_CALL( propagateBounds(scip, conss[c], result, &nbndchg) );  /*lint !e613*/
      if( *result == SCIP_CUTOFF || *result == SCIP_REDUCEDDOM )
      {
         SCIPdebugMessage("enforced by %s\n", *result == SCIP_CUTOFF ? "cutting off node" : "reducing domain");
         return SCIP_OKAY;
      }
   }

   SCIPwarningMessage("could not enforce feasibility by separating or branching; declaring solution with viol %g feasible\n", SCIPconsGetData(maxviolcons)->violation);
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
} /*lint !e715*/


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSOC)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcons;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, conshdlrdata->doscaling, &maxviolcons) );

   if( maxviolcons == NULL )
      *result = SCIP_FEASIBLE;

   *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
} /*lint !e715*/


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckSOC)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          maxviol;
   SCIP_Bool          dolinfeasshift;
   SCIP_SOL*          polishedsol;
   int                c;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result     = SCIP_FEASIBLE;
   maxviol     = 0.0;

   dolinfeasshift = conshdlrdata->linfeasshift && (conshdlrdata->trysolheur != NULL);
   polishedsol = NULL;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( computeViolation(scip, conss[c], sol, conshdlrdata->doscaling) );  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* if feasible, just continue */
      if( !SCIPisGT(scip, consdata->violation, SCIPfeastol(scip)) )
         continue;

      *result = SCIP_INFEASIBLE;

      if( consdata->violation > maxviol )
         maxviol = consdata->violation;

      if( printreason )
      {
         SCIP_Real unscaledviol;

         unscaledviol  = consdata->lhsval;
         if( !SCIPisInfinity(scip, unscaledviol) )
            unscaledviol -= consdata->rhscoeff * (SCIPgetSolVal(scip, sol, consdata->rhsvar) + consdata->rhsoffset);

         SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );  /*lint !e613*/            
         SCIPinfoMessage(scip, NULL, "\tviolation: %g (scaled: %g)\n", unscaledviol, consdata->violation);
      }

      /* if we do linear feasibility shifting, then try to adjust solution */
      if( dolinfeasshift )
      {
         if( SCIPvarGetStatus(consdata->rhsvar) != SCIP_VARSTATUS_MULTAGGR &&
            (  (consdata->rhscoeff > 0.0 && SCIPvarMayRoundUp  (consdata->rhsvar)) ||
               (consdata->rhscoeff < 0.0 && SCIPvarMayRoundDown(consdata->rhsvar)) ) )
         {
            SCIP_Bool success;

            if( polishedsol == NULL )
            {
               if( sol != NULL )
               {
                  SCIP_CALL( SCIPcreateSolCopy(scip, &polishedsol, sol) );
               }
               else
               {
                  SCIP_CALL( SCIPcreateLPSol(scip, &polishedsol, NULL) );
               }
               SCIP_CALL( SCIPunlinkSol(scip, polishedsol) );
            }
            SCIP_CALL( polishSolution(scip, conss[c], polishedsol, &success) );  /*lint !e613*/

            /* disable solution polishing if we failed for this constraint */
            dolinfeasshift = success;
         }
         else /* if locks of the variable are bad or rhs is multi-aggregated, disable solution polishing */
         {
            dolinfeasshift = FALSE;
         }
      }

      /* if solution polishing is off and there is no NLP heuristic or we just check the LP solution,
       * then there is no need to check remaining constraints (NLP heuristic will pick up LP solution anyway) */
      if( !dolinfeasshift && (conshdlrdata->subnlpheur == NULL || sol == NULL))
         break;
   }

   /* if we failed to polish solution, clear the solution */ 
   if( !dolinfeasshift && polishedsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &polishedsol) );
   }

   if( polishedsol != NULL )
   {
      assert(*result == SCIP_INFEASIBLE);
      SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, polishedsol) );
      SCIP_CALL( SCIPfreeSol(scip, &polishedsol) );
   }
   else if( conshdlrdata->subnlpheur != NULL && sol != NULL && *result == SCIP_INFEASIBLE )
   {
      SCIP_CALL( SCIPupdateStartpointHeurSubNlp(scip, conshdlrdata->subnlpheur, sol, maxviol) );
   }

   return SCIP_OKAY;
} /*lint !e715*/


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSOC)
{
   SCIP_RESULT propresult;
   int         c;
   int         nchgbds;

   assert(scip     != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;
   nchgbds = 0;

   for( c = 0; c < nconss && *result != SCIP_CUTOFF; ++c )
   {
      SCIP_CALL( propagateBounds(scip, conss[c], &propresult, &nchgbds) );  /*lint !e613*/
      if( propresult != SCIP_DIDNOTFIND && propresult != SCIP_DIDNOTRUN )
         *result = propresult;
   }

   return SCIP_OKAY;
} /*lint !e715*/


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSOC)
{
   SCIP_CONSHDLRDATA*  conshdlrdata;
   SCIP_CONSDATA*      consdata;
   int                 c;
   SCIP_RESULT         propresult;
   SCIP_Bool           iscutoff;
   SCIP_Bool           isdeleted;

   assert(scip     != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(conshdlr != NULL);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      SCIP_CALL( presolveRemoveFixedVariables(scip, conshdlr, conss[c], ndelconss, nupgdconss, nchgbds, nfixedvars, &iscutoff, &isdeleted) );  /*lint !e613*/
      if( iscutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( isdeleted )
      { /* conss[c] has been deleted */
         *result = SCIP_SUCCESS;
         continue;
      }

      if( conshdlrdata->nauxvars > 0 && !consdata->isapproxadded )
      {
         SCIP_CALL( presolveCreateOuterApprox(scip, consdata->nvars, consdata->vars, consdata->coefs, consdata->offsets, consdata->rhsvar, consdata->rhscoeff, consdata->rhscoeff, consdata->constant, SCIPconsGetName(conss[c]), conss[c], conshdlrdata->nauxvars, conshdlrdata->glineur, naddconss) );  /*lint !e613*/
         consdata->isapproxadded = TRUE;
      }

      SCIP_CALL( propagateBounds(scip, conss[c], &propresult, nchgbds) );  /*lint !e613*/
      switch( propresult )
      {
      case SCIP_DIDNOTRUN:
      case SCIP_DIDNOTFIND:
         break;
      case SCIP_REDUCEDDOM:
         *result = SCIP_SUCCESS;
         break;
      case SCIP_CUTOFF:
         *result = SCIP_CUTOFF;
         SCIPdebugMessage("infeasible in presolve due to propagation for constraint %s\n", SCIPconsGetName(conss[c]));  /*lint !e613*/
         break;
      default:
         SCIPerrorMessage("unexpected result from propagation: %d\n", propresult);
         return SCIP_ERROR;
      } /*lint !e788*/
   }

   return SCIP_OKAY;
} /*lint !e715*/


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropSOC)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of soc constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropSOC NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockSOC)
{
   SCIP_CONSDATA* consdata;
   int            i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(cons     != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("Locking constraint <%s>.\n", SCIPconsGetName(cons));

   /* Changing variables x_i, i <= n, in both directions can lead to an infeasible solution. */
   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   /* Rounding x_{n+1} up will not violate a solution. */
   if( consdata->rhsvar != NULL )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->rhsvar, consdata->rhscoeff > 0.0 ? nlockspos : nlocksneg, consdata->rhscoeff > 0.0 ? nlocksneg : nlockspos) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveSOC)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of soc constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveSOC NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveSOC)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of soc constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveSOC NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableSOC)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of soc constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableSOC NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableSOC)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of soc constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableSOC NULL
#endif


/** variable deletion method of constraint handler */
#define consDelvarsSOC NULL


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintSOC)
{  
   SCIP_CONSDATA* consdata;
   int            i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(cons     != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPinfoMessage(scip, file, "sqrt( ");
   if( consdata->constant != 0.0 )
   {
      SCIPinfoMessage(scip, file, "%.15g", consdata->constant);
   }

   for( i = 0; i < consdata->nvars; ++i )
   {
      SCIPinfoMessage(scip, file, "+ (%.15g*(", consdata->coefs[i]);
      SCIP_CALL( SCIPwriteVarName(scip, file, consdata->vars[i], TRUE) );
      SCIPinfoMessage(scip, file, "%+.15g))^2 ", consdata->offsets[i]);
   }

   SCIPinfoMessage(scip, file, ") <= ");
   if( consdata->rhsvar != NULL )
   {
      SCIPinfoMessage(scip, file, "%.15g*(", consdata->rhscoeff);
      SCIP_CALL( SCIPwriteVarName(scip, file, consdata->rhsvar, TRUE) );
      SCIPinfoMessage(scip, file, "%+.15g)", consdata->rhsoffset);
   }
   else
   {
      SCIPinfoMessage(scip, file, "%.15g", consdata->rhscoeff*consdata->rhsoffset);
   }

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopySOC)
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR**     vars;
   SCIP_VAR*      rhsvar;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(sourcescip != NULL);
   assert(sourceconshdlr != NULL);
   assert(sourcecons != NULL);
   assert(varmap != NULL);
   assert(valid != NULL);
   assert(stickingatnode == FALSE);

   consdata = SCIPconsGetData(sourcecons);
   assert(consdata != NULL);

   *valid = TRUE; 

   SCIP_CALL( SCIPallocBufferArray(sourcescip, &vars, consdata->nvars) );

   /* map variables to active variables of the target SCIP */   
   for( i = 0; i < consdata->nvars && *valid; ++i )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, consdata->vars[i], &vars[i], varmap, consmap, global, valid) );
      assert(!(*valid) || vars[i] != NULL);
   }

   /* map rhs variable to active variable of the target SCIP */   
   if( *valid )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, consdata->rhsvar, &rhsvar, varmap, consmap, global, valid) );
      assert(!(*valid) || rhsvar != NULL);
   }

   /* only create the target constraint, if all variables could be copied */
   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsSOC(scip, cons, name ? name : SCIPconsGetName(sourcecons),
            consdata->nvars, vars, consdata->coefs, consdata->offsets, consdata->constant,
            rhsvar, consdata->rhscoeff, consdata->rhsoffset,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );  /*lint !e644 */
   }

   SCIPfreeBufferArray(sourcescip, &vars);

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
#if 1
static
SCIP_DECL_CONSPARSE(consParseSOC)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_Real* offsets;
   int nvars;
   int varssize;
   SCIP_VAR* rhsvar;
   SCIP_Real rhscoef;
   SCIP_Real rhsoffset;
   SCIP_Real constant;
   SCIP_Real coef;
   SCIP_Real offset;
   char* endptr;

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   /* check that string starts with "sqrt( " */
   if( strncmp(str, "sqrt( ", 6) != 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected 'sqrt( ' at begin of soc constraint string '%s'\n", str);
      *success = FALSE;
      return SCIP_OKAY;
   }
   str += 6;

   *success = TRUE;

   /* check if we have a constant in the beginning */
   if( SCIPstrToRealValue(str, &constant, &endptr) )
      str = endptr;
   else
      constant = 0.0;

   nvars = 0;
   varssize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars,    varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs,   varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &offsets, varssize) );

   /* read '+ (coef*(var+offset))^2' on lhs, as long as possible */
   while( *str != '\0' )
   {
      /* skip whitespace */
      while( isspace((int)*str) )
         ++str;

      /* stop if no more coefficients */
      if( strncmp(str, "+ (", 3) != 0 )
         break;

      str += 3;

      /* parse coef */
      if( !SCIPstrToRealValue(str, &coef, &endptr) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected coefficient at begin of '%s'\n", str);
         *success = FALSE;
         break;
      }
      str = endptr;

      if( strncmp(str, "*(", 2) != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected '*(' at begin of '%s'\n", str);
         *success = FALSE;
         break;
      }
      str += 2;

      /* parse variable name */
      SCIP_CALL( SCIPparseVarName(scip, str, &var, &endptr) );
      if( var == NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable name at '%s'\n", str);
         *success = FALSE;
         break;
      }
      str = endptr;

      /* parse offset */
      if( !SCIPstrToRealValue(str, &offset, &endptr) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected offset at begin of '%s'\n", str);
         *success = FALSE;
         break;
      }
      str = endptr;

      if( strncmp(str, "))^2", 4) != 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected '))^2' at begin of '%s'\n", str);
         *success = FALSE;
         break;
      }
      str += 4;

      if( varssize <= nvars )
      {
         varssize = SCIPcalcMemGrowSize(scip, varssize+1);
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars,    varssize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &coefs,   varssize) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &offsets, varssize) );
      }
      vars[nvars]    = var;
      coefs[nvars]   = coef;
      offsets[nvars] = offset;
      ++nvars;
   }

   if( strncmp(str, ") <=", 4) != 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected ') <=' at begin of '%s'\n", str);
      *success = FALSE;
   }
   str += 4;

   /* read rhs coef*(var+offset) or just a constant */

   /* parse coef */
   if( *success )
   {
      /* skip whitespace */
      while( isspace((int)*str) )
         ++str;

      if( !SCIPstrToRealValue(str, &rhscoef, &endptr) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected coefficient at begin of '%s'\n", str);
         *success = FALSE;
      }
      str = endptr;

      /* skip whitespace */
      while( isspace((int)*str) )
         ++str;
   }

   /* parse *(var+offset) */
   if( *str != '\0' )
   {
      if( *success )
      {
         if( strncmp(str, "*(", 2) != 0 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected '*(' at begin of '%s'\n", str);
            *success = FALSE;
         }
         else
         {
            str += 2;
         }
      }

      /* parse variable name */
      if( *success )
      {
         SCIP_CALL( SCIPparseVarName(scip, str, &rhsvar, &endptr) );
         if( rhsvar == NULL )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable name at '%s'\n", str);
            *success = FALSE;
         }
         else
         {
            str = endptr;
         }
      }

      /* parse offset */
      if( *success )
      {
         if( !SCIPstrToRealValue(str, &rhsoffset, &endptr) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected offset at begin of '%s'\n", str);
            *success = FALSE;
         }
         else
         {
            str = endptr;
         }
      }

      if( *success )
      {
         if( *str != ')' )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected ')' at begin of '%s'\n", str);
            *success = FALSE;
         }
      }
   }
   else if( *success )
   {
      /* only a constant at right hand side */
      rhsoffset = rhscoef;  /*lint !e644*/
      rhscoef = 1.0;
      rhsvar = NULL;
   }

   if( *success )
   {
      assert(!stickingatnode);
      SCIP_CALL( SCIPcreateConsSOC(scip, cons, name, nvars, vars, coefs, offsets, constant, rhsvar, rhscoef, rhsoffset,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable) );  /*lint !e644 */
   }

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &offsets);

   return SCIP_OKAY;
}
#else
#define consParseSOC NULL
#endif

/*
 * constraint specific interface methods
 */

/** creates the handler for second order cone constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSOC(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   conshdlrdata->subnlpheur = NULL;
   conshdlrdata->trysolheur = NULL;

   SCIP_CALL( SCIPincludeEventhdlr(scip, CONSHDLR_NAME"_boundchange",
         "signals a bound change to a second order cone constraint",
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, processVarEvent, NULL) );
   conshdlrdata->eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_boundchange");

   SCIP_CALL( SCIPincludeEventhdlr(scip, CONSHDLR_NAME"_newsolution",
         "handles the event that a new primal solution has been found",
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, processNewSolutionEvent, NULL) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopySOC, consFreeSOC, consInitSOC, consExitSOC,
         consInitpreSOC, consExitpreSOC, consInitsolSOC, consExitsolSOC,
         consDeleteSOC, consTransSOC, consInitlpSOC,
         consSepalpSOC, consSepasolSOC, consEnfolpSOC, consEnfopsSOC, consCheckSOC,
         consPropSOC, consPresolSOC, consRespropSOC, consLockSOC,
         consActiveSOC, consDeactiveSOC,
         consEnableSOC, consDisableSOC,
         consDelvarsSOC, consPrintSOC, consCopySOC, consParseSOC,
         conshdlrdata) );

   if( SCIPfindConshdlr(scip,"quadratic") != NULL )
   {
      /* notify function that upgrades quadratic constraint to SOC's */
      SCIP_CALL( SCIPincludeQuadconsUpgrade(scip, upgradeConsQuadratic, QUADCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );
   }

   /* add soc constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/scaling",
         "whether a constraint should be scaled w.r.t. the current gradient norm when checking for feasibility",
         &conshdlrdata->doscaling,        TRUE,  TRUE,          NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/projectpoint",
         "whether the reference point of a cut should be projected onto the feasible set of the SOC constraint",
         &conshdlrdata->projectpoint,     TRUE,  FALSE,         NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam (scip, "constraints/"CONSHDLR_NAME"/nauxvars",
         "number of auxiliary variables to use when creating a linear outer approx. of a SOC3 constraint; 0 to turn off",
         &conshdlrdata->nauxvars,         FALSE, 0, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/glineur",
         "whether the Glineur Outer Approximation should be used instead of Ben-Tal Nemirovski",
         &conshdlrdata->glineur,          FALSE, TRUE,          NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/minefficacy",
         "minimal efficacy of a cut to be added to LP in separation",
         &conshdlrdata->minefficacy,      FALSE, 0.0001, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/sparsify",
         "whether to sparsify cuts",
         &conshdlrdata->sparsify,         TRUE,  FALSE,         NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/sparsifymaxloss",
         "maximal loss in cut efficacy by sparsification",
         &conshdlrdata->sparsifymaxloss,  TRUE,  0.2, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/sparsifynzgrowth",
         "growth rate of maximal allowed nonzeros in cuts in sparsification",
         &conshdlrdata->sparsifynzgrowth, TRUE,  1.3, 1.000001, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/linfeasshift",
         "whether to try to make solutions feasible in check by shifting the variable on the right hand side",
         &conshdlrdata->linfeasshift,     FALSE, TRUE,          NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "constraints/"CONSHDLR_NAME"/nlpform",
         "which formulation to use when adding a SOC constraint to the NLP (a: automatic, q: nonconvex quadratic form, s: convex sqrt form, e: convex exponential-sqrt form, d: convex division form)",
         &conshdlrdata->nlpform,          FALSE, 'a', "aqsed", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/sepanlpmincont",
         "minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation",
         &conshdlrdata->sepanlpmincont, FALSE, 1.0, 0.0, 2.0, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a second order cone constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables on left hand side of constraint (n) */
   SCIP_VAR**            vars,               /**< array with variables on left hand side (x_i) */
   SCIP_Real*            coefs,              /**< array with coefficients of left hand side variables (alpha_i), or NULL if all 1.0 */
   SCIP_Real*            offsets,            /**< array with offsets of variables (beta_i), or NULL if all 0.0 */
   SCIP_Real             constant,           /**< constant on left hand side (gamma) */
   SCIP_VAR*             rhsvar,             /**< variable on right hand side of constraint (x_{n+1}) */
   SCIP_Real             rhscoeff,           /**< coefficient of variable on right hand side (alpha_{n+1}) */
   SCIP_Real             rhsoffset,          /**< offset of variable on right hand side (beta_{n+1}) */
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
   SCIP_Bool             removable           /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(modifiable == FALSE); /* we do not support column generation */

   /* find the soc constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("SOC constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   assert(vars     != NULL);
   assert(nvars    >= 2);
   assert(constant >= 0.0);
   assert(!SCIPisInfinity(scip, ABS(rhsoffset)));
   assert(!SCIPisInfinity(scip, constant));
   assert(rhsvar == NULL || rhscoeff <= 0.0 || SCIPisGE(scip, local ? SCIPcomputeVarLbLocal(scip, rhsvar) : SCIPcomputeVarLbGlobal(scip, rhsvar), -rhsoffset));
   assert(rhsvar == NULL || rhscoeff >= 0.0 || SCIPisLE(scip, local ? SCIPcomputeVarUbLocal(scip, rhsvar) : SCIPcomputeVarUbGlobal(scip, rhsvar), -rhsoffset));

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->nvars = nvars;
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->vars, vars, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPcaptureVar(scip, vars[i]) );
   }

   if( coefs != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->coefs, coefs, nvars) );
      for( i = 0; i < nvars; ++i )
         if( consdata->coefs[i] < 0.0 )
            consdata->coefs[i] = -consdata->coefs[i];
   }
   else
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->coefs, nvars) );
      for( i = 0; i < nvars; ++i )
         consdata->coefs[i] = 1.0;
   }

   if( offsets != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->offsets, offsets, nvars) );
   }
   else
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->offsets, nvars) );
      BMSclearMemoryArray(consdata->offsets, nvars);
   }

   consdata->constant  = constant;
   consdata->rhsvar    = rhsvar;
   consdata->rhscoeff  = rhscoeff;
   consdata->rhsoffset = rhsoffset;

   if( rhsvar != NULL )
   {
      SCIP_CALL( SCIPcaptureVar(scip, rhsvar) );
   }

   consdata->nlrow = NULL;

   consdata->lhsbndchgeventdatas = NULL;
   consdata->ispropagated        = FALSE;
   consdata->isapproxadded       = FALSE;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, FALSE) );

   if( SCIPisTransformed(scip) )
   {
      SCIP_CONSHDLRDATA* conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, *cons) );
   }

   return SCIP_OKAY;
}

/** Gets the SOC constraint as a nonlinear row representation.
 */
SCIP_RETCODE SCIPgetNlRowSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< a buffer where to store pointer to nonlinear row */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons  != NULL);
   assert(nlrow != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow == NULL )
   {
      SCIP_CALL( createNlRow(scip, SCIPconsGetHdlr(cons), cons) );
   }
   assert(consdata->nlrow != NULL);
   *nlrow = consdata->nlrow;

   return SCIP_OKAY;
}

/** Gets the number of variables on the left hand side of a SOC constraint.
 */
int SCIPgetNLhsVarsSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->nvars;
}

/** Gets the variables on the left hand side of a SOC constraint.
 */
SCIP_VAR** SCIPgetLhsVarsSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->vars;
}

/** Gets the coefficients of the variables on the left hand side of a SOC constraint, or NULL if all are equal to 1.0.
 */
SCIP_Real* SCIPgetLhsCoefsSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->coefs;
}

/** Gets the offsets of the variables on the left hand side of a SOC constraint, or NULL if all are equal to 0.0.
 */
SCIP_Real* SCIPgetLhsOffsetsSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->offsets;
}

/** Gets the constant on the left hand side of a SOC constraint.
 */
SCIP_Real SCIPgetLhsConstantSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->constant;
}

/** Gets the variable on the right hand side of a SOC constraint.
 */
SCIP_VAR* SCIPgetRhsVarSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->rhsvar;
}

/** Gets the coefficient of the variable on the right hand side of a SOC constraint.
 */
SCIP_Real SCIPgetRhsCoefSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->rhscoeff;
}

/** Gets the offset of the variables on the right hand side of a SOC constraint.
 */
SCIP_Real SCIPgetRhsOffsetSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->rhsoffset;
}

/** Adds the constraint to an NLPI problem.
 * Uses nonconvex formulation as quadratic function.
 */
SCIP_RETCODE SCIPaddToNlpiProblemSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SOC constraint */
   SCIP_NLPI*            nlpi,               /**< interface to NLP solver */
   SCIP_NLPIPROBLEM*     nlpiprob,           /**< NLPI problem where to add constraint */
   SCIP_HASHMAP*         scipvar2nlpivar,    /**< mapping from SCIP variables to variable indices in NLPI */
   SCIP_Bool             names               /**< whether to pass constraint names to NLPI */
   )
{
   SCIP_CONSDATA* consdata;
   int            nlininds;
   int*           lininds;
   SCIP_Real*     linvals;
   int            nquadelems;
   SCIP_QUADELEM* quadelems;
   int            j;
   int            lincnt;
   SCIP_Real      lhs;
   SCIP_Real      rhs;
   const char*    name;

   assert(scip     != NULL);
   assert(cons     != NULL);
   assert(nlpi     != NULL);
   assert(nlpiprob != NULL);
   assert(scipvar2nlpivar != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   lhs = -SCIPinfinity(scip);
   rhs = -consdata->constant;

   /* count how length is the linear part, i.e., how many offsets we have */
   nlininds = consdata->rhsoffset != 0.0 ? 1 : 0;
   for( j = 0; j < consdata->nvars; ++j )
      if( consdata->offsets[j] != 0.0 )
         ++nlininds;

   lininds = NULL;
   linvals = NULL;
   if( nlininds )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nlininds) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nlininds) );
   }
   lincnt = 0;

   nquadelems = consdata->nvars + 1;
   SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nquadelems) );

   for( j = 0; j < consdata->nvars; ++j )
   {
      quadelems[j].idx1 = (int) (size_t) SCIPhashmapGetImage(scipvar2nlpivar, consdata->vars[j]);
      quadelems[j].idx2 = quadelems[j].idx1;
      quadelems[j].coef = consdata->coefs[j] * consdata->coefs[j];

      if( consdata->offsets[j] != 0.0 )
      {
         assert(lininds != NULL);
         assert(linvals != NULL);
         lininds[lincnt] = quadelems[j].idx1;
         linvals[lincnt] = 2 * quadelems[j].coef * consdata->offsets[j];
         ++lincnt;

         rhs -= quadelems[j].coef * consdata->offsets[j] * consdata->offsets[j];
      }
   }
   quadelems[consdata->nvars].idx1 = (int) (size_t) SCIPhashmapGetImage(scipvar2nlpivar, consdata->rhsvar);
   quadelems[consdata->nvars].idx2 = quadelems[consdata->nvars].idx1;
   quadelems[consdata->nvars].coef = - consdata->rhscoeff * consdata->rhscoeff;

   if( consdata->rhsoffset != 0.0 )
   {
      assert(lininds != NULL);
      assert(linvals != NULL);
      lininds[lincnt] = quadelems[consdata->nvars].idx1;
      linvals[lincnt] = -2.0 * consdata->rhscoeff * consdata->rhscoeff * consdata->rhsoffset;
      ++lincnt;

      rhs += consdata->rhscoeff * consdata->rhscoeff * consdata->rhsoffset * consdata->rhsoffset;
   }

   assert(lincnt == nlininds);

   name = names ? SCIPconsGetName(cons) : NULL;

   SCIP_CALL( SCIPnlpiAddConstraints(nlpi, nlpiprob, 1,
         &lhs, &rhs,
         &nlininds, &lininds, &linvals,
         &nquadelems, &quadelems,
         NULL, NULL, &name) );

   SCIPfreeBufferArrayNull(scip, &lininds);
   SCIPfreeBufferArrayNull(scip, &linvals);
   SCIPfreeBufferArray(scip, &quadelems);

   return SCIP_OKAY;
}
