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

/**@file   sepa_oddcycle.c
 * @brief  oddcycle separator
 * @author Robert Waniek
 *
 * We separate odd cycle inequalities in the implication graph. Implemented are the classic method
 * by Groetschel, Lovasz, and Schrijver (GLS) and the levelgraph method by Hoffman and Padberg (HP)
 *
 * Odd cycle inequalities are lifted by a heuristic method based on an idea from Alvarez-Valdes,
 * Parreno, and Tamarit.
 *
 * Some part of this code is based on the odd cycle separator of the program colorbitopt by Marc
 * Pfetsch.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_oddcycle.h"
#include "scip/pub_misc.h"
#include "dijkstra/dijkstra_bh.h"


#define SEPA_NAME              "oddcycle"
#define SEPA_DESC              "odd cycle separator"
#define SEPA_PRIORITY            -15000
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE      /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE      /**< should separation method be delayed, if other separators found cuts? */


/* default values for separator settings */
#define DEFAULT_SCALEFACTOR       1000      /**< factor for scaling of the arc-weights in the Dijkstra algorithm */
#define DEFAULT_USECLASSICAL      TRUE      /**< use GLS method, otherwise HP method */
#define DEFAULT_LIFT              FALSE     /**< lift odd cycle cuts */
#define DEFAULT_REPAIRCYCLES      TRUE      /**< try to repair violated cycles in which a variable and its negated appear */
#define DEFAULT_ADDSELFARCS       TRUE      /**< add links between a variable and its negated */
#define DEFAULT_INCLUDETRIANGLES  TRUE      /**< separate triangles (3-cliques) found as 3-cycles or repaired larger cycles */
#define DEFAULT_MULTIPLECUTS      FALSE     /**< even if a variable is already covered by a cut, still try it as start node */
#define DEFAULT_ALLOWMULTIPLECUTS TRUE      /**< even if a variable is already covered by a cut, still allow another cut to cover it too */
#define DEFAULT_LPLIFTCOEF        FALSE     /**< TRUE: choose lifting candidate with highest value of coefficient*lpvalue
                                             *   FALSE: choose lifting candidate with highest coefficient */
#define DEFAULT_RECALCLIFTCOEF    TRUE      /**< whether lifting coefficients should be recomputed */
#define DEFAULT_MAXSEPACUTS       5000      /**< maximal number of oddcycle cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT   5000      /**< maximal number of oddcycle cuts separated per separation round in root node */
#define DEFAULT_PERCENTTESTVARS   0         /**< percent of variables to try the chosen method on [0-100] */
#define DEFAULT_OFFSETTESTVARS    100       /**< offset of variables to try the chosen method on */
#define DEFAULT_MAXCUTSROOT       1         /**< maximal number of oddcycle cuts generated per root of the levelgraph */
#define DEFAULT_SORTSWITCH        3         /**< unsorted (0), maxlp (1), minlp (2), maxfrac (3), minfrac (4) */
#define DEFAULT_MAXREFERENCE      0         /**< minimal weight on an edge (in level graph or Dijkstra graph) */
#define DEFAULT_MAXROUNDS         10        /**< maximal number of rounds pre node */
#define DEFAULT_MAXROUNDSROOT     10        /**< maximal number of rounds in the root node */
#define DEFAULT_MAXNLEVELS        20        /**< maximal number of levels in level graph */
#define DEFAULT_MAXPERNODESLEVEL  100       /**< maximal percentage of nodes allowed in one level of the levelgraph [0-100] */
#define DEFAULT_OFFSETNODESLEVEL  10        /**< additional offset of nodes allowed in one level of the levelgraph */
#define DEFAULT_SORTROOTNEIGHBORS TRUE      /**< sort neighbors of the root in the level graph */
#define DEFAULT_MAXCUTSLEVEL      50        /**< maximal number of cuts produced per level */


/*
 * Data structures
 */

/** Graph structure for level graph
 *
 *  This graph is tailored to the heuristic search for odd holes, @see separateHeur().
 *
 *  This undirected graph is represented by a directed graph with forward and backward arcs. Arcs are
 *  forward if they lead from a level l to level l+1, i.e., away from the root; backward arcs
 *  lead from a level l+1 to level l. This distinction enables a fast construction and search
 *  process. In the latter only forward or backward arcs have to be searched.
 *
 *  Target nodes and weights of the arcs incident to each node (adjacency lists) are stored
 *  consecutively in the arrays targetForward, targetBackward, weightForward, and weightBackward.
 *  The end of each list is marked by a -1 in targetForward and targetBackward.
 */
struct levelGraph
{
   unsigned int          nnodes;             /**< number of nodes */
   unsigned int          nedges;             /**< number of arcs */
   unsigned int          n;                  /**< maximal number of nodes of the level graph */
   unsigned int          m;                  /**< maximal number of arcs of the level graph */
   unsigned int          nlevels;            /**< number of levels completely inserted so far */
   unsigned int*         level;              /**< level number for each node */
   unsigned int          lastF;              /**< index of last storage element (in targetForward, weightForward) for forward direction */
   unsigned int          lastB;              /**< index of last storage element (in targetBackward, weightBackward) for backward direction */
   int*                  beginForward;       /**< index of forward adjacency list (in targetForward, weightForward) for each node */
   int*                  beginBackward;      /**< index of backward adjacency list (in targetBackward, weightBackward) for each node */
   int*                  targetForward;      /**< target nodes of forward arcs */
   int*                  targetBackward;     /**< target nodes of backward arcs */
   unsigned int*         weightForward;      /**< weights of forward arcs */
   unsigned int*         weightBackward;     /**< weights of backwards arcs */
   unsigned int          sizeForward;        /**< size of targetForward and weightForward */
   unsigned int          sizeBackward;       /**< size of targetBackward and weightBackward */
   int*                  beginAdj;           /**< index of list of arcs inside a level (in sourceAdj) for each node
                                              *   (the index points at the first arc starting from this node) */
   unsigned int*         sourceAdj;          /**< source nodes of arcs inside a level */
   unsigned int*         targetAdj;          /**< target nodes of arcs inside a level */
   unsigned int*         weightAdj;          /**< weights of arcs inside a level */
   unsigned int*         levelAdj;           /**< index of the first arc inside a given level */
   unsigned int          sizeAdj;            /**< size of sourceAdj, targetAdj and weightAdj */
};

typedef struct levelGraph LEVELGRAPH;


/** sorting type for starting node or root node iteration order
 *
 *  If the array should be sorted (1-4), the variable array is sorted every round by the chosen
 *  sorttype and the search method tries the nodes in order of the array.  If the array is used
 *  unsorted (0), the search methods tries the nodes in order of the array and stores the last
 *  processed start node or root node and continues from this position in the next separation round.
 */
enum sorttype
{
   UNSORTED              = 0,                /**< variable array is unsorted */
   MAXIMAL_LPVALUE       = 1,                /**< variable array is sorted by maximal lp-value */
   MINIMAL_LPVALUE       = 2,                /**< variable array is sorted by minimal fractionality */
   MAXIMAL_FRACTIONALITY = 3,                /**< variable array is sorted by maximal lp-value */
   MINIMAL_FRACTIONALITY = 4                 /**< variable array is sorted by minimal fractionality */
};
typedef enum sorttype SORTTYPE;


/** separator data */
struct SCIP_SepaData
{
   int                   scale;              /**< factor for scaling of the arc-weights */
   unsigned int          ncuts;              /**< number of cuts, added by the separator so far (in current and past calls) */
   unsigned int          oldncuts;           /**< number of cuts at the start the current separation round */
   int                   nliftedcuts;        /**< number of lifted cuts, added by the separator so far (in current and past calls) */
   SCIP_Bool             useclassical;       /**< use GLS method, otherwise HP method */
   SCIP_Bool             multiplecuts;       /**< an odd cycle cut of length L can be generated L times; forbidding multiple cuts
                                              *   per node might be faster but might miss some cuts in the current round */
   SCIP_Bool             allowmultiplecuts;  /**< allow multiple cuts covering one node */
   SCIP_Bool             liftoddcycles;      /**< TRUE, iff we try to lift odd cycle inequalities */
   SCIP_Bool             addselfarcs;        /**< add arcs between the nodes of a variable and its negated; since not all implications
                                              *   are in the graph, this often finds more cycles */
   SCIP_Bool             repaircycles;       /**< if a variable and its negated appear in a cycle, we can repair the cycle
                                              *   by removing both and reconnecting the remaining nodes of the cycle */
   SCIP_Bool             includetriangles;   /**< handle triangles found as 3-cycles or repaired larger cycles */
   LEVELGRAPH*           levelgraph;         /**< level graph when using HP method, NULL otherwise */
   Dijkstra_Graph*       dijkstragraph;      /**< Dijkstra graph if using method by GLS, NULL otherwise */
   unsigned int*         mapping;            /**< mapping for getting the index of a variable in the sorted variable array */
   SCIP_Bool             lpliftcoef;         /**< TRUE: choose lifting candidate with highest value of coefficient*lpvalue
                                              *   FALSE: choose lifting candidate with highest coefficient */
   SCIP_Bool             recalcliftcoef;     /**< whether lifting coefficients should be recomputed */
   int                   maxsepacuts;        /**< maximal number of oddcycle cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of oddcycle cuts separated per separation round in the root node */
   int                   maxsepacutsround;   /**< maximal number of oddcycle cuts separated per separation round in the current node */
   SORTTYPE              sortswitch;         /*   sorted variable array: unsorted (0), maxlp (1), minlp (2), maxfrac (3), minfrac (4) */
   int                   lastroot;           /**< save root of last GLS-method run */
   SCIP_Bool             sortrootneighbors;  /**< sort neighbors of the root in the level graph */
   int                   percenttestvars;    /**< percentage of variables to try the chosen method on [0-100] */
   int                   offsettestvars;     /**< offset of variables to try the chosen method on */
   int                   maxpernodeslevel;   /**< percentage of nodes allowed in the same level of the level graph [0-100] */
   int                   offsetnodeslevel;   /**< additional offset of nodes allowed in one level of the levelgraph */
   unsigned int          maxlevelsize;       /**< maximal number of nodes allowed in the same level of the level graph */
   int                   maxcutsroot;        /**< maximal number of oddcycle cuts generated per root of the levelgraph */
   int                   maxcutslevel;       /**< maximal number of oddcycle cuts generated per level of the level graph */
   int                   maxrounds;          /**< maximal number of oddcycle separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of oddcycle separation rounds in the root node (-1: unlimited) */
   int                   maxreference;       /**< minimal weight on an edge (in level graph or Dijkstra graph) */
   int                   maxnlevels;         /**< maximal number of levels in level graph */
};



/*
 * Local methods
 */



/*
 * debugging methods
 */

#ifdef SCIP_DEBUG

/** displays cycle of pred data structure w.r.t. variable names of the original problem (including
 *  status: original or negated node in graph)
 */
static
void printCycle(
   SCIP_VAR**            vars,               /**< problem variables */
   unsigned int*         pred,               /**< cycle stored as predecessor list */
   unsigned int          nbinvars,           /**< number of binary problem variables */
   unsigned int          startnode           /**< a node of the cycle */
   )
{
   unsigned int varsindex;
   unsigned int counter;

   assert(vars != NULL);
   assert(pred != NULL);
   assert(nbinvars > 0);
   assert(startnode < 4*nbinvars);

   counter = 0;
   varsindex = startnode;

   /* print start/end node */
   if( varsindex < nbinvars || ( varsindex >= 2*nbinvars && varsindex < 3*nbinvars ) )
   {
      SCIPdebugMessage("+ %s\n",SCIPvarGetName(vars[varsindex%nbinvars]));
   }
   else
   {
      SCIPdebugMessage("- %s\n",SCIPvarGetName(vars[varsindex%nbinvars]));
   }

   /* print inner nodes */
   for( varsindex = pred[startnode]; varsindex != startnode; varsindex = pred[varsindex] )
   {
      if( varsindex < nbinvars || ( varsindex >= 2*nbinvars && varsindex < 3*nbinvars ) )
      {
         SCIPdebugMessage("+ %s\n",SCIPvarGetName(vars[varsindex%nbinvars]));
      }
      else
      {
         SCIPdebugMessage("- %s\n",SCIPvarGetName(vars[varsindex%nbinvars]));
      }
      ++counter;
   }

   /* print start/end node */
   if( varsindex < nbinvars || ( varsindex >= 2*nbinvars && varsindex < 3*nbinvars ) )
   {
      SCIPdebugMessage("+ %s\n",SCIPvarGetName(vars[varsindex%nbinvars]));
   }
   else
   {
      SCIPdebugMessage("- %s\n",SCIPvarGetName(vars[varsindex%nbinvars]));
   }

   ++counter;
   SCIPdebugMessage("original cycle has %d variables.\n",counter);
}
#endif


/*
 * lifting methods
 */

/** using the level graph (if possible) or Dijkstra graph data structure (depending on the used
 *  method) we determine whether two nodes are adjacent
 */
static
SCIP_Bool isNeighbor(
   SCIP_VAR**            vars,               /**< problem variables */
   unsigned int          nbinvars,           /**< number of binary problem variables */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   unsigned int          a,                  /**< node index of first variable */
   unsigned int          b                   /**< node index of second variable */
   )
{
   unsigned int i;

   assert(vars != NULL);
   assert(nbinvars > 2);
   assert(sepadata != NULL);
   assert(sepadata->levelgraph != NULL || sepadata->useclassical);
   assert(sepadata->dijkstragraph != NULL || !(sepadata->useclassical));
   assert(a < 2*nbinvars);
   assert(b < 2*nbinvars);
   assert(a != b);

   /* determine adjacency using the Dijkstra graph */
   if( sepadata->useclassical )
   {
      if( sepadata->dijkstragraph->outcnt[a] == 0 || sepadata->dijkstragraph->outcnt[b] == 0 )
         return FALSE;

      /* @todo later: if helpful: sort head and weight list once */
      for( i = sepadata->dijkstragraph->outbeg[a]; i < sepadata->dijkstragraph->outbeg[a] + sepadata->dijkstragraph->outcnt[a]; ++i )
      {
         if( sepadata->dijkstragraph->head[i] == b + 2*nbinvars )
            return TRUE;
      }
   }
   else    /* determine adjacency using the level graph */
   {
      /* if a and b are contained in the level graph (with their arcs), we can check inside the level graph structure */
      if( (sepadata->levelgraph->beginForward[a] != -1 || sepadata->levelgraph->beginBackward[a] != -1)
         && (sepadata->levelgraph->beginForward[b] != -1 || sepadata->levelgraph->beginBackward[b] != -1) )
      {
         assert(sepadata->levelgraph->level[a] <= sepadata->levelgraph->nlevels);
         assert(sepadata->levelgraph->level[b] <= sepadata->levelgraph->nlevels);

         /* if a and b are not in neighbored levels or the same level, they cannot be adjacent */
         if( sepadata->levelgraph->level[a] > sepadata->levelgraph->level[b] + 1
            || sepadata->levelgraph->level[b] > sepadata->levelgraph->level[a] + 1 )
            return FALSE;

         assert(sepadata->levelgraph->level[a] == sepadata->levelgraph->level[b]
            || sepadata->levelgraph->level[a]+1 == sepadata->levelgraph->level[b]
            || sepadata->levelgraph->level[a] == sepadata->levelgraph->level[b]+1);

         /* first case of adjacent level */
         if( sepadata->levelgraph->level[a] == sepadata->levelgraph->level[b]+1 )
         {
            if( sepadata->levelgraph->beginBackward[a] >= 0 )
            {
               i = (unsigned int) sepadata->levelgraph->beginBackward[a];
               while( sepadata->levelgraph->targetBackward[i] != -1 )
               {
                  if( sepadata->levelgraph->targetBackward[i] == (int)b )
                     return TRUE;
                  ++i;
               }
            }
         }
         else if( sepadata->levelgraph->level[a] == sepadata->levelgraph->level[b]-1 )    /* second case of adjacent level */
         {
            if( sepadata->levelgraph->beginForward[a] >= 0 )
            {
               i = (unsigned int) sepadata->levelgraph->beginForward[a];
               while( sepadata->levelgraph->targetForward[i] != -1 )
               {
                  if( sepadata->levelgraph->targetForward[i] == (int)b )
                     return TRUE;
                  ++i;
               }
            }
         }
         else          /* same level (note that an edge between a and b is stored for a if a < b, otherwise it is stored for b) */
         {
            assert(sepadata->levelgraph->level[a] == sepadata->levelgraph->level[b]);
            assert(sepadata->levelgraph->level[a] > 0); /* root has no neighbor in the same level */

            if( a < b && sepadata->levelgraph->beginAdj[a] >= 0 )
            {
               i = (unsigned int) sepadata->levelgraph->beginAdj[a];
               assert(i >= sepadata->levelgraph->levelAdj[sepadata->levelgraph->level[a]]);

               while( sepadata->levelgraph->sourceAdj[i] == a && i < sepadata->levelgraph->levelAdj[sepadata->levelgraph->level[a]+1] )
               {
                  if( sepadata->levelgraph->targetAdj[i] == b )
                     return TRUE;

                  /* if adjacency list ends we are done and a and b are not adjacent */
                  if( sepadata->levelgraph->sourceAdj[i] == 0 && sepadata->levelgraph->targetAdj[i] == 0 )
                     return FALSE;

                  assert(sepadata->levelgraph->sourceAdj[i] < sepadata->levelgraph->targetAdj[i]);
                  ++i;
               }
            }
            if( b < a && sepadata->levelgraph->beginAdj[b] >= 0 )
            {
               i = (unsigned int) sepadata->levelgraph->beginAdj[b];
               assert(i >= sepadata->levelgraph->levelAdj[sepadata->levelgraph->level[b]]);

               while( sepadata->levelgraph->sourceAdj[i] == b && i < sepadata->levelgraph->levelAdj[sepadata->levelgraph->level[b]+1])
               {
                  if( sepadata->levelgraph->targetAdj[i] == a )
                     return TRUE;

                  /* if adjacency list ends we are done and a and b are not adjacent */
                  if( sepadata->levelgraph->sourceAdj[i] == 0 && sepadata->levelgraph->targetAdj[i] == 0 )
                     return FALSE;

                  assert(sepadata->levelgraph->sourceAdj[i] < sepadata->levelgraph->targetAdj[i]);
                  ++i;
               }
            }
         }
      }
      /* if a or b is not in the levels already completely inserted in the levelgraph, we check
       * their adjacency by the SCIP data structures
       */
      else
      {
         SCIP_Bool originala;
         SCIP_Bool originalb;

         unsigned int nbinimpls;
         SCIP_VAR** implvars;
         SCIP_BOUNDTYPE* impltypes;
#ifndef NDEBUG
         SCIP_Real* implbounds;
#endif
         unsigned int ncliques;
         SCIP_CLIQUE** cliques;
         unsigned int ncliquevars;
         SCIP_VAR** cliquevars;
         SCIP_Bool* cliquevals;

         unsigned int j;

         /* get original variables */
         originala = TRUE;
         if( a >= nbinvars )
         {
            a = a - nbinvars;
            originala = FALSE;
         }
         assert(a < nbinvars);

         originalb = TRUE;
         if( b >= nbinvars )
         {
            b = b - nbinvars;
            originalb = FALSE;
         }
         assert(b < nbinvars);

         /* nodes cannot be connected by trivial observations */
         if( ( SCIPvarGetNBinImpls(vars[a], originala) + SCIPvarGetNCliques(vars[a], originala) == 0 )
            || ( SCIPvarGetNBinImpls(vars[b], originalb) + SCIPvarGetNCliques(vars[b], originalb) == 0 ) )
            return FALSE;
         if( ( SCIPvarGetNBinImpls(vars[b], originalb) == 0 && SCIPvarGetNCliques(vars[a], originala) == 0 )
            || ( SCIPvarGetNBinImpls(vars[a], originala) == 0 && SCIPvarGetNCliques(vars[b], originalb) == 0 ) )
            return FALSE;

         /* @todo later: possible improvement: do this test for implications and cliques separately if this here is time consuming */
         /* one of the nodes seems to have more arcs than the other, we swap them (since adjacency is symmetric) */
         if( SCIPvarGetNBinImpls(vars[a], originala) + 2 * SCIPvarGetNCliques(vars[a], originala) >
            SCIPvarGetNBinImpls(vars[b], originalb) + 2 * SCIPvarGetNCliques(vars[b], originalb) )
         {
            unsigned int temp;
            SCIP_Bool varfixingtemp;

            temp = b;
            varfixingtemp = originalb;
            b = a;
            originalb = originala;
            a = temp;
            originala = varfixingtemp;
         }

         /* check whether there is an implication a = 1 -> b = 0 */
         nbinimpls = (unsigned int) SCIPvarGetNBinImpls(vars[a], originala);
         implvars = SCIPvarGetImplVars(vars[a], originala);
         impltypes = SCIPvarGetImplTypes(vars[a], originala);
#ifndef NDEBUG
         implbounds = SCIPvarGetImplBounds(vars[a], originala);
#endif

         assert(implvars != NULL || nbinimpls == 0);
         assert(impltypes != NULL || nbinimpls == 0);
         assert(implbounds != NULL || nbinimpls == 0);

         for( i = 0; i < nbinimpls; ++i )
         {
            assert( implvars != NULL && impltypes != NULL && implbounds != NULL ); /* for lint */
            if( SCIPvarGetProbindex(vars[b]) == SCIPvarGetProbindex(implvars[i]) )
            {
               if( impltypes[i] == SCIP_BOUNDTYPE_UPPER && originalb == TRUE )
               {
                  assert(implbounds[i] == 0.0);
                  return TRUE;
               }
               if( impltypes[i] == SCIP_BOUNDTYPE_LOWER && originalb == FALSE )
               {
                  assert(implbounds[i] == 1.0);
                  return TRUE;
               }
            }
         }

         /* check whether a and b are contained in a clique */
         ncliques =  (unsigned int) SCIPvarGetNCliques(vars[a], originala);
         cliques = SCIPvarGetCliques(vars[a], originala);

         assert(cliques != NULL || ncliques == 0);

         for( i = 0; i < ncliques; ++i )
         {
            assert( cliques != NULL );  /* for lint */
            ncliquevars = (unsigned int) SCIPcliqueGetNVars(cliques[i]);
            cliquevars = SCIPcliqueGetVars(cliques[i]);
            cliquevals = SCIPcliqueGetValues(cliques[i]);

            assert(cliquevars != NULL || ncliquevars == 0);
            assert(cliquevals != NULL || ncliquevars == 0);

            for( j = 0; j < ncliquevars; ++j )
            {
               assert( cliquevals != NULL &&  cliquevars != NULL ); /* for lint */
               if( SCIPvarGetProbindex(vars[b]) == SCIPvarGetProbindex(cliquevars[j]) )
               {
                  if( (cliquevals[j] == FALSE && originalb == TRUE) || ( cliquevals[j] == TRUE && originalb == FALSE ) )
                     return TRUE;
               }
            }
         }
      }
   }

   return FALSE;
}

/** inside the lifting heuristic we determine the lifting coefficient by counting the length of
 *  chains adjacent to the lifting candidate.
 *
 *  since we have to exclude all chains adjacent to an already lifted node which is not adjacent to
 *  the current lifting candidate we check all chains of the cycle of length three and block them if
 *  they are adjacent.
 */
static
void checkBlocking(
   unsigned int          a,                  /**< first node of the checked cycle chain of length 3 */
   unsigned int          b,                  /**< second node of the checked cycle chain of length 3 */
   unsigned int          c,                  /**< third node of the checked cycle chain of length 3 */
   unsigned int          i,                  /**< current lifting candidate */
   unsigned int*         cycle,              /**< list of cycle nodes in order of the cycle */
   unsigned int          ncyclevars,         /**< number of variables in the odd cycle */
   SCIP_VAR**            vars,               /**< problem variables */
   unsigned int          nbinvars,           /**< number of binary problem variables */
   unsigned int*         lifted,             /**< list of lifted nodes */
   unsigned int*         nlifted,            /**< number of lifted nodes */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_Bool*            myi                 /**< flag array, if cycle node is inner point of a counted chain */
   )
{
   unsigned int k;

   assert(a < ncyclevars);
   assert(b < ncyclevars);
   assert(c < ncyclevars);
   assert(cycle != NULL);
   assert(ncyclevars % 2 == 1);
   assert(ncyclevars > 2);
   assert(ncyclevars <= nbinvars);
   assert(vars != NULL);
   assert(nbinvars > 2);
   assert(lifted != NULL);
   assert(nlifted != NULL);
   assert(myi != NULL);

   k = 0;
   while( (myi[a] || myi[b] || myi[c]) && k < *nlifted )
   {
      /* if all three nodes are adjacent to a node which is already lifted and not adjacent with the
       * current lifting candidate, they cannot be regarded */
      if( !isNeighbor(vars, nbinvars, sepadata, i, lifted[k])
         && isNeighbor(vars, nbinvars, sepadata, cycle[a], lifted[k])
         && isNeighbor(vars, nbinvars, sepadata, cycle[b], lifted[k])
         && isNeighbor(vars, nbinvars, sepadata, cycle[c], lifted[k]) )
      {
         myi[a] = FALSE;
         myi[b] = FALSE;
         myi[c] = FALSE;
      }
      ++k;
   }
}

/** determine the heuristic lifting coefficient by counting the length of the adjacent chains of the
 *  candidate (we have to exclude all chains that are adjacent to an already lifted node which is
 *  not adjacent to the current candidate)
 */
static
unsigned int getCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int          i,                  /**< current lifting candidate */
   unsigned int*         cycle,              /**< list of cycle nodes in order of the cycle */
   unsigned int          ncyclevars,         /**< number of variables in the odd cycle */
   SCIP_VAR**            vars,               /**< problem variables */
   unsigned int          nbinvars,           /**< number of binary problem variables */
   unsigned int*         lifted,             /**< list of lifted nodes */
   unsigned int*         nlifted,            /**< number of lifted nodes */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_Bool*            myi                 /**< flag array, if cycle node is inner point of a counted chain */
   )
{
   int j;
   unsigned int k;
   unsigned int coef;                        /* coefficient of lifting candidate of the current step */
   unsigned int end;

   assert(scip != NULL);
   assert(i < 2*nbinvars);
   assert(cycle != NULL);
   assert(ncyclevars % 2 == 1);
   assert(ncyclevars > 2);
   assert(ncyclevars <= 2*nbinvars);
   assert(vars != NULL);
   assert(nbinvars > 2);
   assert(nlifted != NULL);
   assert(lifted != NULL);

   coef = 0;

   /* get inner nodes of adjacent chains in cycle */
   for( j = 1; j < (int)ncyclevars-1; ++j )
   {
      myi[j] = isNeighbor(vars, nbinvars, sepadata, i, cycle[j-1]) && isNeighbor(vars, nbinvars, sepadata, i, cycle[j])
         && isNeighbor(vars, nbinvars, sepadata, i, cycle[j+1]);
   }

   /* the first and last node of the cycle are treated separately */
   myi[0] = isNeighbor(vars, nbinvars, sepadata, i, cycle[ncyclevars-1]) && isNeighbor(vars, nbinvars, sepadata, i, cycle[0])
      && isNeighbor(vars, nbinvars, sepadata, i, cycle[1]);
   myi[ncyclevars-1] = isNeighbor(vars, nbinvars, sepadata, i, cycle[ncyclevars-2])
      && isNeighbor(vars, nbinvars, sepadata, i, cycle[ncyclevars-1])
      && isNeighbor(vars, nbinvars, sepadata, i, cycle[0]);

   /* consider already lifted nodes that are not adjacent to current lifting candidate and
    * remove all inner cycle nodes that are adjacent to them
    */
   for( j = 1; j < (int)ncyclevars-1; ++j )
   {
      checkBlocking((unsigned int) (j-1), (unsigned int) j, (unsigned int) (j+1), i, cycle, ncyclevars, vars, nbinvars, lifted, nlifted, sepadata, myi);
   }
   checkBlocking(ncyclevars-2, ncyclevars-1, 0, i, cycle, ncyclevars, vars, nbinvars, lifted, nlifted, sepadata, myi);
   checkBlocking(ncyclevars-1, 0, 1, i, cycle, ncyclevars, vars, nbinvars, lifted, nlifted, sepadata, myi);

   /* calculate lifting coefficient */
   k = 0;

   /* first, handle the special case, that the first node of the cycle list is part of a chain */
   if( myi[0] )
   {
      ++k;
      end = ncyclevars-1;
      while( myi[end] && end > 0 )
      {
         ++k;
         --end;
      }
      assert(k == ncyclevars || end > 0);

      /* all cycle nodes build a relevant chain (maximal chain s.t. all inner nodes are in myi) */
      if( end == 0 )
      {
         assert(ncyclevars % 2 == 1);
         coef = (ncyclevars-1)/2;
         return coef;
      }
      assert(!myi[end]);

      /* current nonempty relevant chain cannot be extended */
      if( !myi[1] )
      {
         coef = (unsigned int) SCIPfloor(scip,(k+1.0)/2.0);
         assert(coef <= (ncyclevars-1)/2);
         k = 0;
      }
   }
   else
      end = ncyclevars;

   /* find remaining relevant chains */
   j = 1;
   while( j < (int)end )
   {
      /* skip all nodes that are not inner node */
      while( !myi[j] && j<(int)end )
         ++j;

      /* collect all inner nodes (chain is extended) */
      while( myi[j] && j<(int)end )
      {
         ++k;
         ++j;
      }

      if( k > 0 )
      {
         assert(myi[j-1]);
         coef += (unsigned int) SCIPfloor(scip,(k+1.0)/2.0);
         assert(coef <= (ncyclevars-1)/2);
         k = 0;
      }
   }

   return coef;
}

/** Lifting Heuristic based on an idea by Alvarez-Valdes, Parreno, Tamarit
 *
 *  This method is based on the observation, that a non-cycle node can be lifted into the inequality
 *  with coefficient \f$1\f$ if the node is adjacent to the nodes of a 3-chain on the cycle.
 *
 *  The coefficient can be calculated as
 *  \f$\left\lfloor{\frac{|C|-1}{2}}\right\rfloor\f$
 *  where \f$C\f$ is the chain on the cycle.
 *
 *  If the node is connected to several chains, the coefficients of the chains can be summed up, resulting
 *  in a feasible lifting coefficient.
 *
 *  Additionally further variables can be lifted by considering chains connected to the additional lifting node
 *  which are not connected to already lifted nodes.
 *
 *  This method is a feasible heuristic which gives a valid lifted inequality.
 *  (Furthermore the first lifting coefficient is always smaller or equal to the largest possible lifting coefficient.)
 */
static
SCIP_RETCODE liftOddCycleCut(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int*         nlifted,            /**< number of lifted variables */
   unsigned int*         lifted,             /**< indices of the lifted variables */
   unsigned int*         liftcoef,           /**< lifting coefficients */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   unsigned int          nbinvars,           /**< number of binary problem variables */
   unsigned int          startnode,          /**< a node of the cycle */
   unsigned int*         pred,               /**< predecessor of each node (original and negated) in odd cycle */
   unsigned int          ncyclevars,         /**< number of variables in the odd cycle */
   SCIP_Real*            vals,               /**< values of the variables in the given solution */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   unsigned int* cycle;                      /* storage for cycle and lifted nodes (and their coefficients) */
   unsigned int* coef;
   SCIP_Bool* candList;                      /* lifting candidate list */
   unsigned int i;
   unsigned int j;
   unsigned int negated;
   int bestcand;
   unsigned int liftround;
   SCIP_Bool* myi;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(sepadata->levelgraph != NULL || sepadata->useclassical);
   assert(sepadata->dijkstragraph != NULL || !(sepadata->useclassical));
   assert(vars != NULL);
   assert(nbinvars > 2);
   assert(startnode < 2*nbinvars);
   assert(pred != NULL);
   assert(ncyclevars % 2 == 1);
   assert(ncyclevars > 2);
   assert(ncyclevars <= nbinvars);
   assert(result != NULL);
   assert(nlifted != NULL);
   assert(lifted != NULL);
   assert(liftcoef != NULL);

   /* allocate memory for cycle list */
   SCIP_CALL( SCIPallocBufferArray(scip, &cycle, (int) ncyclevars) );

   /* transform cycle from predecessor list to array in order of appearance in cycle */
   cycle[0] = startnode;
   j = 1;
   i = pred[startnode];
   while( i != startnode )
   {
      cycle[j] = i;
      i = pred[i];
      ++j;
   }
   assert(j == ncyclevars);

   /* allocate memory for coefficients of the lifting candidates (used in every step) */
   SCIP_CALL( SCIPallocBufferArray(scip, &coef, (int) (2*nbinvars)) );

   /* allocate memory candidate list and list of lifted nodes */
   SCIP_CALL( SCIPallocBufferArray(scip, &candList, (int) (2*nbinvars)) );

   /* allocate memory for counting of chains in getCoef() */
   SCIP_CALL( SCIPallocBufferArray(scip, &myi, (int) ncyclevars) );

   if( SCIPisStopped(scip) )
      goto TERMINATE;

   /* initialize candidate list */
   for( i = 0; i < 2*nbinvars; ++i )
      candList[i] = TRUE;

   /* remove cycle variables and their negated from candidate list */
   for( i = 0; i < ncyclevars; ++i )
   {
      candList[cycle[i]] = FALSE;
      if( cycle[i] >= nbinvars )
         negated = cycle[i] - nbinvars;
      else
         negated = cycle[i] + nbinvars;
      assert(negated < 2*nbinvars);
      candList[negated] = FALSE;
   }

   /* no candidates lifted so far */
   *nlifted = 0;
   bestcand = 0;
   liftround = 0;

   /* try lifting as long as we have lifting candidates */
   while( bestcand >= 0 )
   {
      /* in case we use a lifting rule, which does not require the first lifting coefficient of all variables: REMOVE this */
      if( sepadata->recalcliftcoef || liftround == 0 )
      {
         for( i = 0; i < 2*nbinvars; ++i )
         {
            if( candList[i] )
            {
               coef[i] = getCoef(scip, i, cycle, ncyclevars, vars, nbinvars, lifted, nlifted, sepadata, myi);
               assert(coef[i] <= (ncyclevars-1)/2);
               if( coef[i] < 1 )
                  candList[i] = FALSE;
            }
         }
      }
      ++liftround;
      bestcand = -1;
      for( i = 0; i < 2*nbinvars; ++i )
      {
         if( candList[i] )
         {
            /* we want to weight our choice of the lifting node by the value of the current lp solution */
            if( sepadata->lpliftcoef )
            {
               if( bestcand < 0 || coef[i]*vals[i] > coef[bestcand]*vals[bestcand] )
                  bestcand = (int) i;
            }
            /* we only regard the coefficient */
            else
            {
               if( bestcand < 0 || coef[i] > coef[bestcand] )
                  bestcand = (int) i;
            }
         }
      }

      /* there is at least one lifting variable */
      if( bestcand >= 0 )
      {
         if( !(sepadata->recalcliftcoef) )
            coef[i] = getCoef(scip, (unsigned int) bestcand, cycle, ncyclevars, vars, nbinvars, lifted, nlifted, sepadata, myi);
         assert(coef[bestcand] <= (ncyclevars-1)/2);
         candList[bestcand] = FALSE;
         if( coef[bestcand] > 0 )
         {
            if( bestcand >= (int)nbinvars )
               negated = bestcand - nbinvars;
            else
               negated = bestcand + nbinvars;
            assert(negated < 2*nbinvars);

            candList[negated] = FALSE;

            assert(*nlifted < nbinvars-ncyclevars);
            lifted[*nlifted] = (unsigned int) bestcand;
            liftcoef[*nlifted] = coef[bestcand];
            ++(*nlifted);
         }
      }
   }

 TERMINATE:
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &myi);
   SCIPfreeBufferArray(scip, &candList);
   SCIPfreeBufferArray(scip, &coef);
   SCIPfreeBufferArray(scip, &cycle);

   return SCIP_OKAY;
}

/*
 * methods for both techniques
 */

/* add the inequality corresponding to the given odd cycle to the LP (if violated)
 * after lifting it (if requested by user flag)
 */
static
SCIP_RETCODE generateOddCycleCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< given primal solution */
   SCIP_VAR**            vars,               /**< problem variables */
   unsigned int          nbinvars,           /**< number of binary problem variables */
   unsigned int          startnode,          /**< a node of the cycle */
   unsigned int*         pred,               /**< predecessor of each node (original and negated) in odd cycle */
   unsigned int          ncyclevars,         /**< number of variables in the odd cycle */
   SCIP_Bool*            incut,              /**< TRUE iff node is covered already by a cut */
   SCIP_Real*            vals,               /**< values of the variables in the given solution */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   unsigned int i;
   unsigned int negatedcount;
   unsigned int negated;

   /* memory for lifting */
   unsigned int nlifted;   /* number of lifted variables */
   unsigned int* lifted;   /* index of the lifted variables */
   unsigned int* liftcoef; /* lifting coefficient */

   /* memory for cut generation */
   SCIP_ROW* cut;
   char cutname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(vars != NULL);
   assert(startnode < 2*nbinvars);
   assert(pred != NULL);
   assert(ncyclevars % 2 == 1);
   assert(ncyclevars <= nbinvars);
   assert(incut != NULL);
   assert(sepadata != NULL);
   assert(sepadata->levelgraph != NULL || sepadata->useclassical);
   assert(sepadata->dijkstragraph != NULL || !(sepadata->useclassical));
   assert(result != NULL);

   /* debug method that prints out all found cycles */
   SCIPdebug(printCycle(vars,pred,nbinvars,startnode));

   /* cycle contains only one node */
   if( ncyclevars < 3 )
   {
      SCIPdebugMessage("fixing variable\n");
      /* strengthening variable bounds due to single-variable-cycle */
      if( startnode < nbinvars )
      {
         SCIP_CALL( SCIPchgVarUb(scip, vars[startnode], 0.0) );
      }
      else
      {
         negated = startnode - nbinvars;
         assert(negated < nbinvars);
         SCIP_CALL( SCIPchgVarLb(scip, vars[negated], 1.0) );
      }
      *result = SCIP_REDUCEDDOM;
      return SCIP_OKAY;
   }

   /* cycle is a triangle (can be excluded by user) */
   if( ncyclevars < 5 && !(sepadata->includetriangles) )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* lift the cycle inequality */
   nlifted = 0;
   lifted = NULL;
   liftcoef = NULL;
   if( sepadata->liftoddcycles )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &lifted, (int) (nbinvars - ncyclevars)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &liftcoef, (int) (nbinvars - ncyclevars)) );
      SCIP_CALL( liftOddCycleCut(scip, &nlifted, lifted, liftcoef, sepadata, vars, nbinvars, startnode, pred, ncyclevars, vals, result) );
   }
   /* if we don't try to lift, we generate and add the cut as is */

   /* create cut */
   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "oddcycle_%d", sepadata->ncuts);
   SCIP_CALL( SCIPcreateEmptyRow(scip, &cut, cutname, -SCIPinfinity(scip), (ncyclevars-1)/2.0, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
   negatedcount = 0;

   /* add variables of odd cycle to cut inequality */
   i = pred[startnode];
   while( i != startnode )
   {
      assert(i < 2*nbinvars);
      if( i < nbinvars )
      {
         /* inserting original variable */
         SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[i], 1.0) );
         incut[i] = TRUE;
      }
      else
      {
         negated = i - nbinvars;
         assert(negated < nbinvars);

         /* inserting negated variable */
         SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[negated], -1.0) );
         ++negatedcount;
         incut[negated] = TRUE;
      }
      i = pred[i];
   }

   /* insert startnode */
   if( startnode < nbinvars )
   {
      /* inserting original variable */
      SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[startnode], 1.0) );
      incut[i] = TRUE;
   }
   else
   {
      negated = startnode - nbinvars;
      assert(negated < nbinvars);

      /* inserting negated variable */
      SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[negated], -1.0) );
      ++negatedcount;
      incut[negated] = TRUE;
   }

   /* add lifted variables to cut inequality (if existing) */
   for( i = 0; i < nlifted; ++i)
   {
      assert(lifted != NULL);
      assert(liftcoef != NULL);
      if( lifted[i] < nbinvars )
      {
         assert(vars[lifted[i]] != NULL);
         SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[lifted[i]], (SCIP_Real) liftcoef[i]) );
      }
      else
      {
         negated = lifted[i] - nbinvars;
         assert(negated < nbinvars);
         assert(vars[negated] != NULL);
         SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[negated], -1.0*liftcoef[i]) );
         negatedcount += liftcoef[i];
      }
   }

   /* modify right hand side corresponding to number of added negated variables */
   SCIP_CALL( SCIPchgRowRhs(scip, cut, SCIProwGetRhs(cut)-negatedcount) );
   SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

   /* not every odd cycle has to be violated due to incompleteness of the implication graph */
   if( SCIPisCutEfficacious(scip, sol, cut) )
   {
      SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE) );
      SCIP_CALL( SCIPaddPoolCut(scip, cut) );
      ++(sepadata->ncuts);
      if( *result == SCIP_DIDNOTFIND )
         *result = SCIP_SEPARATED;

      assert(*result == SCIP_SEPARATED || *result == SCIP_REDUCEDDOM);
   }

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   if( sepadata->liftoddcycles )
   {
      SCIPfreeBufferArray(scip, &liftcoef);
      SCIPfreeBufferArray(scip, &lifted);
   }
   return SCIP_OKAY;
}

/** check whether the given object is really a cycle without sub-cycles (sub-cycles may be
 *  calculated by the GLS algorithm in case there is no violated odd cycle inequality) and removes
 *  pairs of original and negated variables from the cycle
 */
static
SCIP_RETCODE cleanCycle(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int*         pred,               /**< predecessor list of current cycle segment */
   SCIP_Bool*            incycle,            /**< flag array iff node is in cycle segment */
   SCIP_Bool*            incut,              /**< flag array iff node is already covered by a cut */
   unsigned int          x,                  /**< index of current variable */
   unsigned int          startnode,          /**< index of first variable of cycle segment */
   unsigned int          nbinvars,           /**< number of binary problem variables */
   unsigned int*         ncyclevars,         /**< number of nodes in current cycle segment */
   SCIP_Bool             repaircycles,       /**< user flag if we should try to repair damaged cycles */
   SCIP_Bool             allowmultiplecuts,  /**< user flag if we allow multiple cuts per node */
   SCIP_Bool*            success             /**< FALSE iff an irreparable cycle appears */
   )
{
   unsigned int negx;

   assert(scip != NULL);
   assert(pred != NULL);
   assert(incycle != NULL);
   assert(incut != NULL);
   assert(ncyclevars != NULL);
   assert(*ncyclevars <= nbinvars);
   assert(success != NULL);
   assert(*success);

   assert(x < 2*nbinvars);

   /* skip variable if it is already covered by a cut and we do not allow multiple cuts per node */
   if( incut[x] && !allowmultiplecuts )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* get index of negated variable of current variable */
   if( x < nbinvars )
      negx = x + nbinvars;
   else
      negx = x - nbinvars;
   assert(negx < 2*nbinvars);

   /* given object is not an odd cycle (contains sub-cycle) or contains original and negated
    * variable pair but we should not repair this
    */
   if( incycle[x] || (incycle[negx] && !repaircycles) )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* cycle does not contain original and negated variable pair */
   if( !incycle[negx] )
   {
      assert(!incycle[x]);
      incycle[x] = TRUE;
      ++(*ncyclevars);
      return SCIP_OKAY;
   }

   /* delete original and negated variable and cross-link their neighbors the following way, if possible:
    * suppose the cycle contains segments:
    * startnode - ... - a - neg(x) - c1 - c2 - ... - cn-1 - cn - x - z=pred(x)
    *
    * because of the chain a - neg(x) - x - cn it holds that
    *   a=1 => x=0 => neg(x)=1 => cn=0 and
    *   cn=1 => x=0 => neg(x)=1 => a=0
    * because of the chain z - x - neg(x) - b it holds that
    *   z=1 => x=0 => neg(x)=1 => c1=0 and
    *   c1=1 => x=0 => neg(x)=1 => z=0
    *
    * in addition to that, in our linked list structure we need to relink the chain c1-...-cn in reverse order.
    * so we gain the order: a - cn - cn-1 - ... - c2 - c1 - z
    */

   /* if negated variable is first node in cycle,
    * cross-linking not possible because there is no successor z of neg(x) contained in cycle yet
    */
   if( negx == startnode )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* if original and negated variable are neighbors, cross linking is not possible,
    * but x and neg(x) can simply be removed
    * a - neg(x)=pred[a] - x=pred[neg(x)] - z=pred[x] --> a - z=pred[x]=:pred[a]
    */
   if( pred[negx] == x )
   {
      unsigned int a;

      /* find a */
      a = startnode;
      while( pred[a] != negx )
         a = pred[a];

      /* link a and z */
      pred[a] = pred[x];
   }
   /* cross linking as mentioned above */
   else
   {
      unsigned int a;
      unsigned int z;

      /* memory for chain reverse */
      unsigned int* chain;
      unsigned int nchain;

      unsigned int i;

      /* allocate temporary memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &chain, (int) *ncyclevars) );

      /* find and store a */
      a = startnode;
      while( pred[a] != negx )
         a = pred[a];

      /* store chain */
      i = pred[negx];
      nchain = 0;
      while( i != x )
      {
         chain[nchain] = i;
         ++nchain;
         i = pred[i];
      }
      assert(nchain > 0);

      /* store z */
      z = pred[x];

      /* link a and c1 */
      pred[a] = chain[nchain-1];

      /* link cn and z */
      pred[chain[0]] = z;

      /* reverse the chain */
      for( i = nchain-1; i > 0; --i )
         pred[chain[i]] = chain[i-1];

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &chain);
   }

   /* remove negated variable from cycle */
   assert(!incycle[x] && incycle[negx]);
   incycle[negx] = FALSE;
   --(*ncyclevars);

   return SCIP_OKAY;
}

/*
 * methods for separateHeur()
 */

/** memory reallocation method (the graph is normally very dense, so we dynamically allocate only the memory we need)
 *
 *  Since the array sizes differ the method can be called for each of the three data structure types:
 *  - Forward: sizeForward, targetForward, weightForward
 *  - Backward: sizeBackward, targetBackward, weightBackward
 *  - Adj (inner level edges): sizeAdj, sourceAdj, targetAdj, weightAdj
 */
static
SCIP_RETCODE checkArraySizesHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   LEVELGRAPH*           graph,              /**< LEVELGRAPH data structure */
   unsigned int*         size,               /**< given size */
   int**                 targetArray,        /**< given target array (or NULL if sourceAdjArray and targetAdjArray given) */
   unsigned int**        weightArray,        /**< given weight array */
   unsigned int**        sourceAdjArray,     /**< given sourceAdj array (or NULL if targetArray given) */
   unsigned int**        targetAdjArray,     /**< given targetAdj array (or NULL if targetArray given) */
   SCIP_Bool*            success             /**< FALSE, iff memory reallocation fails */
   )
{
   SCIP_Real memorylimit;
   unsigned int additional;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(size != NULL);
   assert(targetArray != NULL || (sourceAdjArray != NULL && targetAdjArray != NULL));
   assert(weightArray != NULL);
   assert(success != NULL);

   SCIPdebugMessage("reallocating...\n");

   additional = MIN(graph->m + graph->n - *size, *size) * ((int) sizeof(**weightArray));
   if( targetArray != NULL )
   {
      additional += MIN(graph->m + graph->n - *size, *size) * ((int) sizeof(**targetArray));
   }
   else
   {
      additional += MIN(graph->m + graph->n - *size, *size) * ((int) sizeof(**sourceAdjArray));
      additional += MIN(graph->m + graph->n - *size, *size) * ((int) sizeof(**targetAdjArray));
   }

   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;

   /* if memorylimit would be exceeded or any other limit is reached free all data and exit */
   if( memorylimit <= additional/1048576.0 || SCIPisStopped(scip) )
   {
      *success = FALSE;
      SCIPdebugMessage("...memory limit exceeded\n");
      return SCIP_OKAY;
   }

   *size = 2 * (*size);

   SCIP_CALL( SCIPreallocBufferArray(scip, weightArray, (int) MIN(graph->m + graph->n, *size)) );
   if( targetArray != NULL )
   {
      SCIP_CALL( SCIPreallocBufferArray(scip, targetArray, (int) MIN(graph->m + graph->n, *size)) );
   }
   else
   {
      SCIP_CALL( SCIPreallocBufferArray(scip, sourceAdjArray, (int) MIN(graph->m, *size)) );
      SCIP_CALL( SCIPreallocBufferArray(scip, targetAdjArray, (int) MIN(graph->m, *size)) );
   }

   /* if memorylimit is exceeded free all data and exit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
   if( memorylimit <= 0.0 )
   {
      *success = FALSE;
      SCIPdebugMessage("...memory limit exceeded\n");
      return SCIP_OKAY;
   }

   SCIPdebugMessage("...with success\n");

   return SCIP_OKAY;
}

/** Add arc to level graph */
static
SCIP_RETCODE addArc(
   SCIP*                 scip,               /**< SCIP data structure */
   LEVELGRAPH*           graph,              /**< LEVELGRAPH data structure */
   unsigned int          u,                  /**< source node */
   unsigned int          v,                  /**< target node */
   unsigned int          level,              /**< number of current level */
   unsigned int          weight,             /**< weight of the arc */
   unsigned int*         nAdj,               /**< array of numbers of arcs inside levels */
   SCIP_Bool*            success             /**< FALSE, iff memory reallocation fails */
   )
{
   /* arc is a forward arc */
   if( graph->level[v] == level+1 )
   {
      graph->targetForward[graph->lastF] = (int) v;
      graph->weightForward[graph->lastF] = weight;
      ++(graph->lastF);
      ++(graph->nedges);
      if( graph->lastF == graph->sizeForward )
      {
         SCIP_CALL( checkArraySizesHeur(scip, graph, &(graph->sizeForward), &(graph->targetForward),
               &(graph->weightForward), NULL, NULL, success) );
         if( !(*success) )
            return SCIP_OKAY;
      }
   }
   else
   {
      assert(graph->level[v] == level || graph->level[v] == level-1);

      /* arc is a backward arc */
      if( graph->level[v] == level-1 )
      {
         graph->targetBackward[graph->lastB] = (int) v;
         graph->weightBackward[graph->lastB] = weight;
         ++(graph->lastB);
         ++(graph->nedges);

         if( graph->lastB == graph->sizeBackward )
         {
            SCIP_CALL( checkArraySizesHeur(scip, graph, &(graph->sizeBackward), &(graph->targetBackward),
                  &(graph->weightBackward), NULL, NULL, success) );
            if( !(*success) )
               return SCIP_OKAY;
         }
      }
      else       /* arc is in the same level */
      {
         assert(graph->level[v] == level);

         /* add arc only once, i.e., if u < v */
         if( u < v )
         {
            graph->sourceAdj[graph->levelAdj[level+1]+*nAdj] = u;
            graph->targetAdj[graph->levelAdj[level+1]+*nAdj] = v;
            graph->weightAdj[graph->levelAdj[level+1]+*nAdj] = weight;
            ++(*nAdj);
            ++(graph->nedges);

            if( graph->levelAdj[level+1]+*nAdj == graph->sizeAdj )
            {
               SCIP_CALL( checkArraySizesHeur(scip, graph, &(graph->sizeAdj), NULL, &(graph->weightAdj),
                     &(graph->sourceAdj), &(graph->targetAdj), success) );
               if( !(*success) )
                  return SCIP_OKAY;
            }
         }
      }
   }
   return SCIP_OKAY;
}

/** add binary implications of the given node u @see createNextLevel() */
static
SCIP_RETCODE addNextLevelBinImpls(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   SCIP_Real*            vals,               /**< values of the binary variables in the current LP relaxation */
   unsigned int          u,                  /**< current node */
   LEVELGRAPH*           graph,              /**< LEVELGRAPH data structure */
   unsigned int          level,              /**< number of current level */
   SCIP_Bool*            inlevelgraph,       /**< flag array if node is already inserted in level graph */
   unsigned int*         newlevel,           /**< array of nodes of the next level */
   unsigned int*         nnewlevel,          /**< number of nodes of the next level */
   unsigned int*         nAdj,               /**< array of numbers of arcs inside levels */
   SCIP_Bool*            success             /**< FALSE, iff memory reallocation fails */
   )
{
   SCIP_Bool varfixing;
   unsigned int nbinimpls;
   unsigned int nbinvars;
   unsigned int varsidx;
   SCIP_VAR** implvars;
   SCIP_BOUNDTYPE* impltypes;
   unsigned int j;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(graph != NULL);
   assert(graph->targetForward != NULL);
   assert(graph->weightForward != NULL);
   assert(graph->targetBackward != NULL);
   assert(graph->weightBackward != NULL);
   assert(graph->sourceAdj != NULL);
   assert(graph->targetAdj != NULL);
   assert(graph->weightAdj != NULL);
   assert(inlevelgraph != NULL);
   assert(newlevel != NULL);
   assert(nnewlevel != NULL);
   assert(nAdj != NULL);
   assert(success != NULL);

   assert(u < (graph->n));

   nbinvars = (graph->n)/2;

   /* current node signifies a problem variable */
   if( u < nbinvars )
   {
      varfixing = TRUE;
      varsidx = u;
   }
   /* current node signifies a negated variable */
   else
   {
      varfixing = FALSE;
      varsidx = u - nbinvars;
   }
   assert(varsidx < nbinvars);
   assert(!SCIPisFeasIntegral(scip, vals[varsidx]));

   /* get binary implications of the current variable */
   nbinimpls = (unsigned int) SCIPvarGetNBinImpls(vars[varsidx], varfixing);
   implvars = SCIPvarGetImplVars(vars[varsidx], varfixing);
   impltypes = SCIPvarGetImplTypes(vars[varsidx], varfixing);

   assert(implvars != NULL || nbinimpls == 0);
   assert(impltypes != NULL || nbinimpls == 0);

   for( j = 0; j < nbinimpls; ++j )
   {
      unsigned int k;
      unsigned int v;
      unsigned int weight;

      assert( implvars != NULL && impltypes != NULL );
      assert(SCIPvarGetType(implvars[j]) == SCIP_VARTYPE_BINARY);

      k = sepadata->mapping[SCIPvarGetProbindex(implvars[j])];
      assert(k < nbinvars);

      /* skip integral neighbors */
      if( SCIPisFeasIntegral(scip, vals[k]) )
         continue;

      /* consider implication to negated variable (x = 1 -> y >= 1 <=>  x = 1 -> neg(y) <= 0) */
      if( impltypes[j] == SCIP_BOUNDTYPE_LOWER )
         v = k + nbinvars;
      /* x = 1 -> y <= 0 */
      else
      {
         assert(impltypes[j] == SCIP_BOUNDTYPE_UPPER);
         v = k;
      }
      assert(v < (graph->n));

      /* if variable is a new node, it will be assigned to the next level, but if the level contains
       * more nodes than allowed (defined by percent per level plus offset), we skip the rest of the
       * nodes
       */
      if( !inlevelgraph[v] && (*nnewlevel) <= sepadata->maxlevelsize )
      {
         ++(graph->nnodes);
         graph->level[v] = level+1;
         inlevelgraph[v] = TRUE;
         newlevel[*nnewlevel] = v;
         ++(*nnewlevel);
      }
      assert((*nnewlevel) > sepadata->maxlevelsize || inlevelgraph[v]);

      /* calculate arc weight and add arc, if the neighbor node is on the same or a neighbor level */
      if( inlevelgraph[v] && (graph->level[v] == level+1 || graph->level[v] == level || graph->level[v] == level-1))
      {
         SCIP_Real tmp;

         /* set weight of arc (x,y) to 1 - x* -y* */
         if( varfixing )
         {
            /* x = 1 -> y <= 0 */
            if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
            {
               tmp = SCIPfeasCeil(scip, sepadata->scale * (1.0 - vals[varsidx] - vals[k]));
               weight = (unsigned int) MAX(tmp, sepadata->maxreference);
            }
            /* x = 1 -> y >= 1 <-> neg(y) <= 0 */
            else
            {
               tmp = SCIPfeasCeil(scip, sepadata->scale * (1.0 - vals[varsidx] - (1-vals[k])));
               weight = (unsigned int) MAX(tmp, sepadata->maxreference);
            }
         }
         else
         {
            /* x = 0 <-> neg(x) = 1 -> y <= 0 */
            if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
            {
               tmp = SCIPfeasCeil(scip, sepadata->scale * (1.0 - (1 - vals[varsidx]) - vals[k]));
               weight = (unsigned int) MAX(tmp, sepadata->maxreference);
            }
            /* x = 0 <-> neg(x) = 1 -> y >= 1 <-> neg(y) <= 0 */
            else
            {
               tmp = SCIPfeasCeil(scip, sepadata->scale * (1.0 - (1 - vals[varsidx]) - (1-vals[k])));
               weight = (unsigned int) MAX(tmp, sepadata->maxreference);
            }
         }

         /* add arc from current to neighbor node */
         SCIP_CALL( addArc(scip, graph, u, v, level, weight, nAdj, success) );
         if( !(*success) )
            return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}

/** add implications from cliques of the given node u
 *
 *  @see createNextLevel()
 */
static
SCIP_RETCODE addNextLevelCliques(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   SCIP_Real*            vals,               /**< values of the binary variables in the current LP relaxation */
   unsigned int          u,                  /**< current node */
   LEVELGRAPH*           graph,              /**< LEVELGRAPH data structure */
   unsigned int          level,              /**< number of current level */
   SCIP_Bool*            inlevelgraph,       /**< flag array if node is already inserted in level graph */
   unsigned int*         newlevel,           /**< array of nodes of the next level */
   unsigned int*         nnewlevel,          /**< number of nodes of the next level */
   unsigned int*         nAdj,               /**< array of numbers of arcs inside levels */
   SCIP_Bool*            success             /**< FALSE, iff memory reallocation fails */
   )
{
   SCIP_Bool varfixing;
   unsigned int ncliques;
   unsigned int nbinvars;
   unsigned int varsidx;
   SCIP_CLIQUE** cliques;
   unsigned int ncliquevars;
   SCIP_VAR** cliquevars;
   SCIP_Bool* cliquevals;
   unsigned int j;
   unsigned int k;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(graph != NULL);
   assert(graph->targetForward != NULL);
   assert(graph->weightForward != NULL);
   assert(graph->targetBackward != NULL);
   assert(graph->weightBackward != NULL);
   assert(graph->sourceAdj != NULL);
   assert(graph->targetAdj != NULL);
   assert(graph->weightAdj != NULL);
   assert(inlevelgraph != NULL);
   assert(newlevel != NULL);
   assert(nnewlevel != NULL);
   assert(nAdj != NULL);
   assert(success != NULL);

   assert(u < (graph->n));

   nbinvars = (graph->n)/2;

   /* current node signifies a problem variable */
   if( u < nbinvars )
   {
      varfixing = TRUE;
      varsidx = u;
   }
   /* current node signifies a negated variable */
   else
   {
      varfixing = FALSE;
      varsidx = u - nbinvars;
   }
   assert(varsidx < nbinvars);
   assert(!SCIPisFeasIntegral(scip, vals[varsidx]));

   /* get cliques of the current variable */
   ncliques = (unsigned int) SCIPvarGetNCliques(vars[varsidx], varfixing);
   if( ncliques == 0 )
      return SCIP_OKAY;

   cliques = SCIPvarGetCliques(vars[varsidx], varfixing);
   assert(cliques != NULL);

   for( j = 0; j < ncliques; ++j )
   {
      ncliquevars = (unsigned int) SCIPcliqueGetNVars(cliques[j]);
      cliquevars = SCIPcliqueGetVars(cliques[j]);
      cliquevals = SCIPcliqueGetValues(cliques[j]);

      assert(cliquevars != NULL || ncliquevars == 0);
      assert(cliquevals != NULL || ncliquevars == 0);

      for( k = 0; k < ncliquevars; ++k )
      {
         unsigned int l;
         unsigned int v;
         unsigned int weight;

         assert( cliquevars != NULL && cliquevals != NULL );  /* for lint */

         l = sepadata->mapping[SCIPvarGetProbindex(cliquevars[k])];
         assert(l < nbinvars);

         /* skip integral neighbors */
         if( SCIPisFeasIntegral(scip, vals[l]) )
            continue;

         /* consider clique with negated variable (x = 1 -> y >= 1 <=>  x = 1 -> neg(y) <= 0) */
         if( cliquevals[k] == FALSE )
            v = l + nbinvars;
         /* x = 1 -> y <= 0 */
         else
            v = l;
         assert(v < (graph->n));

         /* if variable is a new node, it will be assigned to the next level,
          * but if the level contains more nodes than allowed
          * (defined by percent per level plus offset),
          * we skip the rest of the nodes
          */
         if( !inlevelgraph[v] && (*nnewlevel) <= sepadata->maxlevelsize )
         {
            ++(graph->nnodes);
            graph->level[v] = level+1;
            inlevelgraph[v] = TRUE;
            newlevel[*nnewlevel] = v;
            ++(*nnewlevel);
         }
         assert(*nnewlevel > sepadata->maxlevelsize || inlevelgraph[v]);

         /* calculate arc weight and add arc, if the neighbor node is on the same or a neighbor level */
         if( inlevelgraph[v] && (graph->level[v] == level+1 || graph->level[v] == level || graph->level[v] == level-1))
         {
            SCIP_Real tmp;

            /* set weight of arc (x,y) to 1 - x* -y* */
            if( varfixing )
            {
               /* x = 1 -> y <= 0 */
               if( cliquevals[k] )
               {
                  tmp = SCIPfeasCeil(scip, sepadata->scale * (1 - vals[varsidx] - vals[l]));
                  weight = (unsigned int) MAX(tmp, sepadata->maxreference);
               }
               /* x = 1 -> y >= 1 <-> neg(y) <= 0 */
               else
               {
                  tmp = SCIPfeasCeil(scip, sepadata->scale * (1 - vals[varsidx] - (1-vals[l])));
                  weight = (unsigned int) MAX(tmp, sepadata->maxreference);
               }
            }
            else
            {
               /* x = 0 <-> neg(x) = 1 -> y <= 0 */
               if( !cliquevals[k] )
               {
                  tmp = SCIPfeasCeil(scip, sepadata->scale * (1 - (1 - vals[varsidx]) - vals[l]));
                  weight = (unsigned int) MAX(tmp, sepadata->maxreference);
               }
               /* x = 0 <-> neg(x) = 1 -> y >= 1 <-> neg(y) <= 0 */
               else
               {
                  tmp = SCIPfeasCeil(scip, sepadata->scale * (1 - (1 - vals[varsidx]) - (1-vals[l])));
                  weight = (unsigned int) MAX(tmp, sepadata->maxreference);
               }
            }

            /* add arc from current to neighbor node */
            SCIP_CALL( addArc(scip, graph, u, v, level, weight, nAdj, success) );
            if( !(*success) )
               return SCIP_OKAY;
         }
      }
   }
   return SCIP_OKAY;
}


/** sort level of root neighbors
 *
 *  If we limit the size of nodes of a level, we want to add the best neighbors to the next level.
 *  Since sorting every level is too expensive, we sort the neighbors of the root (if requested).
 *
 *  Create the first level as follows:
 *  - create flag array for binary variables and their negated and set their values FALSE
 *  - iterate over the implication and clique neighbors of the root and set their flag array values to TRUE
 *  - create variable array and insert all variables with flag value TRUE
 *  - sort variable array by maximal fractionality
 *  - add variables from sorted array to levelgraph until first level is full (or all variables are inserted)
 *
 *  Even inserting all variables might help for the following creation of further levels since the neighbors
 *  of nodes with high fractionality often have high fractionalities themselves and would be inserted first
 *  when further levels would have been sorted (which actually is not the case).
 */
static
SCIP_RETCODE insertSortedRootNeighbors(
   SCIP*                 scip,               /**< SCIP data structure */
   LEVELGRAPH*           graph,              /**< LEVELGRAPH data structure */
   unsigned int          nbinvars,           /**< number of binary problem variables */
   unsigned int          ncurlevel,          /**< number of nodes of the current level */
   unsigned int          u,                  /**< source node */
   SCIP_Real*            vals,               /**< values of the binary variables in the current LP relaxation */
   SCIP_VAR**            vars,               /**< problem variables */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   unsigned int*         nnewlevel,          /**< number of nodes of the next level */
   SCIP_Bool*            inlevelgraph,       /**< nodes in new graph corr. to old graph (-1 if unassigned) */
   unsigned int          level,              /**< number of current level */
   unsigned int*         newlevel,           /**< array of nodes of the next level */
   SCIP_Bool*            success             /**< FALSE, iff memory reallocation fails */
   )
{
   /* storage for the neighbors of the root */
   unsigned int root;
   unsigned int nneighbors;
   SCIP_Bool* isneighbor;
   int* neighbors;
   SCIP_Real* neighvals;
   SCIP_Real* sortvals;

   SCIP_Bool varfixing;
   unsigned int varsidx;

   /* storage for implications to the neighbors of the root node */
   unsigned int nbinimpls;
   SCIP_VAR** implvars;
   SCIP_BOUNDTYPE* impltypes;

   /* storage for cliques to the neighbors of the root node */
   unsigned int ncliques;
   SCIP_CLIQUE** cliques;
   unsigned int ncliquevars;
   SCIP_VAR** cliquevars;
   SCIP_Bool* cliquevals;

   unsigned int j;
   unsigned int k;
   unsigned int v;

   /* allocate flag array for neighbor detection */
   SCIP_CALL( SCIPallocBufferArray(scip, &isneighbor, (int) graph->n) );
   BMSclearMemoryArray(isneighbor, graph->n);

   nbinvars = (graph->n)/2;

   assert(ncurlevel == 1);
   root = u;

   /* current node signifies a problem variable */
   if( root < nbinvars )
   {
      varfixing = TRUE;
      varsidx = root;
   }
   /* current node signifies a negated variable */
   else
   {
      varfixing = FALSE;
      varsidx = root - nbinvars;
   }
   assert(varsidx < nbinvars);
   assert(!SCIPisFeasIntegral(scip, vals[varsidx]));

   /* count implications of the root */
   nbinimpls = (unsigned int) SCIPvarGetNBinImpls(vars[varsidx], varfixing);
   if( nbinimpls > 0 )
   {
      unsigned int jidx;

      implvars = SCIPvarGetImplVars(vars[varsidx], varfixing);
      impltypes = SCIPvarGetImplTypes(vars[varsidx], varfixing);
      assert(implvars != NULL);
      assert(impltypes != NULL);

      for( j = 0; j < nbinimpls ; ++j )
      {
         jidx = sepadata->mapping[SCIPvarGetProbindex(implvars[j])];
         assert(jidx < nbinvars);

         if( SCIPisFeasIntegral(scip, vals[jidx]))
            continue;

         if( varfixing == TRUE)
         {
            /* implication x + y <= 1 */
            if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
               isneighbor[jidx] = TRUE;
            /* implication x + neg(y) <= 1 */
            else
            {
               assert(impltypes[j] == SCIP_BOUNDTYPE_LOWER);
               isneighbor[jidx+nbinvars] = TRUE;
            }
         }
         else
         {
            /* implication neg(x) + neg(y) <= 1 */
            if( impltypes[j] == SCIP_BOUNDTYPE_LOWER )
               isneighbor[jidx+nbinvars] = TRUE;
            /* implication neg(x) + y <= 1 */
            else
            {
               assert(impltypes[j] == SCIP_BOUNDTYPE_UPPER);
               isneighbor[jidx] = TRUE;
            }
         }
      }
   }

   /* count cliques of the root */
   ncliques = (unsigned int) SCIPvarGetNCliques(vars[varsidx], varfixing);
   if( ncliques > 0 )
   {
      cliques = SCIPvarGetCliques(vars[varsidx], varfixing);
      assert(cliques != NULL);

      for( j = 0; j < ncliques; ++j )
      {
         ncliquevars = (unsigned int) SCIPcliqueGetNVars(cliques[j]);
         cliquevars = SCIPcliqueGetVars(cliques[j]);
         cliquevals = SCIPcliqueGetValues(cliques[j]);

         assert(cliquevars != NULL || ncliquevars == 0);
         assert(cliquevals != NULL || ncliquevars == 0);

         for( k = 0; k < ncliquevars; ++k )
         {
            unsigned int kidx;

            assert( cliquevars != NULL && cliquevals != NULL ); /* for lint */

            kidx = sepadata->mapping[SCIPvarGetProbindex(cliquevars[k])];
            assert(kidx < nbinvars);

            /* skip integral neighbors */
            if( SCIPisFeasIntegral(scip, vals[kidx]))
               continue;

            /* skip root */
            if( kidx == varsidx )
               continue;

            if( cliquevals[k] == TRUE )
               isneighbor[kidx] = TRUE;
            else
            {
               assert(cliquevals[k] == FALSE);
               isneighbor[kidx+nbinvars] = TRUE;
            }
         }
      }
   }

   /* root cannot be part of the next level */
   assert(!isneighbor[root]);

   nneighbors = 0;
   /* count root neighbors */
   for( j = 0; j < graph->n; ++j )
   {
      if( isneighbor[j] )
         ++nneighbors;
   }

   /* allocate memory for sorting of root neighbors */
   SCIP_CALL( SCIPallocBufferArray(scip, &neighbors, (int) nneighbors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &neighvals, (int) nneighbors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortvals, (int) nneighbors) );

   k = 0;
   for( j = 0; j < graph->n; ++j )
   {
      if( isneighbor[j] )
      {
         assert(j != root);
         neighbors[k] = (int) j;
         neighvals[k] = vals[j];
         assert(!SCIPisFeasIntegral(scip, neighvals[k]));
         ++k;
      }
   }
   assert(k == nneighbors);

   /* calculate fractionality of neighbors */
   for( j = 0; j < nneighbors; ++j )
      sortvals[j] = MIN(1-neighvals[j],neighvals[j]);

   /* sort neighbors by fractionality */
   SCIPsortDownRealInt(sortvals, neighbors, (int) nneighbors);

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &sortvals);

   /* insert sorted neighbors until level size limit is reached (or all neighbors are inserted) */
   for( j = 0; j < nneighbors && (*nnewlevel) <= sepadata->maxlevelsize; ++j )
   {
      SCIP_Real tmp;

      v = (unsigned int) neighbors[j];

      /* only the root is contained in the levelgraph */
      assert(!inlevelgraph[v] || v == root+nbinvars || v == root-nbinvars);

      /* insert neighbor into levelgraph */
      ++(graph->nnodes);
      graph->level[v] = level+1;
      inlevelgraph[v] = TRUE;
      newlevel[*nnewlevel] = v;
      ++(*nnewlevel);

      assert(!SCIPisFeasIntegral(scip, vals[varsidx]));
      assert(!SCIPisFeasIntegral(scip, neighvals[j]));

      graph->targetForward[graph->lastF] = (int) v;
      if( varfixing )
      {
         tmp = SCIPfeasCeil(scip, sepadata->scale * (1.0 - vals[varsidx] - neighvals[j]));
         graph->weightForward[graph->lastF] = (unsigned int) MAX(tmp, sepadata->maxreference);
      }
      else
      {
         assert(varfixing == FALSE);
         tmp = SCIPfeasCeil(scip, sepadata->scale * (1.0 - (1.0-vals[varsidx]) - neighvals[j]));
         graph->weightForward[graph->lastF] = (unsigned int) MAX(tmp, sepadata->maxreference);
      }
      ++(graph->lastF);
      ++(graph->nedges);
      if( graph->lastF == graph->sizeForward )
      {
         SCIP_CALL( checkArraySizesHeur(scip, graph, &(graph->sizeForward), &(graph->targetForward),
               &(graph->weightForward), NULL, NULL, success) );

         if( !(*success) )
            break;
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &neighvals);
   SCIPfreeBufferArray(scip, &neighbors);
   SCIPfreeBufferArray(scip, &isneighbor);

   return SCIP_OKAY;
}

/** Find shortest path from start node to root
 *
 *  We perform a BFS to find the shortest path to the root. D stores the distance to the start
 *  node, P stores the partent nodes in the shortest path tree (-1 if node has not been reached).
 */
static
SCIP_RETCODE findShortestPathToRoot(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   scale,              /**< scaling factor for edge weights */
   LEVELGRAPH*           graph,              /**< LEVELGRAPH data structure */
   unsigned int          startnode,          /**< start node for path search */
   unsigned int*         distance,           /**< distances of searched nodes from root */
   unsigned int*         queue,              /**< node queue for path search */
   SCIP_Bool*            inQueue,            /**< whether node is in the queue */
   int*                  parentTree          /**< parent tree (-1 if no parent) */
   )
{
   unsigned int i;
   int startQueue;
   int endQueue;
   unsigned int u;
   int v;
   unsigned int d;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->beginBackward != NULL);
   assert(graph->targetBackward != NULL);
   assert(graph->weightBackward != NULL);
   assert(distance != NULL);
   assert(queue != NULL);
   assert(inQueue != NULL);
   assert(parentTree != NULL);

   /* initialize distances */
   for( i = 0; i < graph->n; ++i )
   {
      distance[i] = 2*(graph->nnodes)*scale;
      parentTree[i] = -1;
      inQueue[i] = FALSE;
   }
   distance[startnode] = 0;

   /* initialize queue */
   startQueue = 0;
   endQueue = 0;
   queue[0] = startnode;
   v = 0;

   /* as long as queue is not empty */
   while( startQueue <= endQueue )
   {
      /* pop first node from queue */
      u = queue[startQueue];
      ++startQueue;

      /* check adjacent nodes */
      assert(graph->beginBackward[u] >= 0);
      i = (unsigned int) graph->beginBackward[u];
      for( v = graph->targetBackward[i]; v >= 0; v = graph->targetBackward[++i] )
      {
         /* distance to u via current arc: */
         d = distance[u] + graph->weightBackward[i];

         /* if we found a shorter connection */
         if( d < distance[v] )
         {
            distance[v] = d;
            parentTree[v] = (int) u;

            /* insert in queue if not already present */
            if( !inQueue[v] )
            {
               ++endQueue;
               queue[endQueue] = (unsigned int) v;
               inQueue[v] = TRUE;
            }
         }
      }
      /* it is not necessary to stop if we found the root (in this case there are no arcs left) and we stop anyway */
   }
   assert(parentTree[u] != -1);

   return SCIP_OKAY;
}


/** Block shortest path
 *
 *  We traverse the shortest path found by findShortestPathToRoot() and block all neighbors (in the
 *  original graph) of nodes in the path, i.e., we set blocked to TRUE. We do not block neighbors of
 *  the root node, since they have to be used. For the start node we only block nodes on the
 *  previous layers,
 *
 *  @see findShortestPathToRoot()
 */
static
SCIP_RETCODE blockRootPath(
   SCIP*                 scip,               /**< SCIP data structure */
   LEVELGRAPH*           graph,              /**< LEVELGRAPH data structure */
   unsigned int          startnode,          /**< start node */
   SCIP_Bool*            inlevelgraph,       /**< nodes in new graph corr. to old graph (-1 if unassigned) */
   SCIP_Bool*            blocked,            /**< whether node is blocked */
   int*                  parentTree,         /**< parent tree */
   unsigned int          root                /**< root of the current level graph */
   )
{
   unsigned int u;
   unsigned int i;
   int v;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->level != NULL);
   assert(graph->beginForward != NULL);
   assert(graph->targetForward != NULL);
   assert(graph->beginBackward != NULL);
   assert(graph->targetBackward != NULL);
   assert(graph->sourceAdj != NULL);
   assert(graph->targetAdj != NULL);
   assert(inlevelgraph != NULL);
   assert(blocked != NULL);
   assert(parentTree != NULL);

   assert(parentTree[root] >= 0);

   /* follow the path from the predecessor of root to the start node and block all neighbors */
   u = (unsigned int) parentTree[root];
   while( u != startnode )
   {
      /*  block neighbors of u in higher level */
      i = (unsigned int) graph->beginForward[u];
      for( v = graph->targetForward[i]; v >= 0; v = graph->targetForward[++i] )
      {
         assert(inlevelgraph[v]);
         blocked[v] = TRUE;
      }

      /*  block neighbors of u in lower level */
      i = (unsigned int) graph->beginBackward[u];
      for( v = graph->targetBackward[i]; v >= 0; v = graph->targetBackward[++i] )
      {
         assert(inlevelgraph[v]);
         blocked[v] = TRUE;
      }

      /*  block neighbors of u in same level */
      assert(graph->level[u] > 0);
      for( i = graph->levelAdj[graph->level[u]]; i < graph->levelAdj[graph->level[u]+1]; ++i )
      {
         assert(graph->sourceAdj[i] < graph->targetAdj[i]);
         assert(graph->level[graph->sourceAdj[i]] == graph->level[graph->targetAdj[i]]);

         /* remember that these arcs are only stored for one direction */
         if( graph->sourceAdj[i] == u )
         {
            blocked[graph->targetAdj[i]] = TRUE;
         }
         if( graph->targetAdj[i] == u )
         {
            blocked[graph->sourceAdj[i]] = TRUE;
         }
      }

      /* get next node on the path */
      u = (unsigned int) parentTree[u];
   }
   assert(u == startnode);

   /* block nodes adjacent to start node on previous level */
   assert(graph->beginBackward[u] > 0);
   i = (unsigned int) graph->beginBackward[u];
   for( v = graph->targetBackward[i]; v >= 0; v = graph->targetBackward[++i] )
      blocked[v] = TRUE;

   return SCIP_OKAY;
}


/** Find shortest path from root to target node
 *
 *  We perform a BFS to find the shortest path from the root. The only difference to
 *  findShortestPathToRoot() is that we avoid nodes that are blocked.
 */
static
SCIP_RETCODE
findUnblockedShortestPathToRoot(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   scale,              /**< scaling factor for edge weights */
   LEVELGRAPH*           graph,              /**< LEVELGRAPH data structure */
   unsigned int          startnode,          /**< start node for path search */
   unsigned int*         distance,           /**< distances of searched nodes from root */
   unsigned int*         queue,              /**< node queue for path search */
   SCIP_Bool*            inQueue,            /**< whether node is in the queue */
   int*                  parentTreeBackward, /**< parent tree (-1 if no parent) */
   unsigned int          root,               /**< root of the current level graph */
   SCIP_Bool*            blocked             /**< whether nodes can be used */
   )
{
   unsigned int i;
   int startQueue;
   int endQueue;
   unsigned int u;
   int v;
   unsigned int d;
   int* parentTree;
   int* transform;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(graph->beginBackward != NULL);
   assert(graph->targetBackward != NULL);
   assert(graph->weightBackward != NULL);
   assert(distance != NULL);
   assert(queue != NULL);
   assert(inQueue != NULL);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &parentTree, (int) graph->n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transform, (int) graph->n) );

   assert(parentTree != NULL);
   assert(transform != NULL);

   /* initialize distances */
   for( i = 0; i < graph->n; ++i )
   {
      distance[i] = 2*(graph->nnodes)*scale;
      parentTree[i] = -1;
      parentTreeBackward[i] = -1;
      transform[i] = -1;
      inQueue[i] = FALSE;
   }
   distance[startnode] = 0;

   /* initialize queue */
   startQueue = 0;
   endQueue = 0;
   queue[0] = startnode;
   v = 0;

   /* as long as queue is not empty */
   while( startQueue <= endQueue )
   {
      /* pop first node from queue */
      u = queue[startQueue];
      ++startQueue;

      /* check adjacent nodes */
      assert(graph->beginBackward[u] >= 0);
      i = (unsigned int) graph->beginBackward[u];
      for( v = graph->targetBackward[i]; v >= 0; v = graph->targetBackward[++i] )
      {
         if( blocked[v] && v != (int) root)
            continue;

         /* distance to u via current arc: */
         d = distance[u] + graph->weightBackward[i];

         /* if we found a shorter connection */
         if( d < distance[v] )
         {
            distance[v] = d;
            parentTree[v] = (int) u;

            /* insert in queue if not already present */
            if( !inQueue[v] )
            {
               ++endQueue;
               queue[endQueue] = (unsigned int) v;
               inQueue[v] = TRUE;
            }
         }
      }
      /* it is not necessary to stop if we found the root (in this case there are no arcs left) and we stop anyway */
   }

   /* reverse order such that it is a path from the root */
   v = (int) root;
   transform[0] = (int) root;
   i = 1;
   while(parentTree[v] >= 0)
   {
      transform[i] = parentTree[v];
      ++i;
      v = parentTree[v];
   }
   --i;
   while(i > 0)
   {
      parentTreeBackward[transform[i]] = transform[i-1];
      --i;
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &transform);
   SCIPfreeBufferArray(scip, &parentTree);

   return SCIP_OKAY;
}

/** create next level of level graph for odd cycle separation
 *
 *  @see separateHeur()
 */
static
SCIP_RETCODE createNextLevel(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   SCIP_Real*            vals,               /**< values of the binary variables in the current LP relaxation */
   LEVELGRAPH*           graph,              /**< LEVELGRAPH data structure */
   unsigned int          level,              /**< number of current level */
   SCIP_Bool*            inlevelgraph,       /**< flag array if node is already inserted in level graph */
   unsigned int*         curlevel,           /**< array of nodes of the current level */
   unsigned int          ncurlevel,          /**< number of nodes of the current level */
   unsigned int*         newlevel,           /**< array of nodes of the next level */
   unsigned int*         nnewlevel,          /**< number of nodes of the next level */
   SCIP_Bool*            success             /**< FALSE, iff memory reallocation fails */
   )
{
   unsigned int i;
   unsigned int nbinvars;
   unsigned int nAdj;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(vals != NULL);
   assert(graph != NULL);
   assert(graph->level != NULL);
   assert(graph->beginForward != NULL);
   assert(graph->targetForward != NULL);
   assert(graph->weightForward != NULL);
   assert(graph->beginBackward != NULL);
   assert(graph->targetBackward != NULL);
   assert(graph->weightBackward != NULL);
   assert(graph->beginAdj != NULL);
   assert(graph->levelAdj != NULL);
   assert(graph->sourceAdj != NULL);
   assert(graph->targetAdj != NULL);
   assert(graph->weightAdj != NULL);
   assert(inlevelgraph != NULL);
   assert(curlevel != NULL);
   assert(newlevel != NULL);
   assert(success != NULL);

   *nnewlevel = 0;
   nAdj = 0;
   assert((graph->n) % 2 == 0);
   nbinvars = (graph->n)/2;

   /* for every node in current level add its implications and assign its neighbors to the next
    * level, if neighbor is not already existing in the level graph
    */
   for( i = 0; i < ncurlevel; ++i )
   {
      unsigned int negated;
      unsigned int u;

      /* get node */
      u = curlevel[i];
      assert(u < (graph->n));
      assert(graph->level[u] == level);
      assert(graph->beginForward[u] < 0);
      assert(graph->beginBackward[u] < 0);
      assert(graph->beginAdj[u] < 0);
      assert(inlevelgraph[u]);

      /* get negated */
      if( u < nbinvars )
         negated = u + nbinvars;
      else
         negated = u - nbinvars;
      assert(negated < (graph->n));
      assert(negated < nbinvars || u < nbinvars);
      assert(negated >= nbinvars || u >= nbinvars);

      /* initialize adjacency lists for node u */
      graph->beginForward[u] = (int) graph->lastF;
      graph->beginBackward[u] = (int) graph->lastB;
      graph->beginAdj[u] = (int) (graph->levelAdj[level+1] + nAdj);

      /* if we want to add arcs between a variable and its negated */
      if( sepadata->addselfarcs )
      {
         /* add negated variable, if not existing in the levelgraph,
          * but if the level contains more nodes than allowed
          * (defined by percent per level plus offset),
          * we skip the rest of the nodes
          */
         if( !inlevelgraph[negated] && (*nnewlevel) <= sepadata->maxlevelsize )
         {
            ++(graph->nnodes);
            graph->level[negated] = level+1;
            inlevelgraph[negated] = TRUE;
            newlevel[*nnewlevel] = negated;
            ++(*nnewlevel);
         }
         assert( *nnewlevel > sepadata->maxlevelsize || inlevelgraph[negated] );

         /* add self-arc if negated variable is on a neighbored level */
         if( inlevelgraph[negated] && ((graph->level[negated] == level - 1)
               || (graph->level[negated] == level) || (graph->level[negated] == level + 1)) )
         {
            /* add arc from u to its negated variable */
            SCIP_CALL( addArc(scip, graph, u, negated, level, 0, &nAdj, success) );
            if( !(*success) )
               return SCIP_OKAY;
         }
      }

      /* insert level of sorted root neighbors (if requested) */
      if( graph->nlevels == 0 && sepadata->sortrootneighbors )
      {
         SCIP_CALL( insertSortedRootNeighbors(scip, graph, nbinvars, ncurlevel, u, vals, vars,
               sepadata, nnewlevel, inlevelgraph, level, newlevel, success) );
      }
      else
      {
         /* add arc from u to all other neighbors of variable implication graph */
         SCIP_CALL( addNextLevelBinImpls(scip, sepadata, vars, vals, u, graph, level, inlevelgraph,
               newlevel, nnewlevel, &nAdj, success) );

         if( !(*success) )
            return SCIP_OKAY;

         SCIP_CALL( addNextLevelCliques(scip, sepadata, vars, vals, u, graph, level, inlevelgraph,
               newlevel, nnewlevel, &nAdj, success) );
      }
      if( !(*success) )
         return SCIP_OKAY;

      /* every node has a backward arc */
      assert(graph->lastB > (unsigned int) graph->beginBackward[u] || graph->nlevels == 0 );

      /* root has outgoing arcs otherwise we would have skipped it */
      assert(graph->lastF > 0);

      /* close adjacency lists */
      graph->targetForward[graph->lastF] = -1;
      ++(graph->lastF);
      if( graph->lastF == graph->sizeForward )
      {
         SCIP_CALL( checkArraySizesHeur(scip, graph, &(graph->sizeForward), &(graph->targetForward),
               &(graph->weightForward), NULL, NULL, success) );

         if( !(*success) )
            return SCIP_OKAY;
      }
      graph->targetBackward[graph->lastB] = -1;
      ++(graph->lastB);
      if( graph->lastB == graph->sizeBackward )
      {
         SCIP_CALL( checkArraySizesHeur(scip, graph, &(graph->sizeBackward), &(graph->targetBackward),
               &(graph->weightBackward), NULL, NULL, success) );

         if( !(*success) )
            return SCIP_OKAY;
      }

      /* terminate adjacency list with 0 for current level lifting */
      graph->sourceAdj[graph->levelAdj[level+1]+nAdj] = 0;
      graph->targetAdj[graph->levelAdj[level+1]+nAdj] = 0;
   }
   graph->levelAdj[level+2] = graph->levelAdj[level+1]+nAdj;

   return SCIP_OKAY;
}

/** The heuristic method for finding odd cycles by Hoffman, Padberg uses a level graph
 *  which is constructed as follows:
 *  First we choose a node (i.e. a variable of the problem or its negated) as root
 *  and assign it to level 0 (and no other node is assigned to level 0).
 *  All neighbors of the root are assigned to level 1 and the arcs between are added.
 *
 *  In general:
 *  All neighbors of nodes in level i that are assigned to level i+1, if they do not already appear in levels <= i.
 *  All arcs between nodes in level i and their neighbors are added.
 *
 *  In the construction we only take nodes that are contained in the fractional graph,
 *  i.e., their current LP values are not integral.
 *
 *  Since SCIP stores implications between original and negated variables,
 *  the level graph has at most twice the number of fractional binary variables of the problem.
 *
 *  Since the implication graph of SCIP is (normally) incomplete,
 *  it is possible to use arcs between an original variable and its negated
 *  to obtain more cycles which are valid but not found due to missing links.
 */
static
SCIP_RETCODE separateHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_SOL*             sol,                /**< given primal solution */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   /* memory for variable data */
   SCIP_VAR** varstemp;                      /* variables of the current SCIP (unsorted) */
   SCIP_VAR** vars;                          /* variables of the current SCIP (sorted if requested) */
   SCIP_Real* vals;                          /* LP-values of the variables (and negated variables) */
   unsigned int nbinvars;                    /* number of nodecandidates for implicationgraph */
   SCIP_Bool* incut;                         /* flag array for check if a variable is already covered by a cut */

   /* storage for levelgraph */
   LEVELGRAPH graph;
   unsigned int* curlevel;
   unsigned int* newlevel;
   unsigned int ncurlevel;
   unsigned int nnewlevel;
   SCIP_Bool* inlevelgraph;

   /* storage for path finding */
   unsigned int* queue;
   SCIP_Bool* inQueue;
   int* parentTree;
   int* parentTreeBackward;
   unsigned int* distance;
   SCIP_Bool* blocked;

   /* counter and limits  */
   unsigned int maxroots;                    /* maximum of level graph roots */
   unsigned int rootcounter;                 /* counter of tried roots */
   unsigned int ncutsroot;                   /* counter of cuts per root */
   unsigned int ncutslevel;                  /* counter of cuts per level */

   unsigned int i;
   unsigned int j;
   unsigned int k;
   int temp;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(result != NULL);

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &varstemp, NULL, &temp, NULL, NULL, NULL) );

   assert(varstemp != NULL || temp == 0);

   if( temp == 0 )
      return SCIP_OKAY;

   nbinvars = (unsigned int) temp;
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, (int) (2 * nbinvars)) );

   vars = NULL;
   /* duplicate variable data array for sorting (if requested) */
   if( sepadata->sortswitch != UNSORTED )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, varstemp, temp) );
   }

   switch( sepadata->sortswitch )
   {
   case UNSORTED :
      /* if no sorting is requested, we use the normal variable array */
      vars = varstemp;
      break;
   case MAXIMAL_LPVALUE :
      assert(vars != NULL);

      /* store lp-values */
      for( i = 0; i < nbinvars; ++i )
         vals[i] = SCIPgetSolVal(scip, sol, vars[i]);

      /* sort by lp-value, maximal first */
      SCIPsortDownRealPtr(vals, (void**)vars, (int) nbinvars);
      break;

   case MINIMAL_LPVALUE :
      assert(vars != NULL);

      /* store lp-values */
      for( i = 0; i < nbinvars; ++i )
         vals[i] = SCIPgetSolVal(scip, sol, vars[i]);

      /* sort by lp-value, minimal first */
      SCIPsortRealPtr(vals, (void**)vars, (int) nbinvars);
      break;

   case MAXIMAL_FRACTIONALITY  :
      assert(vars != NULL);

      /* store lp-values and determine fractionality */
      for( i = 0; i < nbinvars; ++i )
      {
         vals[i] = SCIPgetSolVal(scip, sol, vars[i]);
         vals[i] = MIN(1 - vals[i], vals[i]);
      }

      /* sort by fractionality, maximal first */
      SCIPsortDownRealPtr(vals, (void**)vars, (int) nbinvars);
      break;

   case MINIMAL_FRACTIONALITY :
      assert(vars != NULL);

      /* store lp-values and determine fractionality */
      for( i = 0; i < nbinvars; ++i )
      {
         vals[i] = SCIPgetSolVal(scip, sol, vars[i]);
         vals[i] = MIN(1 - vals[i], vals[i]);
      }

      /* sort by fractionality, minimal first */
      SCIPsortRealPtr(vals, (void**)vars, (int) nbinvars);
      break;

   default :
      SCIPerrorMessage("invalid sortswitch value\n");
      SCIPABORT();
   }
   assert(vars != NULL);

   /* create mapping for getting the index of a variable via its probindex to the index in the sorted variable array */
   SCIP_CALL( SCIPallocBufferArray(scip, &(sepadata->mapping), (int) nbinvars) );
   for( i = 0; i < nbinvars; ++i )
      sepadata->mapping[SCIPvarGetProbindex(vars[i])] = i;

   graph.n = 2*nbinvars;

   /* the implication graph is redundant and therefore more implications and clique arcs may occur than should be possible
    * @todo later: filtering of edges which were already added
    */
   /* graph.m = nbinvars*(2*nbinvars-1); */ /* = 2*nbinvars*(2*nbinvars-1)/2 */
   graph.m = INT_MAX;

   /* set sizes for graph memory storage */
   graph.sizeForward = 100*graph.n;
   graph.sizeBackward = 100*graph.n;
   graph.sizeAdj = 100*graph.n;

   /* allocate memory for level graph structure */
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.level, (int) graph.n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.beginForward, (int) graph.n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.beginBackward, (int) graph.n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.targetForward, (int) MIN(graph.sizeForward, graph.m)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.targetBackward, (int) MIN(graph.sizeBackward, graph.m)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.weightForward, (int) MIN(graph.sizeForward, graph.m)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.weightBackward, (int) MIN(graph.sizeBackward, graph.m)) );

   SCIP_CALL( SCIPallocBufferArray(scip, &curlevel, (int) graph.n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newlevel, (int) graph.n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.beginAdj, (int) graph.n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.sourceAdj, (int) MIN(graph.sizeAdj,graph.m)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.targetAdj, (int) MIN(graph.sizeAdj,graph.m)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.weightAdj, (int) MIN(graph.sizeAdj,graph.m)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.levelAdj, (int) graph.n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inlevelgraph, (int) graph.n) );

   SCIP_CALL( SCIPallocBufferArray(scip, &queue, (int) graph.n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inQueue, (int) graph.n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &parentTree, (int) graph.n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &parentTreeBackward, (int) graph.n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &distance, (int) graph.n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blocked, (int) graph.n) );

   SCIP_CALL( SCIPallocBufferArray(scip, &incut, (int) (2 * nbinvars)) );

   /* initialize LP value and cut flag for all variables */
   BMSclearMemoryArray(incut, 2*nbinvars);
   for( i = 0; i < nbinvars; ++i )
      vals[i] = SCIPgetSolVal(scip, sol, vars[i]);

   for( i = nbinvars; i < 2*nbinvars; ++i )
      vals[i] = 1 - vals[i-nbinvars];

   /* determine the number of level graph roots */
   maxroots = (unsigned int) SCIPceil(scip, sepadata->offsettestvars + (0.02 * nbinvars * sepadata->percenttestvars));
   sepadata->maxlevelsize = (unsigned int) SCIPceil(scip, sepadata->offsetnodeslevel + 0.01 * sepadata->maxpernodeslevel * graph.n);
   rootcounter = 0;

   /* check each node as root */
   for( i = (unsigned int) sepadata->lastroot; i < graph.n && rootcounter < maxroots
           && sepadata->ncuts - sepadata->oldncuts < (unsigned int) sepadata->maxsepacutsround
           && !SCIPisStopped(scip) ; ++i )
   {
      /* skip node if it is already covered by a cut and if we do not want to search cycles starting
       * with a node already covered by a cut
       */
      if( incut[i] && !(sepadata->multiplecuts) )
         continue;

      /* skip variable if its LP-value is not fractional */
      if( SCIPisFeasIntegral(scip, vals[i]) )
         continue;

      /* consider original and negated variable pair and skip variable if there is only one edge leaving the pair */
      if( (SCIPvarGetNBinImpls(vars[i % nbinvars], TRUE) + SCIPvarGetNBinImpls(vars[i % nbinvars], FALSE) < 2)
         && (SCIPvarGetNCliques(vars[i % nbinvars], TRUE) + SCIPvarGetNCliques(vars[i % nbinvars], FALSE) < 1) )
         continue;

      /* skip variable having too less implics and cliques itself */
      if( i < nbinvars )
      {
         if( SCIPvarGetNBinImpls(vars[i % nbinvars], TRUE ) < 1 && SCIPvarGetNCliques(vars[i % nbinvars], TRUE ) < 1 )
            continue;

         if( !(sepadata->addselfarcs) && SCIPvarGetNBinImpls(vars[i % nbinvars], TRUE ) < 2
            && SCIPvarGetNCliques(vars[i % nbinvars], TRUE ) < 1 )
            continue;
      }
      else
      {
         if( SCIPvarGetNBinImpls(vars[i % nbinvars], FALSE) < 1 && SCIPvarGetNCliques(vars[i % nbinvars], FALSE) < 1 )
            continue;

         if( !(sepadata->addselfarcs) && SCIPvarGetNBinImpls(vars[i % nbinvars], FALSE) < 2
            && SCIPvarGetNCliques(vars[i % nbinvars], FALSE) < 1 )
            continue;
      }

      /* node is actually considered as root node for the level graph */
      ++rootcounter;
      ncutsroot = 0;

      /* initialize graph */
      for( j = 0; j < graph.n; ++j)
      {
         graph.beginForward[j] = -1;
         graph.beginBackward[j] = -1;
         graph.beginAdj[j] = -1;
         inlevelgraph[j] = FALSE;
         blocked[j] = FALSE;
      }
      graph.lastF = 0;
      graph.lastB = 0;
      graph.nlevels = 0;
      graph.nedges = 0;

      /* insert root (first level contains root only) */
      inlevelgraph[i] = TRUE;
      graph.level[i] = 0;
      graph.levelAdj[0] = 0;
      graph.nnodes = 1;
      curlevel[0] = i;
      ncurlevel = 1;

      /* there are no arcs inside the root level */
      graph.levelAdj[graph.nlevels+1] = 0;

      /* create new levels until there are not more nodes for a new level */
      do
      {
         SCIP_Bool success;
         success = TRUE;

         /* all neighbors of nodes in level i that are assigned to level i+1,
            if they do not already appear in levels <= i. */
         SCIP_CALL( createNextLevel(scip, sepadata, vars, vals, &graph, graph.nlevels, inlevelgraph,
               curlevel, ncurlevel, newlevel, &nnewlevel, &success) );

         if( !success )
            goto TERMINATE;

         /* search for odd holes */
         if( graph.nlevels > 0 && (sepadata->includetriangles || graph.nlevels > 1) )
         {
            unsigned int maxcutslevel;

            ncutslevel = 0;

            /* calculate maximal cuts in this level due to cut limitations (per level, per root, per separation round) */
            maxcutslevel = (unsigned int) sepadata->maxcutslevel;
            maxcutslevel = (unsigned int) MIN(maxcutslevel, ncutsroot-sepadata->maxcutsroot);
            maxcutslevel = (unsigned int) MIN(maxcutslevel, sepadata->maxsepacutsround + sepadata->oldncuts - sepadata->ncuts);

            /* for each cross edge in this level find both shortest paths to root (as long as no limits are reached) */
            for( j = graph.levelAdj[graph.nlevels+1]; j < graph.levelAdj[graph.nlevels+2]
                    && ncutslevel < maxcutslevel && !SCIPisStopped(scip); ++j )
            {
               unsigned int ncyclevars;
               unsigned int u;

               /* storage for cut generation */
               unsigned int* pred; /* predecessor list */
               SCIP_Bool* incycle; /* flag array for check of double variables in found cycle */

               assert(graph.sourceAdj[j] < graph.targetAdj[j]);

               /* find shortest path from source to root and update weight of cycle */
               SCIP_CALL( findShortestPathToRoot(scip, sepadata->scale, &graph, graph.sourceAdj[j], distance, queue, inQueue, parentTree) );

#ifndef NDEBUG
               /* check that this path ends in the root node */
               u = i;
               k = 1;
               while( u != graph.sourceAdj[j] )
               {
                  assert(parentTree[u] != -1 && k <= graph.n);
                  u = (unsigned int) parentTree[u];
                  ++k;
               }
#endif

               /* block all nodes that are adjacent to nodes of the first path */
               for( k = 0; k < graph.nnodes; ++k )
                  blocked[k] = FALSE;
               SCIP_CALL( blockRootPath(scip, &graph, graph.sourceAdj[j], inlevelgraph, blocked, parentTree, i) );

               /* if the target is block, no violated odd hole can be found */
               if( blocked[graph.targetAdj[j]] )
                  continue;

               /* find shortest path from root to target node avoiding blocked nodes */
               SCIP_CALL( findUnblockedShortestPathToRoot(scip, sepadata->scale, &graph,
                     graph.targetAdj[j], distance, queue, inQueue, parentTreeBackward, i, blocked) );

               /* no odd cycle cut found */
               if( parentTreeBackward[graph.targetAdj[j]] < 0 )
                  continue;

               /* allocate and initialize predecessor list and flag array representing odd cycle */
               SCIP_CALL( SCIPallocBufferArray(scip, &pred, (int) (2 * nbinvars)) );
               SCIP_CALL( SCIPallocBufferArray(scip, &incycle, (int) (2 * nbinvars)) );
               for( k = 0; k < 2*nbinvars; ++k )
               {
                  pred[k] = DIJKSTRA_UNUSED;
                  incycle[k] = FALSE;
               }
               ncyclevars = 0;
               success = TRUE;

               /* check cycle for x-neg(x)-sub-cycles and clean them
                *  (note that a variable cannot appear twice in a cycle since it is only once in the graph)
                * convert parentTreeBackward and parentTree to pred&incycle structure for generateOddCycleCut
                */
               u = graph.targetAdj[j];

               /* add path to root to cycle */
               while( success && u != i )
               {
                  /* insert u in predecessor list */
                  pred[u] = (unsigned int) parentTreeBackward[u];

                  /* remove pairs of original and negated variable from cycle */
                  SCIP_CALL( cleanCycle(scip, pred, incycle, incut, u, graph.targetAdj[j], nbinvars, &ncyclevars,
                        sepadata->repaircycles, sepadata->allowmultiplecuts, &success) );

                  assert(parentTreeBackward[u] >= 0 || u == i);

                  /* select next node on path */
                  u = (unsigned int) parentTreeBackward[u];
               }

               /* add path from root to cycle */
               while( success && u != graph.sourceAdj[j] )
               {
                  /* insert u in predecessor list */
                  pred[u] = (unsigned int) parentTree[u];

                  /* remove pairs of original and negated variable from cycle */
                  SCIP_CALL( cleanCycle(scip, pred, incycle, incut, u, graph.targetAdj[j], nbinvars, &ncyclevars,
                        sepadata->repaircycles, sepadata->allowmultiplecuts, &success) );

                  /* select next node on path */
                  u = (unsigned int) parentTree[u];
               }
               assert(!success || u == graph.sourceAdj[j]);

               /* close the cycle */
               if( success )
               {
                  pred[u] = graph.targetAdj[j];

                  /* remove pairs of original and negated variable from cycle */
                  SCIP_CALL( cleanCycle(scip, pred, incycle, incut, u, graph.targetAdj[j], nbinvars, &ncyclevars,
                        sepadata->repaircycles, sepadata->allowmultiplecuts, &success) );
               }

               /* generate cut (if cycle is valid) */
               if(success)
               {
                  unsigned int oldncuts;
                  oldncuts = sepadata->ncuts;

                  sepadata->levelgraph = &graph;
                  SCIP_CALL( generateOddCycleCut(scip, sol, vars, nbinvars, graph.targetAdj[j], pred, ncyclevars,
                        incut, vals, sepadata, result) );
#ifndef NDEBUG
                  sepadata->levelgraph = NULL;
#endif
                  if(oldncuts < sepadata->ncuts)
                  {
                     ++ncutsroot;
                     ++ncutslevel;
                  }
               }

               /* free temporary memory */
               SCIPfreeBufferArray(scip, &incycle);
               SCIPfreeBufferArray(scip, &pred);
            }
         }

         /* copy new level to current one */
         ++(graph.nlevels);
         for( j = 0; j < nnewlevel; ++j )
            curlevel[j] = newlevel[j];
         ncurlevel = nnewlevel;
      }
      /* stop level creation loop if new level is empty or any limit is reached */
      while( nnewlevel > 0 && !SCIPisStopped(scip)
         && graph.nlevels < (unsigned int) sepadata->maxnlevels
         && ncutsroot < (unsigned int) sepadata->maxcutsroot
         && sepadata->ncuts - sepadata->oldncuts < (unsigned int) sepadata->maxsepacutsround);
   }

   /* store the last tried root (when running without sorting the variable array, we don't want
    * to always check the same variables and therefore start next time where we stopped last time)
    */
   if( sepadata->sortswitch == UNSORTED )
   {
      if( i == graph.n )
         sepadata->lastroot = 0;
      else
         sepadata->lastroot = (int) i;
   }

 TERMINATE:
   /* free memory */
   SCIPfreeBufferArray(scip, &incut);

   SCIPfreeBufferArray(scip, &blocked);
   SCIPfreeBufferArray(scip, &distance);
   SCIPfreeBufferArray(scip, &parentTreeBackward);
   SCIPfreeBufferArray(scip, &parentTree);
   SCIPfreeBufferArray(scip, &inQueue);
   SCIPfreeBufferArray(scip, &queue);

   SCIPfreeBufferArray(scip, &inlevelgraph);
   SCIPfreeBufferArray(scip, &graph.levelAdj);
   SCIPfreeBufferArray(scip, &graph.weightAdj);
   SCIPfreeBufferArray(scip, &graph.targetAdj);
   SCIPfreeBufferArray(scip, &graph.sourceAdj);
   SCIPfreeBufferArray(scip, &graph.beginAdj);
   SCIPfreeBufferArray(scip, &newlevel);
   SCIPfreeBufferArray(scip, &curlevel);

   SCIPfreeBufferArray(scip, &graph.weightBackward);
   SCIPfreeBufferArray(scip, &graph.weightForward);
   SCIPfreeBufferArray(scip, &graph.targetBackward);
   SCIPfreeBufferArray(scip, &graph.targetForward);
   SCIPfreeBufferArray(scip, &graph.beginBackward);
   SCIPfreeBufferArray(scip, &graph.beginForward);
   SCIPfreeBufferArray(scip, &graph.level);

   SCIPfreeBufferArray(scip, &(sepadata->mapping));
   if( sepadata->sortswitch != UNSORTED )
   {
      SCIPfreeBufferArray(scip, &vars);
   }
   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}

/* methods for separateGLS() */

/** memory reallocation method (the graph is normally very dense, so we dynamically allocate only the memory we need) */
static
SCIP_RETCODE checkArraySizesGLS(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int          maxarcs,            /**< maximal size of graph->head and graph->weight */
   unsigned int*         arraysize,          /**< current size of graph->head and graph->weight */
   Dijkstra_Graph*       graph,              /**< Dijkstra Graph data structure */
   SCIP_Bool*            success             /**< FALSE, iff memory reallocation fails */
   )
{
   SCIP_Real memorylimit;
   unsigned int additional;
   unsigned int j;
   unsigned int oldarraysize;

   assert(scip != NULL);
   assert(arraysize != NULL);
   assert(graph != NULL);
   assert(graph->head != NULL);
   assert(graph->weight != NULL);
   assert(success != NULL);

   SCIPdebugMessage("reallocating graph->head and graph->weight...\n");

   additional = (MIN(maxarcs, 2 * (*arraysize)) - (*arraysize)) * ((int) sizeof(*(graph->head)));
   additional += (MIN(maxarcs, 2 * (*arraysize)) - (*arraysize)) * ((int) sizeof(*(graph->weight)));

   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;

   /* if memorylimit would be exceeded or any other limit is reached free all data and exit */
   if( memorylimit <= additional/1048576.0 || SCIPisStopped(scip) )
   {
      *success = FALSE;
      SCIPdebugMessage("...memory limit exceeded\n");
      return SCIP_OKAY;
   }

   oldarraysize = *arraysize;
   *arraysize = 2*(*arraysize);

   SCIP_CALL( SCIPreallocBufferArray(scip, &(graph->head), (int) MIN(maxarcs, (*arraysize))) );
   SCIP_CALL( SCIPreallocBufferArray(scip, &(graph->weight), (int) MIN(maxarcs, (*arraysize))) );

   /* if memorylimit exceeded, leave the separator */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );

   if( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;

   if( memorylimit <= 0.0 )
   {
      SCIPdebugMessage("...memory limit exceeded - freeing all arrays\n");
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* initialize new segments of graph as empty graph */
   for( j = oldarraysize; j < MIN(maxarcs,(*arraysize)); ++j )
   {
      (graph->head)[j] = DIJKSTRA_UNUSED;
      (graph->weight)[j] = DIJKSTRA_UNUSED;
   }

   SCIPdebugMessage("...with success\n");

   return SCIP_OKAY;
}

/** add binary implications of the given node */
static
SCIP_RETCODE addGLSBinImpls(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   unsigned int          varsidx,            /**< index of current variable inside the problem variables */
   unsigned int          dijkindex,          /**< index of current variable inside the Dijkstra Graph */
   SCIP_Real*            vals,               /**< value of the variables in the given solution */
   unsigned int          nbinvars,           /**< number of binary problem variables */
   unsigned int          nbinimpls,          /**< number of binary implications of the current node */
   Dijkstra_Graph*       graph,              /**< Dijkstra Graph data structure */
   unsigned int*         narcs,              /**< current number of arcs inside the Dijkstra Graph */
   unsigned int          maxarcs,            /**< maximal number of arcs inside the Dijkstra Graph */
   SCIP_Bool             original,           /**< TRUE, iff variable is a problem variable */
   SCIP_Bool*            emptygraph,         /**< TRUE, iff there is no arc in the implication graph of the binary variables of SCIP */
   unsigned int*         arraysize,          /**< current size of graph->head and graph->weight */
   SCIP_Bool*            success             /**< FALSE, iff memory reallocation fails */
   )
{
   SCIP_VAR** implvars;         /* implications of the current variable             */
   SCIP_BOUNDTYPE* impltypes;   /* type of the implications of the current variable */
   unsigned int m;
   SCIP_VAR* neighbor;          /* current neighbor of the current variable         */
   unsigned int neighindex;
#ifndef NDEBUG
   SCIP_Real* implbounds;
#endif

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(vars != NULL);
   assert(graph != NULL);
   assert(graph->head != NULL);
   assert(graph->weight != NULL);
   assert(narcs != NULL);
   assert(emptygraph != NULL);
   assert(arraysize != NULL);
   assert(success != NULL);

   /* get implication data */
   implvars = SCIPvarGetImplVars(vars[varsidx], original);
   impltypes = SCIPvarGetImplTypes(vars[varsidx], original);
#ifndef NDEBUG
   implbounds = SCIPvarGetImplBounds(vars[varsidx], original);
#endif

   assert(implvars != NULL || nbinimpls == 0);
   assert(impltypes != NULL || nbinimpls == 0);
   assert(implbounds != NULL || nbinimpls == 0);

   /* add all implications to the graph */
   for( m = 0; m < nbinimpls; ++m )
   {
      SCIP_Real tmp;

      assert( implvars != NULL && impltypes != NULL && implbounds != NULL ); /* for lint */
      assert(SCIPvarGetType(implvars[m]) == SCIP_VARTYPE_BINARY);

      neighbor = implvars[m];
      neighindex = sepadata->mapping[SCIPvarGetProbindex(neighbor)];
      assert(neighindex < nbinvars);

      /* we use only variables with fractional LP-solution values */
      if( SCIPisFeasIntegral(scip, vals[neighindex]) )
         continue;

      /* forward direction (the backward is created at the occurrence of the current variable in the cliquevars of the neighbor) */
      /* add implication for x==1 */
      if( original )
      {
         /* implication to y=0 (I->III) */
         if( impltypes[m] == SCIP_BOUNDTYPE_UPPER )
         {
            assert(implbounds[m] == 0.0);

            tmp = SCIPfeasCeil(scip, sepadata->scale * ( 1 - vals[varsidx] - vals[neighindex] ));
            graph->weight[*narcs] = (unsigned int) MAX(0.0, tmp);
            graph->head[*narcs] = neighindex + 2 * nbinvars;
         }
         /* implication to y=1 (I->IV) */
         else
         {
            assert(impltypes[m] == SCIP_BOUNDTYPE_LOWER && implbounds[m] == 1.0 );

            tmp = SCIPfeasCeil(scip, sepadata->scale * ( 1 - vals[varsidx] - (1 - vals[neighindex]) ));
            graph->weight[*narcs] = (unsigned int) MAX(0.0, tmp);
            graph->head[*narcs] = neighindex + 3 * nbinvars;
         }
      }
      else       /* add implication for x==0 */
      {
         /* implication to y=0 (II->III) */
         if( impltypes[m] == SCIP_BOUNDTYPE_UPPER )
         {
            assert(implbounds[m] == 0.0);

            tmp = SCIPfeasCeil(scip, sepadata->scale * ( 1 - (1 - vals[varsidx]) - vals[neighindex] ));
            graph->weight[*narcs] = (unsigned int) MAX(0.0, tmp);
            graph->head[*narcs] = neighindex + 2 * nbinvars;
         }
         /* implication to y=1 (II->IV) */
         else
         {
            assert(impltypes[m] == SCIP_BOUNDTYPE_LOWER && implbounds[m] == 1.0 );

            tmp = SCIPfeasCeil(scip, sepadata->scale * (1 - (1 - vals[varsidx]) - (1 - vals[neighindex]) ));
            graph->weight[*narcs] = (unsigned int) MAX(0.0, tmp);
            graph->head[*narcs] = neighindex + 3 * nbinvars;
         }
      }

      /* update minimum and maximum weight values */
      if( graph->weight[*narcs] < graph->min_weight )
         graph->min_weight = graph->weight[*narcs];

      if( graph->weight[*narcs] > graph->max_weight )
         graph->max_weight = graph->weight[*narcs];

      assert(graph->head[*narcs] >= 2*nbinvars);
      assert(graph->head[*narcs] < 4*nbinvars);
      ++(*narcs);

      if( *arraysize == *narcs )
      {
         SCIP_CALL( checkArraySizesGLS(scip, maxarcs, arraysize, graph, success) );

         if( !(*success) )
            return SCIP_OKAY;
      }
      assert((*narcs) < maxarcs);
      ++(graph->outcnt[dijkindex]);

      *emptygraph = FALSE;
   }

   return SCIP_OKAY;
}

/** add implications from cliques of the given node */
static
SCIP_RETCODE addGLSCliques(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_VAR**            vars,               /**< problem variables */
   unsigned int          varsidx,            /**< index of current variable inside the problem variables */
   unsigned int          dijkindex,          /**< index of current variable inside the Dijkstra Graph */
   SCIP_Real*            vals,               /**< value of the variables in the given solution */
   unsigned int          nbinvars,           /**< number of binary problem variables */
   unsigned int          ncliques,           /**< number of cliques of the current node */
   Dijkstra_Graph*       graph,              /**< Dijkstra Graph data structure */
   unsigned int*         narcs,              /**< current number of arcs inside the Dijkstra Graph */
   unsigned int          maxarcs,            /**< maximal number of arcs inside the Dijkstra Graph */
   SCIP_Bool             original,           /**< TRUE, iff variable is a problem variable */
   SCIP_Bool*            emptygraph,         /**< TRUE, iff there is no arc in the implication graph of the binary variables of SCIP */
   unsigned int*         arraysize,          /**< current size of graph->head and graph->weight */
   SCIP_Bool*            success             /**< FALSE, iff memory reallocation fails */
   )
{
   SCIP_VAR* neighbor;                       /* current neighbor of the current variable */
   unsigned int neighindex;
   SCIP_CLIQUE** cliques;                    /* cliques of the current variable (with x==0/1) */
   unsigned int ncliquevars;                 /* number of variables of a clique */
   SCIP_VAR** cliquevars;                    /* variables of a clique */
   SCIP_Bool* cliquevals;                    /* is the cliquevariable fixed to TRUE or to FALSE */
   unsigned int k;
   unsigned int m;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(vars != NULL);
   assert(graph != NULL);
   assert(graph->head != NULL);
   assert(graph->weight != NULL);
   assert(narcs != NULL);
   assert(emptygraph != NULL);
   assert(arraysize != NULL);
   assert(success != NULL);

   /* if current variable has cliques of current clique-type */
   cliques = SCIPvarGetCliques(vars[varsidx], original);
   assert(cliques != NULL || ncliques == 0);

   for( k = 0; k < ncliques; ++k )
   {
      assert( cliques != NULL ); /* for lint */

      /* get clique data */
      cliquevars = SCIPcliqueGetVars(cliques[k]);
      ncliquevars = (unsigned int) SCIPcliqueGetNVars(cliques[k]);
      cliquevals = SCIPcliqueGetValues(cliques[k]);

      assert(cliquevars != NULL || ncliquevars == 0);
      assert(cliquevals != NULL || ncliquevars == 0);

      /* add arcs for all fractional variables in clique */
      for( m = 0; m < ncliquevars; ++m )
      {
         SCIP_Real tmp;

         assert( cliquevars != NULL && cliquevals != NULL ); /* for lint */
         neighbor = cliquevars[m];

         neighindex = sepadata->mapping[SCIPvarGetProbindex(neighbor)];
         assert(neighindex < nbinvars);

         /* ignore current variable */
         if( neighindex == varsidx )
            continue;

         /* we use only variables with fractional LP-solution values */
         if( SCIPisFeasIntegral(scip, vals[neighindex]) )
            continue;

         /* forward direction (the backward is created at the occurrence of the current variable in the cliquevars of the neighbor) */
         /* x==1 */
         if( original )
         {
            /* implication to y=0 (I->III) */
            if( cliquevals[m] )
            {
               tmp = SCIPfeasCeil(scip, sepadata->scale * (1 - vals[varsidx] - vals[neighindex] ));
               graph->weight[*narcs] = (unsigned int) MAX(0.0, tmp);
               graph->head[*narcs] = neighindex + 2 * nbinvars;
            }
            /* implication to y=1 (I->IV) (cliquevals[m] == FALSE) */
            else
            {
               tmp = SCIPfeasCeil(scip, sepadata->scale * (1 - vals[varsidx] - (1 - vals[neighindex]) ));
               graph->weight[*narcs] = (unsigned int) MAX(0.0, tmp);
               graph->head[*narcs] = neighindex + 3 * nbinvars;
            }
         }
         /* x==0 */
         else
         {
            /* implication to y=0 (II->III) */
            if( cliquevals[m] )
            {
               tmp = SCIPfeasCeil(scip, sepadata->scale * (1 - (1 - vals[varsidx]) - vals[neighindex] ));
               graph->weight[*narcs] = (unsigned int) MAX(0.0, tmp);
               graph->head[*narcs] = neighindex + 2 * nbinvars;
            }
            /* implication to y=1 (II->IV) (cliquevals[m] == FALSE) */
            else
            {
               tmp = SCIPfeasCeil(scip, sepadata->scale * (1 - (1 - vals[varsidx]) - (1-vals[neighindex]) )) ;
               graph->weight[*narcs] = (unsigned int) MAX(0.0, tmp);
               graph->head[*narcs] = neighindex + 3 * nbinvars;
            }
         }

         /* update minimum and maximum weight values */
         if( graph->weight[*narcs] < graph->min_weight )
            graph->min_weight = graph->weight[*narcs];

         if( graph->weight[*narcs] > graph->max_weight )
            graph->max_weight = graph->weight[*narcs];

         ++(*narcs);
         if( *arraysize == *narcs )
         {
            SCIP_CALL( checkArraySizesGLS(scip, maxarcs, arraysize, graph, success) );

            if( !(*success) )
               return SCIP_OKAY;
         }
         assert((*narcs) < maxarcs);
         ++(graph->outcnt[dijkindex]);

         *emptygraph = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** The classical method for finding odd cycles by Groetschel, Lovasz, Schrijver uses a bipartite graph
 *  which contains in each partition a node for every node in the original graph.
 *  All arcs uv of the original graph are copied to arcs from u of the first partition to v' of the second partition
 *  and from u' of the second partition to v of the first partition.
 *  A Dijkstra algorithm is used to find a path from a node x to its copy x', if existing.
 *  The nodes in the original graph corresponding to the nodes on the path form an odd cycle.
 *
 *  Since SCIP stores implications between original and negated variables,
 *  our original graph has at most twice the number of binary variables of the problem.
 *  By creating the bipartite graph we gain 4 segments of the graph:
 *
 *  I   - nodes of the original variables in the first bipartition \n
 *  II  - nodes of the negated variables in the first bipartition \n
 *  III - nodes of the original variables in the second bipartition \n
 *  IV  - nodes of the negated variables in the second bipartition
 *
 *  The length of every segment if the number of binary variables in the original problem.
 *
 *  Since the implication graph of SCIP is (normally) incomplete,
 *  it is possible to use arcs between an original variable and its negated
 *  to obtain more cycles which are valid but not found due to missing links.
 */
static
SCIP_RETCODE separateGLS(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_SOL*             sol,                /**< given primal solution */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_Bool emptygraph;                     /* flag if graph contains an arc */

   SCIP_Real* vals;                          /* values of the variables in the given solution */
   SCIP_Bool* incut;

   unsigned int i;
   unsigned int j;

   SCIP_VAR** varstemp;                      /* variables of the current SCIP (unsorted) */
   SCIP_VAR** vars;                          /* variables of the current SCIP (sorted if requested) */
   unsigned int nbinvars;                    /* number of binary problem variables */
   SCIP_Bool original;                       /* flag if the current variable is original or negated */
   int ntempvars;

   unsigned int nbinimpls;                   /* number of binary implications of the current variable */
   unsigned int ncliques;                    /* number of cliques of the current variable */

   Dijkstra_Graph graph;                     /* Dijkstra graph data structure */
   unsigned int arraysize;                   /* current size of graph->head and graph->weight */
   unsigned int narcs;                       /* number of arcs in the Dijkstra graph */
   unsigned int maxarcs;                     /* maximum number of arcs in the Dijkstra graph */
   unsigned int maxstarts;                   /* maximum number of start nodes */
   unsigned int startcounter;                /* counter of tried start nodes */

   unsigned int startnode;                   /* start node for Dijkstra algorithm */
   unsigned int endnode;                     /* target node for Dijkstra algorithm */
   unsigned long long* dist;                 /* distance matrix for Dijkstra algorithm */
   unsigned int* pred;                       /* predecessor list for found cycle */
   unsigned int* entry;                      /* storage for Dijkstra algorithm */
   unsigned int* order;                      /* storage for Dijkstra algorithm */
   unsigned int dijkindex;
   SCIP_Bool success;                        /* flag for check for several errors */

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(result != NULL);

   success = TRUE;
   emptygraph = TRUE;

   SCIP_CALL( SCIPgetVarsData(scip, &varstemp, NULL, &ntempvars, NULL, NULL, NULL) );
   assert(varstemp != NULL || ntempvars == 0);

   if( ntempvars == 0 )
      return SCIP_OKAY;

   nbinvars = (unsigned int) ntempvars;

   /* initialize flag array to avoid multiple cuts per variable, if requested by user-flag */
   SCIP_CALL( SCIPallocBufferArray(scip, &incut, (int) (2 * nbinvars)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, (int) (2 * nbinvars)) );

   /* duplicate variable data array for sorting (if requested) */
   vars = NULL;
   if( sepadata->sortswitch != UNSORTED )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &vars, varstemp, ntempvars) );
   }

   switch( sepadata->sortswitch )
   {
   case UNSORTED :
      /* if no sorting is requested, we use the normal variable array */
      vars = varstemp;
      break;

   case MAXIMAL_LPVALUE :
      assert(vars != NULL);

      /* store lp-values */
      for( i = 0; i < nbinvars; ++i )
         vals[i] = SCIPgetSolVal(scip, sol, vars[i]);

      /* sort by lp-value, maximal first */
      SCIPsortDownRealPtr(vals, (void**)vars, (int) nbinvars);
      break;

   case MINIMAL_LPVALUE :
      assert(vars != NULL);

      /* store lp-values */
      for( i = 0; i < nbinvars; ++i )
         vals[i] = SCIPgetSolVal(scip, sol, vars[i]);

      /* sort by lp-value, minimal first */
      SCIPsortRealPtr(vals, (void**)vars, (int) nbinvars);
      break;

   case MAXIMAL_FRACTIONALITY  :
      assert(vars != NULL);

      /* store lp-values and determine fractionality */
      for( i = 0; i < nbinvars; ++i )
      {
         vals[i] = SCIPgetSolVal(scip, sol, vars[i]);
         vals[i] = MIN(1-vals[i], vals[i]);
      }

      /* sort by fractionality, maximal first */
      SCIPsortDownRealPtr(vals, (void**)vars, (int) nbinvars);
      break;

   case MINIMAL_FRACTIONALITY :
      assert(vars != NULL);

      /* store lp-values and determine fractionality */
      for( i = 0; i < nbinvars; ++i )
      {
         vals[i] = SCIPgetSolVal(scip, sol, vars[i]);
         vals[i] = MIN(1-vals[i], vals[i]);
      }

      /* sort by fractionality, minimal first */
      SCIPsortRealPtr(vals, (void**)vars, (int) nbinvars);
      break;

   default :
      SCIPerrorMessage("invalid sortswitch value\n");
      SCIPABORT();
   }
   assert(vars != NULL);

   /* create mapping for getting the index of a variable via its probindex to the index in the sorted variable array */
   SCIP_CALL( SCIPallocBufferArray(scip, &(sepadata->mapping), (int) nbinvars) );
   BMSclearMemoryArray(incut, 2 * nbinvars);

   /* initialize LP value and cut flag for all variables */
   for( i = 0; i < nbinvars; ++i )
   {
      assert( 0 <= SCIPvarGetProbindex(vars[i]) && SCIPvarGetProbindex(vars[i]) < (int) nbinvars );  /* since binary variables are first */
      sepadata->mapping[SCIPvarGetProbindex(vars[i])] = i;
      vals[i] = SCIPgetSolVal(scip, sol, vars[i]); /* need to get new values, since they might be corrupted */
   }

   for( i = nbinvars; i < 2*nbinvars; ++i )
      vals[i] = 1 - vals[i - nbinvars];

   /* initialize number of nodes in  Dijkstra graph (2*2*n nodes in a mirrored bipartite graph with negated variables) */
   graph.nodes = 4*nbinvars;

   /* initialize number of arcs in  Dijkstra graph
    * (nbinvars-1 possible arcs per node (it is not possible to be linked to variable and negated)
    * + 1 self-arc (arc to negated variable)
    * + 1 dummy arc for Dijkstra data structure
    * = nbinvars+1 arcs per node
    * * graph.nodes
    * = (nbinvars+1)*graph.nodes
    * + graph.nodes => separating entries for arclist)
    */
   graph.arcs = (nbinvars+1) * graph.nodes;

   /* the implication graph is redundant and therefore more implications and clique arcs may occur than should be possible
    * @todo later: filtering of edges which were already added,  maxarcs should be graph.arcs rather than INT_MAX;
    */
   maxarcs = INT_MAX;

   /* allocate memory for Dijkstra graph arrays */
   arraysize = 100 * graph.nodes;
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.outbeg, (int) graph.nodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.outcnt, (int) graph.nodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.head, (int) MIN(maxarcs, arraysize)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &graph.weight, (int) MIN(maxarcs, arraysize)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dist, (int) graph.nodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pred, (int) graph.nodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &entry, (int) graph.nodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &order, (int) graph.nodes) );

   /* initialize Dijkstra graph as empty graph */
   for( i = 0; i < MIN(arraysize,maxarcs); ++i )
   {
      graph.head[i] = DIJKSTRA_UNUSED;
      graph.weight[i] = DIJKSTRA_UNUSED;
   }
   graph.min_weight = DIJKSTRA_FARAWAY;
   graph.max_weight = 0;
   narcs = 0;

#ifndef NDEBUG
   for( i = 0; i < graph.nodes; ++i )
   {
      graph.outbeg[i] = 0;
      graph.outcnt[i] = 0;
   }
#endif

   /* add arcs from first to second partition to Dijkstra graph (based on the original fractional implication graph) */
   for( dijkindex = 0; dijkindex < 2*nbinvars; ++dijkindex )
   {
      graph.outbeg[dijkindex] = narcs;
      graph.outcnt[dijkindex] = 0;

      /* decide if we have original or negated variable */
      if( dijkindex < nbinvars )
      {
         i = dijkindex;
         original = TRUE;
      }
      else
      {
         i = dijkindex - nbinvars;
         original = FALSE;
      }
      assert(i < nbinvars);

      /* if the variable has a fractional value we add it to the graph */
      if( !SCIPisFeasIntegral(scip, vals[i]) )
      {
         nbinimpls = (unsigned int) SCIPvarGetNBinImpls(vars[i], original);
         ncliques =  (unsigned int) SCIPvarGetNCliques(vars[i], original);

         /* insert arcs for binary implications (take var => getImpl(Bin) => getImplVar => add forward-arc) */
         /* add implications of implication-type "original" if current variable has them */
         if( nbinimpls >= 1 )
         {
            /* implications from x = 1/0 to y = 0/1 (I/II -> III/IV) */
            SCIP_CALL( addGLSBinImpls(scip, sepadata, vars, i, dijkindex, vals, nbinvars, nbinimpls, &graph,
                  &narcs, maxarcs, original, &emptygraph, &arraysize, &success) );

            if( !success )
               goto TERMINATE;
         }

         /* insert arcs for cliques (take var => getCliques => take cliquevar => add forward-arc) */
         /* add clique arcs of clique-type "original" if current variable has them */
         if( ncliques >= 1 )
         {
            /* x==1/0 -> y==0/1 (I/II -> III/IV) */
            SCIP_CALL( addGLSCliques(scip, sepadata, vars, i, dijkindex, vals, nbinvars, ncliques, &graph,
                  &narcs, maxarcs, original, &emptygraph, &arraysize, &success) );

            if( !success )
               goto TERMINATE;
         }
      }

      /* add link to copy of negated variable (useful if/because the implication graph is incomplete) */
      if( sepadata->addselfarcs && graph.outcnt[dijkindex] > 0 )
      {
         /* I -> IV */
         if( original )
         {
            assert(dijkindex < nbinvars);
            graph.head[narcs] = dijkindex + 3*nbinvars;
         }
         /* II -> III */
         else
         {
            assert(dijkindex >= nbinvars && dijkindex < 2*nbinvars);
            graph.head[narcs] = dijkindex + nbinvars;
         }
         graph.weight[narcs] = 0;

         /* update minimum and maximum weight values */
         if( graph.weight[narcs] < graph.min_weight )
            graph.min_weight = graph.weight[narcs];

         if( graph.weight[narcs] > graph.max_weight )
            graph.max_weight = graph.weight[narcs];

         ++narcs;
         if( arraysize == narcs )
         {
            SCIP_CALL( checkArraySizesGLS(scip, maxarcs, &arraysize, &graph, &success) );
            if( !success )
               goto TERMINATE;
         }
         assert(narcs < maxarcs);
         ++(graph.outcnt[dijkindex]);
      }

      /* add separating arc */
      graph.head[narcs] = DIJKSTRA_UNUSED;
      graph.weight[narcs] = DIJKSTRA_UNUSED;
      ++narcs;
      if( arraysize == narcs )
      {
         SCIP_CALL( checkArraySizesGLS(scip, maxarcs, &arraysize, &graph, &success) );
         if( !success )
            goto TERMINATE;
      }
      assert(narcs < maxarcs);
   }

   /* if the graph is empty, there is nothing to do */
   if( emptygraph )
      goto TERMINATE;

   /* add arcs from second to first partition to Dijkstra graph */
   for( i = 0; i < 2*nbinvars; ++i )
   {
      graph.outbeg[2*nbinvars+i] = narcs;
      graph.outcnt[2*nbinvars+i] = 0;

      /* copy all arcs to head from the second to the first bipartition */
      for( j = graph.outbeg[i]; j < graph.outbeg[i] + graph.outcnt[i]; ++j )
      {
         /* there are only arcs from first bipartition to the second */
         assert(graph.head[j] >= 2*nbinvars && graph.head[j] < 4*nbinvars);

         /* the backward arcs head from III->I or IV->II */
         graph.head[narcs] = graph.head[j]-2*nbinvars;
         graph.weight[narcs] = graph.weight[j];
         ++narcs;
         if( arraysize == narcs )
         {
            SCIP_CALL( checkArraySizesGLS(scip, maxarcs, &arraysize, &graph, &success) );

            if( !success )
               goto TERMINATE;
         }
         assert(narcs < maxarcs);
         ++(graph.outcnt[2*nbinvars+i]);
      }

      /* add separating arc */
      graph.head[narcs] = DIJKSTRA_UNUSED;
      graph.weight[narcs] = DIJKSTRA_UNUSED;
      ++narcs;

      if( arraysize == narcs )
      {
         SCIP_CALL( checkArraySizesGLS(scip, maxarcs, &arraysize, &graph, &success) );

         if( !success )
            goto TERMINATE;
      }
      assert(narcs < maxarcs);
   }

   SCIPdebugMessage("--- graph successfully created (%u nodes, %u arcs) ---\n", graph.nodes, narcs);

   /* graph is now prepared for Dijkstra methods */
   assert(Dijsktra_graphIsValid(&graph));

   /* determine the number of start nodes */
   maxstarts = (unsigned int) SCIPceil(scip, sepadata->offsettestvars + (0.02 * nbinvars * sepadata->percenttestvars));
   startcounter = 0;

   /* separate odd cycle inequalities by GLS method */
   for( i = (unsigned int) sepadata->lastroot; i < 2 * nbinvars
           && startcounter < maxstarts
           && sepadata->ncuts - sepadata->oldncuts < (unsigned int) sepadata->maxsepacutsround
           && !SCIPisStopped(scip); ++i )
   {
      SCIP_Bool* incycle;                    /* flag array if variable is contained in the found cycle */
      unsigned int* pred2;                   /* temporary predecessor list for backprojection of found cycle */
      unsigned int ncyclevars;               /* cycle length */
      SCIP_Bool edgedirection;               /* partitionindicator for backprojection from bipartite graph to original graph:
                                              * is the current edge a backwards edge, i.e., it goes from second to first partition? */

      /* skip isolated node */
      if( graph.head[graph.outbeg[i]] == DIJKSTRA_UNUSED )
         continue;

      /* if node has only one arc, there is no odd cycle containing this node
       * (but there are invalid odd circuits containing the only neighbor twice)
       */
      if( graph.head[graph.outbeg[i]+1] == DIJKSTRA_UNUSED )
         continue;

      /* skip node if it is already covered by a cut and
       * we do not want to search cycles starting with a node already covered by a cut
       */
      if( incut[i] && !(sepadata->multiplecuts) )
         continue;

      startcounter++;

      /* search shortest path from node to its counter part in the other partition */
      startnode = i;
      endnode = i + 2*nbinvars;

      (void) graph_dijkstra_bh(&graph, startnode, dist, pred, entry, order);

      /* no odd cycle cut found */
      if( dist[endnode] == DIJKSTRA_FARAWAY )
         continue;

      /* detect cycle including:
       * project bipartitioned graph to original graph of variables and their negated
       * (pred&incycle-structure for generateOddCycleCut)
       * check cycles for double variables and try to clean variable-negated-sub-cycles if existing
       */

      /* allocate and initialize predecessor list and flag array representing odd cycle */
      SCIP_CALL( SCIPallocBufferArray(scip, &pred2, (int) (2 * nbinvars)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &incycle, (int) (2 * nbinvars)) );

      for( j = 0; j < 2*nbinvars; ++j )
      {
         pred2[j] = DIJKSTRA_UNUSED;
         incycle[j] = FALSE;
      }

      ncyclevars = 0;
      edgedirection = TRUE;
      success = TRUE;

      /* construct odd cycle in implication graph from shortest path on bipartite graph */
      for( dijkindex = endnode; dijkindex != startnode && success; dijkindex = pred[dijkindex], edgedirection = !edgedirection )
      {
         if( edgedirection )
         {
            /* check that current node is in second partition and next node is in first partition */
            assert(dijkindex >= 2*nbinvars && dijkindex < 4*nbinvars);
            assert(pred[dijkindex] < 2*nbinvars);

            pred2[dijkindex - 2*nbinvars] = pred[dijkindex];

            /* check whether the object found is really a cycle without sub-cycles
             * (sub-cycles may occur in case there is not violated odd cycle inequality)
             * and remove pairs of original and negated variable from cycle
             */
            SCIP_CALL( cleanCycle(scip, pred2, incycle, incut, dijkindex-2*nbinvars, endnode-2*nbinvars, nbinvars,
                  &ncyclevars, sepadata->repaircycles, sepadata->allowmultiplecuts, &success) );
         }
         else
         {
            /* check that current node is in first partition and next node is in second partition */
            assert(dijkindex < 2*nbinvars);
            assert(pred[dijkindex] >= 2*nbinvars && pred[dijkindex] < 4*nbinvars);

            pred2[dijkindex] = pred[dijkindex] - 2*nbinvars;

            /* check whether the object found is really a cycle without sub-cycles
             * (sub-cycles may occur in case there is not violated odd cycle inequality)
             * and remove pairs of original and negated variable from cycle
             */
            SCIP_CALL( cleanCycle(scip, pred2, incycle, incut, dijkindex, endnode-2*nbinvars, nbinvars, &ncyclevars,
                  sepadata->repaircycles, sepadata->allowmultiplecuts, &success) );
         }
      }

      if( success )
      {
         /* generate cut */
         sepadata->dijkstragraph = &graph;
         SCIP_CALL( generateOddCycleCut(scip, sol, vars, nbinvars, startnode, pred2, ncyclevars, incut, vals, sepadata, result) );
#ifndef NDEBUG
         sepadata->dijkstragraph = NULL;
#endif
      }

      /* free temporary memory */
      SCIPfreeBufferArray(scip, &incycle);
      SCIPfreeBufferArray(scip, &pred2);
   }

   /* store the last tried root (when running without sorting the variable array, we don't want
    * to always check the same variables and therefore start next time where we stopped last time)
    */
   if( sepadata->sortswitch == UNSORTED )
   {
      if( i == 2*nbinvars )
         sepadata->lastroot = 0;
      else
         sepadata->lastroot = (int) i;
   }

 TERMINATE:
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &order);
   SCIPfreeBufferArray(scip, &entry);
   SCIPfreeBufferArray(scip, &pred);
   SCIPfreeBufferArray(scip, &dist);
   SCIPfreeBufferArray(scip, &graph.weight);
   SCIPfreeBufferArray(scip, &graph.head);
   SCIPfreeBufferArray(scip, &graph.outcnt);
   SCIPfreeBufferArray(scip, &graph.outbeg);

   SCIPfreeBufferArray(scip, &(sepadata->mapping));

   /* remove variable data array if duplicated */
   if( sepadata->sortswitch != UNSORTED )
   {
      SCIPfreeBufferArray(scip, &vars);
   }

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &incut);

   return SCIP_OKAY;
}


/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyOddcycle)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaOddcycle(scip) );

   return SCIP_OKAY;
}


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeOddcycle)
{
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeMemory(scip, &sepadata);
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
static
SCIP_DECL_SEPAINIT(sepaInitOddcycle)
{
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   sepadata->ncuts = 0;
   sepadata->oldncuts = 0;
   sepadata->nliftedcuts = 0;
   sepadata->lastroot = 0;
   sepadata->levelgraph = NULL;
   sepadata->dijkstragraph = NULL;

   return SCIP_OKAY;
}


/** deinitialization method of separator (called before transformed problem is freed) */
#define sepaExitOddcycle NULL


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#define sepaInitsolOddcycle NULL


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#define sepaExitsolOddcycle NULL


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpOddcycle)
{
   SCIP_SEPADATA* sepadata;
   int depth;
   int ncalls;

   *result = SCIP_DIDNOTRUN;

   /* get separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   depth = SCIPgetDepth(scip);
   ncalls = SCIPsepaGetNCallsAtNode(sepa);

   /* only call separator a given number of rounds at each b&b node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;


   /* only call separator if enough binary variables are present */
   if( SCIPgetNBinVars(scip) < 3 || (!(sepadata->includetriangles) && SCIPgetNBinVars(scip) < 5))
   {
      SCIPdebugMessage("skipping separator: not enough binary variables\n");
      return SCIP_OKAY;
   }

   /* only call separator if enough fractional variables are present */
   if( SCIPgetNLPBranchCands(scip) < 3 || (!(sepadata->includetriangles) && SCIPgetNLPBranchCands(scip) < 5))
   {
      SCIPdebugMessage("skipping separator: not enough fractional variables\n");
      return SCIP_OKAY;
   }

   /* only call separator if enough implications (but not all implications should come from cliques) are present */
   if( SCIPgetNImplications(scip) < 1 || SCIPgetNImplications(scip) + SCIPgetNCliques(scip) < 3 )
   {
      SCIPdebugMessage("skipping separator: not enough implications present\n");
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;
   sepadata->oldncuts = sepadata->ncuts;

   if( depth == 0 )
      sepadata->maxsepacutsround = sepadata->maxsepacutsroot;
   else
      sepadata->maxsepacutsround = sepadata->maxsepacuts;

   /* perform the actual separation routines */
   if( sepadata->useclassical )
   {
      SCIPdebugMessage("using GLS method for finding odd cycles\n");
      SCIP_CALL( separateGLS(scip, sepadata, NULL, result) );
   }
   else
   {
      SCIPdebugMessage("using level graph heuristic for finding odd cycles\n");
      SCIP_CALL( separateHeur(scip, sepadata, NULL, result) );
   }

   if( sepadata->ncuts > 0 )
   {
      SCIPdebugMessage("added %u of %d cuts (%.2f percent lifted)\n", sepadata->ncuts - sepadata->oldncuts,
         sepadata->maxsepacutsround, (sepadata->nliftedcuts*100.0)/(1.0*sepadata->ncuts));
   }
   else
   {
      SCIPdebugMessage("no cuts added (%d allowed)\n", sepadata->maxsepacutsround);
   }
   SCIPdebugMessage("total sepatime: %.2f - total number of added cuts: %u\n", SCIPsepaGetTime(sepa), sepadata->ncuts);

   return SCIP_OKAY;
}

/** arbitrary primal solution separation method of separator */
#define sepaExecsolOddcycle NULL


/*
 * separator specific interface methods
 */

/** creates the oddcycle separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaOddcycle(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* create oddcycle separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaCopyOddcycle, sepaFreeOddcycle, sepaInitOddcycle, sepaExitOddcycle, sepaInitsolOddcycle,
         sepaExitsolOddcycle, sepaExeclpOddcycle, sepaExecsolOddcycle, sepadata) );

   /* add oddcycle separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/oddcycle/useclassical",
         "should classical search method by Groetschel, Lovasz, Schrijver be used? Otherwise use levelgraph method by Hoffman, Padberg.",
         &sepadata->useclassical, FALSE, DEFAULT_USECLASSICAL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/oddcycle/liftoddcycles",
         "should odd cycle cuts be lifted?",
         &sepadata->liftoddcycles, FALSE, DEFAULT_LIFT, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/maxsepacuts",
         "maximal number of oddcycle cuts separated per separation round",
         &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/maxsepacutsroot",
         "maximal number of oddcycle cuts separated per separation round in the root node",
         &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );

   /* add advanced parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/maxrounds",
         "maximal number of oddcycle separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/maxroundsroot",
         "maximal number of oddcycle separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/scalingfactor",
         "factor for scaling of the arc-weights",
         &sepadata->scale, TRUE, DEFAULT_SCALEFACTOR, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/oddcycle/addselfarcs",
         "add links between a variable and its negated",
         &sepadata->addselfarcs, TRUE, DEFAULT_ADDSELFARCS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/oddcycle/repaircycles",
         "try to repair violated cycles with double appearance of a variable",
         &sepadata->repaircycles, TRUE, DEFAULT_REPAIRCYCLES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/oddcycle/includetriangles",
         "separate triangles found as 3-cycles or repaired larger cycles",
         &sepadata->includetriangles, TRUE, DEFAULT_INCLUDETRIANGLES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/oddcycle/multiplecuts",
         "even if a variable is already covered by a cut, still try it as start node for a cycle search",
         &sepadata->multiplecuts, TRUE, DEFAULT_MULTIPLECUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/oddcycle/allowmultiplecuts",
         "even if a variable is already covered by a cut, still allow another cut to cover it too",
         &sepadata->allowmultiplecuts, TRUE, DEFAULT_ALLOWMULTIPLECUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/oddcycle/lpliftcoef",
         "choose lifting candidate by coef*lpvalue or only by coef",
         &sepadata->lpliftcoef, TRUE, DEFAULT_LPLIFTCOEF, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/oddcycle/recalcliftcoef",
         "calculate lifting coefficient of every candidate in every step (or only if its chosen)",
         &sepadata->recalcliftcoef, TRUE, DEFAULT_RECALCLIFTCOEF, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/sortswitch",
         "use sorted variable array (unsorted(0),maxlp(1),minlp(2),maxfrac(3),minfrac(4))",
         (int*) &sepadata->sortswitch, TRUE, DEFAULT_SORTSWITCH, 0, 4, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/oddcycle/sortrootneighbors",
         "sort level of the root neighbors by fractionality (maxfrac)",
         &sepadata->sortrootneighbors, TRUE, DEFAULT_SORTROOTNEIGHBORS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/percenttestvars",
         "percentage of variables to try the chosen method on [0-100]",
         &sepadata->percenttestvars, TRUE, DEFAULT_PERCENTTESTVARS, 0, 100, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/offsettestvars",
         "offset of variables to try the chosen method on (additional to the percentage of testvars)",
         &sepadata->offsettestvars, TRUE, DEFAULT_OFFSETTESTVARS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/maxpernodeslevel",
         "percentage of nodes allowed in the same level of the level graph [0-100]",
         &sepadata->maxpernodeslevel, TRUE, DEFAULT_MAXPERNODESLEVEL, 0, 100, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/offsetnodeslevel",
         "offset of nodes allowed in the same level of the level graph (additional to the percentage of levelnodes)",
         &sepadata->offsetnodeslevel, TRUE, DEFAULT_OFFSETNODESLEVEL, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/maxnlevels",
         "maximal number of levels in level graph",
         &sepadata->maxnlevels, TRUE, DEFAULT_MAXNLEVELS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/maxcutsroot",
         "maximal number of oddcycle cuts generated per chosen variable as root of the level graph",
         &sepadata->maxcutsroot, TRUE, DEFAULT_MAXCUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/maxcutslevel",
         "maximal number of oddcycle cuts generated in every level of the level graph",
         &sepadata->maxcutslevel, TRUE, DEFAULT_MAXCUTSLEVEL, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "separating/oddcycle/maxreference",
         "minimal weight on an edge (in level graph or bipartite graph)",
         &sepadata->maxreference, TRUE, DEFAULT_MAXREFERENCE, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
