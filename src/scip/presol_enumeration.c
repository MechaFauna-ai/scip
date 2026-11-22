/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_enumeration.c
 * @ingroup DEFPLUGINS_PRESOL
 * @brief  enumeration presolver: enumerates all assignments of small pure binary constraints in probing mode
 * @author Max Buckley
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <limits.h>
#include <string.h>

#include "scip/presol_enumeration.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_presol.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_presol.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_timing.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"

#define PRESOL_NAME            "enumeration"
#define PRESOL_DESC            "enumerates all assignments of small pure binary constraints in probing mode"
#define PRESOL_PRIORITY         -2000000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_MAXROWSIZE             8 /**< maximal number of unfixed binary variables of an enumerated constraint */
#define DEFAULT_MAXROWS              400 /**< maximal number of constraints that are enumerated */
#define DEFAULT_MAXFAILS               6 /**< maximal number of consecutive constraints without reduction */
#define DEFAULT_MAXOVERLAP           0.5 /**< maximal fraction of variables shared with an already selected constraint */
#define DEFAULT_MAXPROPROUNDS         -1 /**< maximal number of propagation rounds per fixing (-1: no limit) */
#define DEFAULT_MAXCALLS              10 /**< maximal number of calls per presolving run */
#define MAXROWSIZE_LIMIT              16 /**< hard upper bound for maxrowsize (2^16 assignments) */

/** presolver data */
struct SCIP_PresolData
{
   int                   maxrowsize;         /**< maximal number of unfixed binary variables of an enumerated constraint */
   int                   maxrows;            /**< maximal number of constraints that are enumerated */
   int                   maxfails;           /**< maximal number of consecutive constraints without reduction */
   SCIP_Real             maxoverlap;         /**< maximal fraction of variables shared with an already selected constraint */
   int                   maxproprounds;      /**< maximal number of propagation rounds per fixing (-1: no limit) */
   int                   maxcalls;           /**< maximal number of calls per presolving run */
   int                   ncalls;             /**< number of calls in the current presolving run */
   SCIP_Longint          lastchanges;        /**< problem changes counter at the last call (-1: not called yet) */
};

/** data collected while enumerating the assignments of one constraint */
typedef struct EnumData
{
   SCIP_VAR**            rowvars;            /**< unfixed binary variables of the constraint */
   int                   nrowvars;           /**< number of row variables */
   unsigned int*         sols;               /**< assignments of the row variables of the feasible leaves (bit i = value of rowvars[i]) */
   int                   nsols;              /**< number of feasible leaves */
   SCIP_Real*            worstlb;            /**< per problem variable: minimum lower bound over leaves where it was raised */
   SCIP_Real*            worstub;            /**< per problem variable: maximum upper bound over leaves where it was lowered */
   int*                  nlbleaves;          /**< per problem variable: number of leaves in which its lower bound was raised */
   int*                  nubleaves;          /**< per problem variable: number of leaves in which its upper bound was lowered */
   int*                  touched;            /**< problem indices of variables with nonzero counters */
   int                   ntouched;           /**< number of touched variables */
   SCIP_Real*            leaflb;             /**< per problem variable: lower bound in the current leaf (if stamped) */
   SCIP_Real*            leafub;             /**< per problem variable: upper bound in the current leaf (if stamped) */
   int*                  leafstamp;          /**< per problem variable: leaf in which leaflb/leafub were set */
   int*                  leafvars;           /**< problem indices of variables with bound changes in the current leaf */
   int                   nleafvars;          /**< number of variables with bound changes in the current leaf */
   int                   nvars;              /**< number of problem variables */
   int                   maxproprounds;      /**< maximal number of propagation rounds per fixing */
} ENUMDATA;

/** records the bound changes of all probing nodes on the path to the current leaf */
static
void recordLeaf(
   SCIP*                 scip,               /**< SCIP data structure */
   ENUMDATA*             data                /**< enumeration data */
   )
{
   SCIP_NODE* node;
   unsigned int sol;
   int stamp;
   int i;

   stamp = data->nsols + 1;
   data->nleafvars = 0;

   /* collect the tightest bounds of the variables changed on the probing path */
   for( node = SCIPgetCurrentNode(scip); node != NULL && SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE;
        node = SCIPnodeGetParent(node) )
   {
      SCIP_DOMCHG* domchg;
      int nboundchgs;
      int b;

      domchg = SCIPnodeGetDomchg(node);
      if( domchg == NULL )
         continue;

      nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
      for( b = 0; b < nboundchgs; ++b )
      {
         SCIP_BOUNDCHG* boundchg;
         SCIP_VAR* var;
         int idx;

         boundchg = SCIPdomchgGetBoundchg(domchg, b);
         var = SCIPboundchgGetVar(boundchg);
         idx = SCIPvarGetProbindex(var);
         if( idx < 0 || idx >= data->nvars )
            continue;

         if( data->leafstamp[idx] != stamp )
         {
            data->leafstamp[idx] = stamp;
            data->leaflb[idx] = SCIPvarGetLbLocal(var);
            data->leafub[idx] = SCIPvarGetUbLocal(var);
            data->leafvars[data->nleafvars++] = idx;
         }
      }
   }

   /* merge into the worst-case bounds over all leaves */
   for( i = 0; i < data->nleafvars; ++i )
   {
      SCIP_VAR* var;
      int idx;

      idx = data->leafvars[i];
      var = SCIPgetVars(scip)[idx];

      if( data->nlbleaves[idx] == 0 && data->nubleaves[idx] == 0 )
         data->touched[data->ntouched++] = idx;

      if( SCIPisGT(scip, data->leaflb[idx], SCIPvarGetLbGlobal(var)) )
      {
         data->worstlb[idx] = (data->nlbleaves[idx] == 0) ? data->leaflb[idx] : MIN(data->worstlb[idx], data->leaflb[idx]);
         ++data->nlbleaves[idx];
      }
      if( SCIPisLT(scip, data->leafub[idx], SCIPvarGetUbGlobal(var)) )
      {
         data->worstub[idx] = (data->nubleaves[idx] == 0) ? data->leafub[idx] : MAX(data->worstub[idx], data->leafub[idx]);
         ++data->nubleaves[idx];
      }
   }

   /* store the assignment of the row variables */
   sol = 0;
   for( i = 0; i < data->nrowvars; ++i )
   {
      assert(SCIPisFeasEQ(scip, SCIPvarGetLbLocal(data->rowvars[i]), SCIPvarGetUbLocal(data->rowvars[i])));
      if( SCIPvarGetLbLocal(data->rowvars[i]) > 0.5 )
         sol |= (1u << i);
   }
   data->sols[data->nsols++] = sol;
}

/** enumerates the assignments of the row variables from position pos on (depth-first, in probing mode) */
static
SCIP_RETCODE enumerate(
   SCIP*                 scip,               /**< SCIP data structure */
   ENUMDATA*             data,               /**< enumeration data */
   int                   pos                 /**< position of the next row variable */
   )
{
   SCIP_VAR* var;
   int val;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   if( pos == data->nrowvars )
   {
      recordLeaf(scip, data);
      return SCIP_OKAY;
   }

   var = data->rowvars[pos];

   /* variable fixed by propagation: only one branch */
   if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5 )
   {
      SCIP_CALL( enumerate(scip, data, pos + 1) );
      return SCIP_OKAY;
   }

   for( val = 0; val <= 1; ++val )
   {
      SCIP_Bool cutoff;
      int depth;

      depth = SCIPgetProbingDepth(scip);
      SCIP_CALL( SCIPnewProbingNode(scip) );
      SCIP_CALL( SCIPfixVarProbing(scip, var, (SCIP_Real)val) );
      SCIP_CALL( SCIPpropagateProbing(scip, data->maxproprounds, &cutoff, NULL) );

      if( !cutoff )
      {
         SCIP_CALL( enumerate(scip, data, pos + 1) );
      }

      SCIP_CALL( SCIPbacktrackProbing(scip, depth) );

      if( SCIPisStopped(scip) )
         break;
   }

   return SCIP_OKAY;
}

/** collects the unfixed variables of a constraint if all of them are binary and there are between 2 and maxrowsize */
static
SCIP_RETCODE getBinaryRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   maxrowsize,         /**< maximal number of unfixed variables */
   SCIP_VAR**            buffer,             /**< buffer for the constraint variables */
   int                   buffersize,         /**< size of buffer */
   SCIP_VAR**            rowvars,            /**< array to store the unfixed binary variables (size maxrowsize) */
   int*                  nrowvars,           /**< pointer to store the number of unfixed binary variables, 0 if unsuitable */
   SCIP_Bool*            overflow            /**< pointer to store whether the buffer was too small */
   )
{
   SCIP_Bool success;
   int nconsvars;
   int i;

   *nrowvars = 0;
   *overflow = FALSE;

   SCIP_CALL( SCIPgetConsNVars(scip, cons, &nconsvars, &success) );
   if( !success || nconsvars < 2 )
      return SCIP_OKAY;
   if( nconsvars > buffersize )
   {
      *overflow = TRUE;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPgetConsVars(scip, cons, buffer, buffersize, &success) );
   if( !success )
      return SCIP_OKAY;

   for( i = 0; i < nconsvars; ++i )
   {
      SCIP_VAR* var = buffer[i];
      int j;

      /* the enumeration covers both values, so a negated variable can be replaced by its active counterpart */
      if( SCIPvarIsNegated(var) )
         var = SCIPvarGetNegationVar(var);

      if( SCIPvarGetLbLocal(var) > 0.5 || SCIPvarGetUbLocal(var) < 0.5 )
         continue;
      if( !SCIPvarIsBinary(var) || !SCIPvarIsActive(var) || *nrowvars >= maxrowsize )
      {
         *nrowvars = 0;
         return SCIP_OKAY;
      }

      /* skip duplicates (e.g., a variable and its negation are both reported as the active variable) */
      for( j = 0; j < *nrowvars; ++j )
         if( rowvars[j] == var )
            break;
      if( j == *nrowvars )
         rowvars[(*nrowvars)++] = var;
   }

   if( *nrowvars < 2 )
      *nrowvars = 0;

   return SCIP_OKAY;
}

/** applies the reductions found by enumerating one constraint */
static
SCIP_RETCODE applyReductions(
   SCIP*                 scip,               /**< SCIP data structure */
   ENUMDATA*             data,               /**< enumeration data */
   int*                  nfixedvars,         /**< pointer to count fixed variables */
   int*                  naggrvars,          /**< pointer to count aggregated variables */
   int*                  nchgbds,            /**< pointer to count changed bounds */
   SCIP_Bool*            cutoff,             /**< pointer to store whether infeasibility was detected */
   SCIP_Bool*            found               /**< pointer to store whether a reduction was applied */
   )
{
   SCIP_VAR** vars;
   unsigned int allone;
   unsigned int allzero;
   int i;
   int j;

   *found = FALSE;
   vars = SCIPgetVars(scip);

   /* worst-case bounds: a bound is valid if it is tightened in every feasible leaf */
   for( i = 0; i < data->ntouched && !(*cutoff); ++i )
   {
      SCIP_VAR* var;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;
      int idx;

      idx = data->touched[i];
      var = vars[idx];

      if( data->nlbleaves[idx] == data->nsols )
      {
         SCIP_CALL( SCIPtightenVarLb(scip, var, data->worstlb[idx], FALSE, &infeasible, &tightened) );
         *cutoff = *cutoff || infeasible;
         if( tightened )
         {
            *found = TRUE;
            if( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) )
               ++(*nfixedvars);
            else
               ++(*nchgbds);
         }
      }
      if( !(*cutoff) && data->nubleaves[idx] == data->nsols )
      {
         SCIP_CALL( SCIPtightenVarUb(scip, var, data->worstub[idx], FALSE, &infeasible, &tightened) );
         *cutoff = *cutoff || infeasible;
         if( tightened )
         {
            *found = TRUE;
            if( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) )
               ++(*nfixedvars);
            else
               ++(*nchgbds);
         }
      }
   }

   if( *cutoff )
      return SCIP_OKAY;

   /* clique: at most one row variable is one in every feasible assignment */
   {
      SCIP_Bool atmostone = TRUE;
      SCIP_Bool exactlyone = TRUE;

      for( i = 0; i < data->nsols; ++i )
      {
         unsigned int s = data->sols[i];
         if( (s & (s - 1)) != 0 )
            atmostone = exactlyone = FALSE;
         if( s == 0 )
            exactlyone = FALSE;
      }

      if( atmostone && data->nrowvars >= 2 )
      {
         SCIP_Bool infeasible;
         int nbdchgs = 0;

         SCIP_CALL( SCIPaddClique(scip, data->rowvars, NULL, data->nrowvars, exactlyone, &infeasible, &nbdchgs) );
         *cutoff = *cutoff || infeasible;
         *nchgbds += nbdchgs;
      }
   }

   if( *cutoff )
      return SCIP_OKAY;

   /* pairs of row variables that are equal or complementary in every feasible assignment */
   allone = 0;
   allzero = 0;
   for( i = 0; i < data->nsols; ++i )
   {
      allone |= data->sols[i];
      allzero |= ~data->sols[i];
   }

   for( i = 0; i < data->nrowvars && !(*cutoff); ++i )
   {
      /* skip variables that are fixed in all assignments (handled by the worst-case bounds) */
      if( !(allone & (1u << i)) || !(allzero & (1u << i)) )
         continue;

      for( j = i + 1; j < data->nrowvars && !(*cutoff); ++j )
      {
         SCIP_Bool equal = TRUE;
         SCIP_Bool compl = TRUE;
         int s;

         if( !(allone & (1u << j)) || !(allzero & (1u << j)) )
            continue;

         for( s = 0; s < data->nsols && (equal || compl); ++s )
         {
            SCIP_Bool vi = (data->sols[s] >> i) & 1u;
            SCIP_Bool vj = (data->sols[s] >> j) & 1u;
            if( vi != vj )
               equal = FALSE;
            else
               compl = FALSE;
         }

         if( (equal || compl) && SCIPvarIsActive(data->rowvars[i]) && SCIPvarIsActive(data->rowvars[j])
            && !SCIPdoNotAggr(scip) )
         {
            SCIP_Bool infeasible;
            SCIP_Bool redundant;
            SCIP_Bool aggregated;

            /* equal: x_i - x_j = 0, complementary: x_i + x_j = 1 */
            SCIP_CALL( SCIPaggregateVars(scip, data->rowvars[i], data->rowvars[j], 1.0, equal ? -1.0 : 1.0,
                  equal ? 0.0 : 1.0, &infeasible, &redundant, &aggregated) );
            *cutoff = *cutoff || infeasible;
            if( aggregated )
            {
               ++(*naggrvars);
               *found = TRUE;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */

/** copy method for presolver plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyEnumeration)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   SCIP_CALL( SCIPincludePresolEnumeration(scip) );

   return SCIP_OKAY;
}

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeEnumeration)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** presolving initialization method of presolver (called when presolving is about to begin) */
static
SCIP_DECL_PRESOLINITPRE(presolInitpreEnumeration)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   presoldata->ncalls = 0;
   presoldata->lastchanges = -1;

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecEnumeration)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
   SCIP_CONS** conss;
   SCIP_CONS** candconss;
   SCIP_Real* candscores;
   SCIP_VAR** buffer;
   SCIP_VAR*** selrows;
   int* selnrows;
   ENUMDATA data;
   SCIP_Bool cutoff;
   int buffersize;
   int nconss;
   int ncands;
   int nsel;
   int nfails;
   int nenum;
   int nleaves;
   int oldnfixed;
   int oldnaggr;
   int oldnchgbds;
   int c;
   int i;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   if( SCIPinProbing(scip) || SCIPgetNBinVars(scip) == 0 || SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* like HiGHS, run late in presolving and again only after the problem changed (small pure binary rows appear as other
    * presolvers fix variables), with a limit on the number of calls per presolving run
    */
   if( presoldata->ncalls >= presoldata->maxcalls )
      return SCIP_OKAY;
   {
      SCIP_Longint changes = (SCIP_Longint)nfixedvars[0] + naggrvars[0] + nchgbds[0] + ndelconss[0] + nupgdconss[0]
         + nchgcoefs[0] + nchgsides[0];
      if( changes == presoldata->lastchanges )
         return SCIP_OKAY;
      presoldata->lastchanges = changes;
   }
   ++presoldata->ncalls;

   nconss = SCIPgetNConss(scip);
   if( nconss == 0 )
      return SCIP_OKAY;
   conss = SCIPgetConss(scip);

   *result = SCIP_DIDNOTFIND;

   buffersize = 4 * presoldata->maxrowsize + 16;
   SCIP_CALL( SCIPallocBufferArray(scip, &buffer, buffersize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candconss, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candscores, nconss) );

   /* collect candidate constraints and score them by the locks and cliques of their variables */
   ncands = 0;
   for( c = 0; c < nconss; ++c )
   {
      SCIP_VAR* rowvars[MAXROWSIZE_LIMIT];
      SCIP_Bool overflow;
      SCIP_Real score;
      int nrowvars;

      if( !SCIPconsIsActive(conss[c]) || SCIPconsIsModifiable(conss[c]) )
         continue;

      SCIP_CALL( getBinaryRow(scip, conss[c], presoldata->maxrowsize, buffer, buffersize, rowvars, &nrowvars, &overflow) );
      if( nrowvars == 0 )
         continue;

      score = 0.0;
      for( i = 0; i < nrowvars; ++i )
         score += SCIPvarGetNLocksDown(rowvars[i]) + SCIPvarGetNLocksUp(rowvars[i])
            + SCIPvarGetNCliques(rowvars[i], TRUE) + SCIPvarGetNCliques(rowvars[i], FALSE);

      candconss[ncands] = conss[c];
      candscores[ncands] = score / nrowvars;
      ++ncands;
   }

   if( ncands == 0 )
      goto FREECANDS;

   SCIPsortDownRealPtr(candscores, (void**)candconss, ncands);

   /* select up to maxrows constraints, skipping those that overlap too much with an already selected one */
   SCIP_CALL( SCIPallocBufferArray(scip, &selrows, presoldata->maxrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &selnrows, presoldata->maxrows) );
   nsel = 0;
   for( c = 0; c < ncands && nsel < presoldata->maxrows; ++c )
   {
      SCIP_VAR* rowvars[MAXROWSIZE_LIMIT];
      SCIP_Bool overflow;
      SCIP_Bool similar;
      int nrowvars;
      int s;

      SCIP_CALL( getBinaryRow(scip, candconss[c], presoldata->maxrowsize, buffer, buffersize, rowvars, &nrowvars, &overflow) );
      if( nrowvars == 0 )
         continue;

      similar = FALSE;
      for( s = 0; s < nsel && !similar; ++s )
      {
         int noverlap = 0;
         int a;
         int b;

         for( a = 0; a < nrowvars; ++a )
            for( b = 0; b < selnrows[s]; ++b )
               if( rowvars[a] == selrows[s][b] )
                  ++noverlap;

         if( noverlap > presoldata->maxoverlap * MIN(nrowvars, selnrows[s]) )
            similar = TRUE;
      }
      if( similar )
         continue;

      SCIP_CALL( SCIPduplicateBufferArray(scip, &selrows[nsel], rowvars, nrowvars) );
      selnrows[nsel] = nrowvars;
      ++nsel;
   }

   /* enumeration data, indexed by problem index */
   data.nvars = SCIPgetNVars(scip);
   data.maxproprounds = presoldata->maxproprounds;
   SCIP_CALL( SCIPallocBufferArray(scip, &data.sols, 1 << presoldata->maxrowsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &data.worstlb, data.nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &data.worstub, data.nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &data.nlbleaves, data.nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &data.nubleaves, data.nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &data.touched, data.nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &data.leaflb, data.nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &data.leafub, data.nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &data.leafstamp, data.nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &data.leafvars, data.nvars) );

   cutoff = FALSE;
   nfails = 0;
   nenum = 0;
   nleaves = 0;
   oldnfixed = *nfixedvars; oldnaggr = *naggrvars; oldnchgbds = *nchgbds;
   for( c = 0; c < nsel && !cutoff && nfails < presoldata->maxfails && !SCIPisStopped(scip); ++c )
   {
      SCIP_VAR* rowvars[MAXROWSIZE_LIMIT];
      SCIP_Bool found;
      int nrowvars;

      /* drop variables that were fixed or aggregated by earlier reductions */
      nrowvars = 0;
      for( i = 0; i < selnrows[c]; ++i )
      {
         SCIP_VAR* var = selrows[c][i];
         if( SCIPvarIsActive(var) && SCIPvarGetLbGlobal(var) < 0.5 && SCIPvarGetUbGlobal(var) > 0.5 )
            rowvars[nrowvars++] = var;
      }
      if( nrowvars < 2 )
         continue;

      data.rowvars = rowvars;
      data.nrowvars = nrowvars;
      data.nsols = 0;
      data.ntouched = 0;

      SCIP_CALL( SCIPstartProbing(scip) );
      SCIP_CALL( enumerate(scip, &data, 0) );
      SCIP_CALL( SCIPendProbing(scip) );
      ++nenum;
      nleaves += data.nsols;

      if( SCIPisStopped(scip) )
         found = FALSE;
      else if( data.nsols == 0 )
      {
         SCIPdebugMsg(scip, "enumeration: no feasible assignment of %d binary variables\n", nrowvars);
         cutoff = TRUE;
         found = TRUE;
      }
      else
      {
         SCIP_CALL( applyReductions(scip, &data, nfixedvars, naggrvars, nchgbds, &cutoff, &found) );
      }

      /* reset the per-variable counters of this row */
      for( i = 0; i < data.ntouched; ++i )
      {
         data.nlbleaves[data.touched[i]] = 0;
         data.nubleaves[data.touched[i]] = 0;
      }

      if( found )
      {
         nfails = 0;
         *result = SCIP_SUCCESS;
      }
      else
         ++nfails;
   }

   if( cutoff )
      *result = SCIP_CUTOFF;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
      "   (%.1fs) enumeration: %d candidates, %d selected, %d enumerated (%d feasible assignments): %d fixings, %d aggregations, %d bound changes%s\n",
      SCIPgetSolvingTime(scip), ncands, nsel, nenum, nleaves, *nfixedvars - oldnfixed, *naggrvars - oldnaggr,
      *nchgbds - oldnchgbds, cutoff ? ", infeasible" : "");

   SCIPfreeBufferArray(scip, &data.leafvars);
   SCIPfreeBufferArray(scip, &data.leafstamp);
   SCIPfreeBufferArray(scip, &data.leafub);
   SCIPfreeBufferArray(scip, &data.leaflb);
   SCIPfreeBufferArray(scip, &data.touched);
   SCIPfreeBufferArray(scip, &data.nubleaves);
   SCIPfreeBufferArray(scip, &data.nlbleaves);
   SCIPfreeBufferArray(scip, &data.worstub);
   SCIPfreeBufferArray(scip, &data.worstlb);
   SCIPfreeBufferArray(scip, &data.sols);

   for( i = nsel - 1; i >= 0; --i )
      SCIPfreeBufferArray(scip, &selrows[i]);
   SCIPfreeBufferArray(scip, &selnrows);
   SCIPfreeBufferArray(scip, &selrows);

FREECANDS:
   SCIPfreeBufferArray(scip, &candscores);
   SCIPfreeBufferArray(scip, &candconss);
   SCIPfreeBufferArray(scip, &buffer);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** creates the enumeration presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolEnumeration(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecEnumeration, presoldata) );
   assert(presol != NULL);

   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyEnumeration) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeEnumeration) );
   SCIP_CALL( SCIPsetPresolInitpre(scip, presol, presolInitpreEnumeration) );
   presoldata->ncalls = 0;
   presoldata->lastchanges = -1;

   SCIP_CALL( SCIPaddIntParam(scip, "presolving/" PRESOL_NAME "/maxrowsize",
         "maximal number of unfixed binary variables of an enumerated constraint",
         &presoldata->maxrowsize, TRUE, DEFAULT_MAXROWSIZE, 2, MAXROWSIZE_LIMIT, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "presolving/" PRESOL_NAME "/maxrows",
         "maximal number of constraints that are enumerated",
         &presoldata->maxrows, TRUE, DEFAULT_MAXROWS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "presolving/" PRESOL_NAME "/maxfails",
         "maximal number of consecutive constraints without reduction",
         &presoldata->maxfails, TRUE, DEFAULT_MAXFAILS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "presolving/" PRESOL_NAME "/maxoverlap",
         "maximal fraction of variables shared with an already selected constraint",
         &presoldata->maxoverlap, TRUE, DEFAULT_MAXOVERLAP, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "presolving/" PRESOL_NAME "/maxproprounds",
         "maximal number of propagation rounds per fixing (-1: no limit)",
         &presoldata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "presolving/" PRESOL_NAME "/maxcalls",
         "maximal number of calls per presolving run",
         &presoldata->maxcalls, TRUE, DEFAULT_MAXCALLS, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
