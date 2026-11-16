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

/**@file   sepa_impliedbounds.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  implied bounds separator
 * @author Kati Wolter
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/pub_implics.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sepa.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_cut.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sepa.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_var.h"
#include "scip/sepa_impliedbounds.h"


#define SEPA_NAME              "impliedbounds"
#define SEPA_DESC              "implied bounds separator"
#define SEPA_PRIORITY               -50
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define RELCUTCOEFMAXRANGE          1.0 /**< maximal allowed range of cut coefficients, relative to 1/feastol */
#define DEFAULT_USETWOSIZECLIQUES  TRUE /**< should violated inequalities for cliques with 2 variables be separated? */
#define DEFAULT_USECLIQUEAGGR      TRUE /**< should implied bounds be aggregated over cliques of binary variables? */
#define DEFAULT_MAXCLIQUESIZE      1000 /**< maximal size of a clique used for aggregating implied bounds */

/** separator-specific data for the implied bounds separator */
struct SCIP_SepaData
{
   SCIP_Bool             usetwosizecliques;  /**< should violated inequalities for cliques with 2 variables be separated? */
   SCIP_Bool             usecliqueaggr;      /**< should implied bounds be aggregated over cliques of binary variables? */
   int                   maxcliquesize;      /**< maximal size of a clique used for aggregating implied bounds */
};

/*
 * Local methods
 */

/** adds given cut with two variables, if it is violated */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_Real             val1,               /**< given coefficient of first variable */
   SCIP_VAR*             var1,               /**< given first variable */
   SCIP_Real             solval1,            /**< current LP solution value of first variable */
   SCIP_Real             val2,               /**< given coefficient of second variable */
   SCIP_VAR*             var2,               /**< given second variable */
   SCIP_Real             solval2,            /**< current LP solution value of second variable */
   SCIP_Real             rhs,                /**< given right hand side of the cut to add */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to update number of cuts added */
   )
{
   SCIP_Real activity;

   assert(ncuts != NULL);
   assert(cutoff != NULL);
   *cutoff = FALSE;

   /* calculate activity of cut */
   activity = val1 * solval1 + val2 * solval2;
   /*SCIPdebugMsg(scip, " -> %g<%s>[%g] + %g<%s>[%g] <= %g (act: %g)\n",
     val1, SCIPvarGetName(var1), solval1, val2, SCIPvarGetName(var2), solval2, rhs, activity);*/

   /* check, if cut is violated */
   if( SCIPisEfficacious(scip, activity - rhs) )
   {
      SCIP_ROW* cut;
      char cutname[SCIP_MAXSTRLEN];

      /* create cut */
      (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "implbd%" SCIP_LONGINT_FORMAT "_%d", SCIPgetNLPs(scip), *ncuts);
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
      SCIP_CALL( SCIPaddVarToRow(scip, cut, var1, val1) );
      SCIP_CALL( SCIPaddVarToRow(scip, cut, var2, val2) );
      SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
      /* set cut rank: for implied bounds we always set to 1 */
      SCIProwChgRank(cut, 1);

#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, " -> found cut (activity = %g): ", activity);
      SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif

      SCIP_CALL( SCIPaddPoolCut(scip, cut) );
      (*ncuts)++;

      /* release cut */
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   }

   return SCIP_OKAY;
}

/** returns the bound of an implication (v == fixing) -> y <= b (y >= b), or SCIP_INVALID if no such implication exists */
static
SCIP_Real getImplBound(
   SCIP_VAR*             var,                /**< binary variable v */
   SCIP_Bool             fixing,             /**< fixing value of v */
   SCIP_VAR*             implvar,            /**< implied variable y */
   SCIP_BOUNDTYPE        impltype            /**< type of implied bound */
   )
{
   SCIP_VAR** implvars;
   SCIP_BOUNDTYPE* impltypes;
   int nimpls;
   int pos;
   int k;

   nimpls = SCIPvarGetNImpls(var, fixing);
   if( nimpls == 0 )
      return SCIP_INVALID;

   implvars = SCIPvarGetImplVars(var, fixing);

   /* implications are sorted by the implied variable, a variable may have a lower and an upper bound implication */
   if( !SCIPsortedvecFindPtr((void**)implvars, SCIPvarComp, (void*)implvar, nimpls, &pos) )
      return SCIP_INVALID;

   impltypes = SCIPvarGetImplTypes(var, fixing);
   for( k = MAX(pos - 1, 0); k <= MIN(pos + 1, nimpls - 1); ++k )
   {
      if( implvars[k] == implvar && impltypes[k] == impltype )
         return SCIPvarGetImplBounds(var, fixing)[k];
   }

   return SCIP_INVALID;
}

/** returns the key of an implication (x == fixing) -> y (<= or >=) for the set of covered implications */
static
void* coveredKey(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< binary variable x */
   SCIP_Bool             fixing,             /**< value of x */
   SCIP_VAR*             implvar,            /**< implied variable y */
   SCIP_BOUNDTYPE        impltype            /**< type of implication */
   )
{
   size_t nvars = (size_t)SCIPgetNVars(scip);

   return (void*)((((size_t)SCIPvarGetProbindex(var) * 2 + (fixing ? 1 : 0)) * nvars
         + (size_t)SCIPvarGetProbindex(implvar)) * 2 + (impltype == SCIP_BOUNDTYPE_UPPER ? 1 : 0) + 1); /*lint !e571*/
}

/** aggregates the implications on a variable y over the cliques containing a given literal
 *
 *  If the literals l_1, ..., l_k form a clique, i.e., at most one of them can be one, and l_i == 1 implies y <= p_i
 *  with p_i < ub (global upper bound of y), then
 *
 *      y + sum_i (ub - p_i) l_i <= ub
 *
 *  is valid and dominates each of the single implied bound cuts y + (ub - p_i) l_i <= ub. The case y >= p_i is
 *  analogous. A literal is either a binary variable x or its negation 1 - x.
 */
static
SCIP_RETCODE separateCliqueAggregation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            solvals,            /**< array with solution values of all problem variables */
   SCIP_HASHSET*         doneset,            /**< set of already processed (clique, variable, boundtype) triples */
   SCIP_HASHSET*         covered,            /**< set of implications covered by added aggregated cuts */
   SCIP_VAR*             var,                /**< binary variable x */
   SCIP_Bool             fixing,             /**< value of x of the implication */
   SCIP_VAR*             implvar,            /**< implied variable y */
   SCIP_BOUNDTYPE        impltype,           /**< type of implication */
   SCIP_VAR**            cutvars,            /**< buffer for cut variables, at least maxcliquesize + 1 entries */
   SCIP_Real*            cutvals,            /**< buffer for cut coefficients, at least maxcliquesize + 1 entries */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to update number of cuts added */
   )
{
   SCIP_CLIQUE** cliques;
   SCIP_Real bound;
   SCIP_Real sign;
   SCIP_Real ysolval;
   int ncliques;
   int nvars;
   int c;

   ncliques = SCIPvarGetNCliques(var, fixing);
   if( ncliques == 0 )
      return SCIP_OKAY;

   bound = impltype == SCIP_BOUNDTYPE_UPPER ? SCIPvarGetUbGlobal(implvar) : SCIPvarGetLbGlobal(implvar);
   if( SCIPisInfinity(scip, REALABS(bound)) )
      return SCIP_OKAY;

   /* write cut as sign * y + sum_i d_i l_i <= sign * bound */
   sign = impltype == SCIP_BOUNDTYPE_UPPER ? 1.0 : -1.0;
   ysolval = solvals[SCIPvarGetProbindex(implvar)];
   nvars = SCIPgetNVars(scip);
   cliques = SCIPvarGetCliques(var, fixing);

   for( c = 0; c < ncliques && !(*cutoff); ++c )
   {
      SCIP_VAR** clqvars;
      SCIP_Bool* clqvals;
      SCIP_Real activity;
      SCIP_Real rhs;
      size_t key;
      int nclqvars;
      int nterms;
      int i;

      nclqvars = SCIPcliqueGetNVars(cliques[c]);
      if( nclqvars <= 2 || nclqvars > sepadata->maxcliquesize )
         continue;

      /* skip clique if it was already used for this implied variable and bound type */
      key = ((size_t)SCIPcliqueGetId(cliques[c]) * (size_t)nvars + (size_t)SCIPvarGetProbindex(implvar)) * 2
         + (impltype == SCIP_BOUNDTYPE_UPPER ? 1 : 0) + 1;
      if( SCIPhashsetExists(doneset, (void*)key) )
         continue;
      SCIP_CALL( SCIPhashsetInsert(doneset, SCIPblkmem(scip), (void*)key) );

      clqvars = SCIPcliqueGetVars(cliques[c]);
      clqvals = SCIPcliqueGetValues(cliques[c]);

      cutvars[0] = implvar;
      cutvals[0] = sign;
      activity = sign * ysolval;
      rhs = sign * bound;
      nterms = 0;

      for( i = 0; i < nclqvars; ++i )
      {
         SCIP_Real implbound;
         SCIP_Real coef;
         SCIP_Real val;

         if( SCIPvarGetProbindex(clqvars[i]) < 0 )
            continue;

         implbound = getImplBound(clqvars[i], clqvals[i], implvar, impltype);
         if( implbound == SCIP_INVALID ) /*lint !e777*/
            continue;

         coef = sign * (bound - implbound);
         if( !SCIPisFeasPositive(scip, coef) || coef * SCIPfeastol(scip) > RELCUTCOEFMAXRANGE )
            continue;

         val = solvals[SCIPvarGetProbindex(clqvars[i])];

         /* literal x contributes coef * x, literal 1 - x contributes coef - coef * x */
         if( clqvals[i] )
         {
            cutvars[nterms + 1] = clqvars[i];
            cutvals[nterms + 1] = coef;
            activity += coef * val;
         }
         else
         {
            cutvars[nterms + 1] = clqvars[i];
            cutvals[nterms + 1] = -coef;
            activity -= coef * val;
            rhs -= coef;
         }
         ++nterms;
      }

      /* with only one implication, the cut is a usual implied bound cut */
      if( nterms <= 1 || !SCIPisEfficacious(scip, activity - rhs) )
         continue;

      {
         SCIP_ROW* cut;
         char cutname[SCIP_MAXSTRLEN];

         (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "implbdclq%" SCIP_LONGINT_FORMAT "_%d", SCIPgetNLPs(scip), *ncuts);
         SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPaddVarsToRow(scip, cut, nterms + 1, cutvars, cutvals) );
         SCIProwChgRank(cut, 1);

         if( SCIPisCutEfficacious(scip, sol, cut) )
         {
#ifdef SCIP_DEBUG
            SCIPdebugMsg(scip, " -> found clique implied bound cut (activity = %g): ", activity);
            SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif
            SCIP_CALL( SCIPaddPoolCut(scip, cut) );
            (*ncuts)++;

            /* the implications of the clique members on y are dominated by this cut */
            for( i = 1; i <= nterms; ++i )
            {
               SCIP_CALL( SCIPhashsetInsert(covered, SCIPblkmem(scip),
                     coveredKey(scip, cutvars[i], cutvals[i] > 0.0, implvar, impltype)) );
            }
         }

         SCIP_CALL( SCIPreleaseRow(scip, &cut) );
      }
   }

   return SCIP_OKAY;
}

/** searches and adds implied bound cuts that are violated by the given solution value array */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            solvals,            /**< array with solution values of all problem variables */
   SCIP_VAR**            fracvars,           /**< array of fractional variables */
   SCIP_Real*            fracvals,           /**< solution values of fractional variables */
   int                   nfracs,             /**< number of fractional variables */
   SCIP_Bool*            cutoff,             /**< whether a cutoff has been detected */
   int*                  ncuts               /**< pointer to store the number of generated cuts */
   )
{
   SCIP_CLIQUE** cliques;
   SCIP_SEPADATA* sepadata;
   SCIP_HASHSET* covered;
   int ncliques;
   int i;

   assert(solvals != NULL);
   assert(fracvars != NULL || nfracs == 0);
   assert(fracvals != NULL || nfracs == 0);
   assert(cutoff != NULL);
   assert(ncuts != NULL);

   *cutoff = FALSE;
   *ncuts = 0;
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPdebugMsg(scip, "searching for implied bound cuts\n");

   /* aggregate implied bounds over cliques first; the single implied bound cuts that are dominated by an added
    * aggregated cut are skipped below
    */
   covered = NULL;
   if( sepadata->usecliqueaggr && nfracs > 0 )
   {
      SCIP_HASHSET* doneset;
      SCIP_VAR** cutvars;
      SCIP_Real* cutvals;

      /* prepare clean clique data */
      SCIP_CALL( SCIPcleanupCliques(scip, cutoff) );

      if( *cutoff )
         return SCIP_OKAY;

      SCIP_CALL( SCIPhashsetCreate(&doneset, SCIPblkmem(scip), 64) );
      SCIP_CALL( SCIPhashsetCreate(&covered, SCIPblkmem(scip), 64) );
      SCIP_CALL( SCIPallocBufferArray(scip, &cutvars, sepadata->maxcliquesize + 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &cutvals, sepadata->maxcliquesize + 1) );

      for( i = 0; i < nfracs && !(*cutoff); i++ )
      {
         int fixing;

         assert(fracvars != NULL);

         if( SCIPvarGetType(fracvars[i]) != SCIP_VARTYPE_BINARY || SCIPvarIsImpliedIntegral(fracvars[i]) )
            continue;

         for( fixing = 0; fixing <= 1 && !(*cutoff); ++fixing )
         {
            SCIP_BOUNDTYPE* impltypes;
            SCIP_VAR** implvars;
            int nimpl;
            int j;

            if( SCIPvarGetNCliques(fracvars[i], (SCIP_Bool)fixing) == 0 )
               continue;

            nimpl = SCIPvarGetNImpls(fracvars[i], (SCIP_Bool)fixing);
            implvars = SCIPvarGetImplVars(fracvars[i], (SCIP_Bool)fixing);
            impltypes = SCIPvarGetImplTypes(fracvars[i], (SCIP_Bool)fixing);

            for( j = 0; j < nimpl && !(*cutoff); j++ )
            {
               if( SCIPvarGetProbindex(implvars[j]) < 0 )
                  continue;

               SCIP_CALL( separateCliqueAggregation(scip, sepa, sepadata, sol, solvals, doneset, covered, fracvars[i],
                     (SCIP_Bool)fixing, implvars[j], impltypes[j], cutvars, cutvals, cutoff, ncuts) );
            }
         }
      }

      SCIPfreeBufferArray(scip, &cutvals);
      SCIPfreeBufferArray(scip, &cutvars);
      SCIPhashsetFree(&doneset, SCIPblkmem(scip));

      if( *cutoff )
      {
         SCIPhashsetFree(&covered, SCIPblkmem(scip));
         return SCIP_OKAY;
      }
   }

   /* search binary variables for violated implications */
   for( i = 0; i < nfracs; i++ )
   {
      SCIP_BOUNDTYPE* impltypes;
      SCIP_Real* implbounds;
      SCIP_VAR** implvars;
      int nimpl;
      int j;

      assert(fracvars != NULL);
      assert(fracvals != NULL);

      /* only process binary variables */
      if( SCIPvarGetType(fracvars[i]) != SCIP_VARTYPE_BINARY || SCIPvarIsImpliedIntegral(fracvars[i]) )
         continue;

      /* get implications of x == 1 */
      nimpl = SCIPvarGetNImpls(fracvars[i], TRUE);
      implvars = SCIPvarGetImplVars(fracvars[i], TRUE);
      impltypes = SCIPvarGetImplTypes(fracvars[i], TRUE);
      implbounds = SCIPvarGetImplBounds(fracvars[i], TRUE);

      /*SCIPdebugMsg(scip, "%d implications for <%s>[%g] == 1\n", nimpl, SCIPvarGetName(fracvars[i]), fracvals[i]);*/

      /* try to add cuts for implications of x == 1
       *    x == 1 -> y <= p:  y <= ub + x * (p - ub)  <==>  y + (ub - p) * x <=  ub
       *    x == 1 -> y >= p:  y >= lb + x * (p - lb)  <==> -y + (p - lb) * x <= -lb
       * with lb (ub) global lower (upper) bound of y
       */
      for( j = 0; j < nimpl; j++ )
      {
         SCIP_Real solval;

         assert(implvars != NULL);
         assert(impltypes != NULL);
         assert(implbounds != NULL);

         /* consider only implications with active implvar */
         if( SCIPvarGetProbindex(implvars[j]) < 0 )
            continue;

         /* skip implications covered by an aggregated cut */
         if( covered != NULL && SCIPhashsetExists(covered, coveredKey(scip, fracvars[i], TRUE, implvars[j], impltypes[j])) )
            continue;

         solval = solvals[SCIPvarGetProbindex(implvars[j])];
         if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
         {
            SCIP_Real ub;

            /* implication x == 1 -> y <= p */
            ub = SCIPvarGetUbGlobal(implvars[j]);

            /* consider only nonredundant and numerical harmless implications */
            if( SCIPisLE(scip, implbounds[j], ub) && (ub - implbounds[j]) * SCIPfeastol(scip) <= RELCUTCOEFMAXRANGE )
            {
               /* add cut if violated */
               SCIP_CALL( addCut(scip, sepa, 1.0, implvars[j], solval, (ub - implbounds[j]), fracvars[i], fracvals[i],
                     ub, cutoff, ncuts) );
               if ( *cutoff )
                  goto TERMINATE;
            }
         }
         else
         {
            SCIP_Real lb;

            /* implication x == 1 -> y >= p */
            lb = SCIPvarGetLbGlobal(implvars[j]);
            assert(impltypes[j] == SCIP_BOUNDTYPE_LOWER);

            /* consider only nonredundant and numerical harmless implications */
            if( SCIPisGE(scip, implbounds[j], lb) && (implbounds[j] - lb) * SCIPfeastol(scip) <= RELCUTCOEFMAXRANGE )
            {
               /* add cut if violated */
               SCIP_CALL( addCut(scip, sepa, -1.0, implvars[j], solval, (implbounds[j] - lb), fracvars[i], fracvals[i],
                     -lb, cutoff, ncuts) );
               if ( *cutoff )
                  goto TERMINATE;
            }
         }
      }

      /* get implications of x == 0 */
      nimpl = SCIPvarGetNImpls(fracvars[i], FALSE);
      implvars = SCIPvarGetImplVars(fracvars[i], FALSE);
      impltypes = SCIPvarGetImplTypes(fracvars[i], FALSE);
      implbounds = SCIPvarGetImplBounds(fracvars[i], FALSE);

      /*SCIPdebugMsg(scip, "%d implications for <%s>[%g] == 0\n", nimpl, SCIPvarGetName(fracvars[i]), fracvals[i]);*/

      /* try to add cuts for implications of x == 0
       *    x == 0 -> y <= p:  y <= p + x * (ub - p)  <==>  y + (p - ub) * x <=  p
       *    x == 0 -> y >= p:  y >= p + x * (lb - p)  <==> -y + (lb - p) * x <= -p
       * with lb (ub) global lower (upper) bound of y
       */
      for( j = 0; j < nimpl; j++ )
      {
         SCIP_Real solval;

         /* consider only implications with active implvar */
         if( SCIPvarGetProbindex(implvars[j]) < 0 )
            continue;

         /* skip implications covered by an aggregated cut */
         if( covered != NULL && SCIPhashsetExists(covered, coveredKey(scip, fracvars[i], FALSE, implvars[j], impltypes[j])) )
            continue;

         solval = solvals[SCIPvarGetProbindex(implvars[j])];
         if( impltypes[j] == SCIP_BOUNDTYPE_UPPER )
         {
            SCIP_Real ub;

            /* implication x == 0 -> y <= p */
            ub = SCIPvarGetUbGlobal(implvars[j]);

            /* consider only nonredundant and numerical harmless implications */
            if( SCIPisLE(scip, implbounds[j], ub) && (ub - implbounds[j]) * SCIPfeastol(scip) < RELCUTCOEFMAXRANGE )
            {
               /* add cut if violated */
               SCIP_CALL( addCut(scip, sepa, 1.0, implvars[j], solval, (implbounds[j] - ub), fracvars[i], fracvals[i],
                     implbounds[j], cutoff, ncuts) );
               if ( *cutoff )
                  goto TERMINATE;
            }
         }
         else
         {
            SCIP_Real lb;

            /* implication x == 0 -> y >= p */
            lb = SCIPvarGetLbGlobal(implvars[j]);
            assert(impltypes[j] == SCIP_BOUNDTYPE_LOWER);

            /* consider only nonredundant and numerical harmless implications */
            if( SCIPisGE(scip, implbounds[j], lb) && (implbounds[j] - lb) * SCIPfeastol(scip) < RELCUTCOEFMAXRANGE )
            {
               /* add cut if violated */
               SCIP_CALL( addCut(scip, sepa, -1.0, implvars[j], solval, (lb - implbounds[j]), fracvars[i], fracvals[i],
                     -implbounds[j], cutoff, ncuts) );
               if ( *cutoff )
                  goto TERMINATE;
            }
         }
      }
   }

TERMINATE:
   if( covered != NULL )
      SCIPhashsetFree(&covered, SCIPblkmem(scip));

   if( *cutoff )
      return SCIP_OKAY;

   /* stop separation here if cliques should not be separated */
   if( ! sepadata->usetwosizecliques )
      return SCIP_OKAY;

   /* prepare clean clique data */
   SCIP_CALL( SCIPcleanupCliques(scip, cutoff) );

   if( *cutoff )
      return SCIP_OKAY;

   cliques = SCIPgetCliques(scip);
   ncliques = SCIPgetNCliques(scip);

   /* loop over cliques of size 2 which are essentially implications and add cuts if they are violated */
   for( i = 0; i < ncliques; ++i )
   {
      SCIP_CLIQUE* clique;
      SCIP_VAR** clqvars;
      SCIP_Bool* clqvals;
      SCIP_Real rhs;

      clique = cliques[i];
      /* only consider inequality cliques of size 2 */
      if( SCIPcliqueGetNVars(clique) != 2 || SCIPcliqueIsEquation(clique) )
         continue;

      /* get variables and values of the clique */
      clqvars = SCIPcliqueGetVars(clique);
      clqvals = SCIPcliqueGetValues(clique);

      /* clique variables should never be equal after clean up */
      assert(clqvars[0] != clqvars[1]);

      /* calculate right hand side of clique inequality, which is initially 1 and decreased by 1 for every occurence of
       * a negated variable in the clique
       */
      rhs = 1.0;
      if( ! clqvals[0] )
         rhs -= 1.0;
      if( ! clqvals[1] )
         rhs -= 1.0;

      /* Basic clique inequality is
       *
       *       cx * x + (1-cx) (1-x) + cy * y + (1-cy) * (1-y) <= 1,
       *
       * where x and y are the two binary variables in the clique and cx and cy are their clique values, where a
       * clique value of 0 means that the negation of the variable should be part of the inequality.
       * Hence, exactly one of the two possible terms for x and y has a nonzero coefficient
       */
      SCIP_CALL( addCut(scip, sepa,
            clqvals[0] ? 1.0 : -1.0, clqvars[0], SCIPgetSolVal(scip, sol, clqvars[0]),
            clqvals[1] ? 1.0 : -1.0, clqvars[1], SCIPgetSolVal(scip, sol, clqvars[1]),
            rhs, cutoff, ncuts) );

      /* terminate if cutoff was found */
      if( *cutoff )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyImpliedbounds)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);

   SCIP_STRINGEQ( SCIPsepaGetName(sepa), SEPA_NAME, SCIP_INVALIDCALL );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaImpliedbounds(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeImpliedbounds)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);
   assert(sepa != NULL);

   SCIP_STRINGEQ( SCIPsepaGetName(sepa), SEPA_NAME, SCIP_INVALIDCALL );

   /* get separation data and free it */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   SCIPfreeBlockMemory(scip, &sepadata);

   /* reset data pointer to NULL */
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpImpliedbounds)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_VAR** fracvars;
   SCIP_Real* solvals;
   SCIP_Real* fracvals;
   SCIP_Bool cutoff;
   int nvars;
   int nbinvars;
   int nfracs;
   int ncuts;

   assert(sepa != NULL);
   assert(scip != NULL);

   *result = SCIP_DIDNOTRUN;

   /* gets active problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );
   if( nbinvars == 0 )
      return SCIP_OKAY;

   /* get fractional problem variables */
   /* todo try out also separating fractional implicit integer variables */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &fracvars, &fracvals, NULL, &nfracs, NULL, NULL) );
   if( nfracs == 0 )
      return SCIP_OKAY;

   /* get solution values for all variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );
   SCIP_CALL( SCIPgetVarSols(scip, nvars, vars, solvals) );

   /* call the cut separation */
   SCIP_CALL( separateCuts(scip, sepa, NULL, solvals, fracvars, fracvals, nfracs, &cutoff, &ncuts) );

   /* adjust result code */
   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( ncuts > 0 )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &solvals);

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolImpliedbounds)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_VAR** fracvars;
   SCIP_Real* solvals;
   SCIP_Real* fracvals;
   SCIP_Bool cutoff;
   int nvars;
   int nbinvars;
   int nfracs;
   int ncuts;
   int i;

   assert(sepa != NULL);
   assert(scip != NULL);

   *result = SCIP_DIDNOTRUN;

   /* gets active problem variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );
   if( nbinvars == 0 )
      return SCIP_OKAY;

   /* get solution values for all variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, vars, solvals) );

   /* get binary problem variables that are fractional in given solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &fracvars, nbinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fracvals, nbinvars) );
   nfracs = 0;
   for( i = 0; i < nbinvars; ++i )
   {
      if( !SCIPisFeasIntegral(scip, solvals[i]) )
      {
         fracvars[nfracs] = vars[i];
         fracvals[nfracs] = solvals[i];
         nfracs++;
      }
   }

   /* call the cut separation */
   ncuts = 0;
   cutoff = FALSE;

   if( nfracs > 0 )
   {
      SCIP_CALL( separateCuts(scip, sepa, sol, solvals, fracvars, fracvals, nfracs, &cutoff, &ncuts) );
   }

   /* adjust result code */
   if ( cutoff )
      *result = SCIP_CUTOFF;
   else if ( ncuts > 0 )
      *result = SCIP_SEPARATED;
   else
      *result = SCIP_DIDNOTFIND;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &fracvals);
   SCIPfreeBufferArray(scip, &fracvars);
   SCIPfreeBufferArray(scip, &solvals);

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the impliedbounds separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaImpliedbounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create impliedbounds separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   assert(sepadata != NULL);

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpImpliedbounds, sepaExecsolImpliedbounds,
         sepadata) );
   assert(sepa != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyImpliedbounds) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeImpliedbounds) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/impliedbounds/usetwosizecliques",
         "should violated inequalities for cliques with 2 variables be separated?",
         &sepadata->usetwosizecliques, TRUE, DEFAULT_USETWOSIZECLIQUES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/impliedbounds/usecliqueaggr",
         "should implied bounds of variables be aggregated over cliques of binary variables?",
         &sepadata->usecliqueaggr, TRUE, DEFAULT_USECLIQUEAGGR, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/impliedbounds/maxcliquesize",
         "maximal size of a clique used for aggregating implied bounds",
         &sepadata->maxcliquesize, TRUE, DEFAULT_MAXCLIQUESIZE, 3, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
