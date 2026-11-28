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

/**@file   sepa_path.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  directed continuous-path aggregation and path-mixing cut separator
 *
 * Separates cMIR, flow cover, and multi-row path-mixing cuts by walking directed
 * paths through tight LP rows that cancel continuous variables away from their
 * simple or variable bounds (following the HighsPathSeparator design):
 *
 * -# Tight LP rows (including pure-continuous flow-conservation equations that
 *    are skipped by fractionality-filtered aggregation heuristics) are classified
 *    into <=, >=, and == rows.
 * -# Equality rows containing a single continuous variable away from its bounds
 *    are registered as automatic substitutions and always eliminated first.
 * -# Remaining continuous variables away from their bounds are indexed by directed
 *    in-arc and out-arc row lists so every aggregation step preserves a valid <=
 *    inequality without introducing negative slacks.
 * -# Along each directed path (up to maxpathlen rows), single-row cMIR and flow
 *    cover cuts are tried; if none succeed on a path of length >= 2, a VBD-aware
 *    multi-row Günlük-Pochet path-mixing cut is attempted across the cumulative
 *    path inequalities.
 * -# Candidate cuts are filtered with SCIPselectCutsHybrid; global cuts go to the
 *    cut pool and local cuts to the separation storage.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <math.h>
#include <string.h>
#include "blockmemshell/memory.h"
#include "scip/cuts.h"
#include "scip/cutsel_hybrid.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sepa.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_cut.h"
#include "scip/scip_exact.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_randnumgen.h"
#include "scip/scip_sepa.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_var.h"
#include "scip/sepa_path.h"

#define SEPA_NAME              "path"
#define SEPA_DESC              "directed continuous-path aggregation and path-mixing cut separator"
#define SEPA_PRIORITY             -3500
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE
#define SEPA_DELAY                FALSE

#define DEFAULT_MAXROUNDS             1 /**< maximal number of path separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        15 /**< maximal number of path separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of path cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     200 /**< maximal number of path cuts separated per separation round in root node */
#define DEFAULT_MAXCUTCANDS        1000 /**< maximal number of cut candidates collected per separation round */
#define DEFAULT_MAXTRIES             50 /**< maximal number of start rows tried per separation round in tree */
#define DEFAULT_MAXTRIESROOT        500 /**< maximal number of start rows tried per separation round in root */
#define DEFAULT_MAXPATHLEN            6 /**< maximal path length for directed row aggregation */
#define DEFAULT_MAXAGGDENSITY      0.20 /**< maximal density of aggregated row */
#define DEFAULT_DENSITYOFFSET       100 /**< additional nonzeros allowed on top of density */
#define DEFAULT_USEPATHMIX         TRUE /**< should multi-row path-mixing cuts be separated? */
#define DEFAULT_TRYFLOWCOVER       TRUE /**< should flow cover cuts be tried along paths? */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if no longer tight? */
#define DEFAULT_INITSEED         0x51A7 /**< initial seed for random arc selection and cut selection */

#define BOUNDSWITCH                 0.5 /**< threshold for bound switching in MIR transformation */
#define POSTPROCESS                TRUE /**< apply postprocessing after MIR calculation */
#define VARTYPEUSEVBDS                2 /**< allow VBD substitution on continuous variables */
#define MINFRAC                    0.05 /**< minimal fractionality of rhs to produce MIR cut for */
#define MAXFRAC                   0.999 /**< maximal fractionality of rhs to produce MIR cut for */
#define MAXTESTDELTA                  5 /**< maximal number of deltas tested by cMIR heuristic */
#define MAXSUBPASSES                 10 /**< maximal substitution passes per path step */
#define MAXPATHLENHARD               16 /**< hard cap on stack path array */
#define MAXAGGRNNZHARD              300 /**< hard cap on nonzeros in aggregated path rows */
#define MAXARCDEG                    30 /**< skip continuous variables appearing in too many rows */
#define MAXWORKPERROUND          500000 /**< per-round effort cap in aggregated nonzero evaluations */

/* row classification */
#define ROWTYPE_UNUSABLE             -2
#define ROWTYPE_GEQ                  -1
#define ROWTYPE_EQ                    0
#define ROWTYPE_LEQ                   1

/* consistent bound transformation types for path-mixing */
#define PMBND_NONE                    0
#define PMBND_SIMPLE_LB               1
#define PMBND_SIMPLE_UB               2
#define PMBND_VLB                     3
#define PMBND_VUB                     4

/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_Real             maxaggdensity;      /**< maximal density of aggregated row */
   int                   maxrounds;          /**< maximal number of separation rounds per node */
   int                   maxroundsroot;      /**< maximal number of separation rounds at root */
   int                   maxsepacuts;        /**< maximal number of cuts selected per round */
   int                   maxsepacutsroot;    /**< maximal number of cuts selected per round at root */
   int                   maxcutcands;        /**< maximal number of cut candidates collected per round */
   int                   maxtries;           /**< maximal start rows per round in tree */
   int                   maxtriesroot;       /**< maximal start rows per round at root */
   int                   maxpathlen;         /**< maximal path length */
   int                   densityoffset;      /**< additional nonzeros allowed on top of density */
   int                   initseed;           /**< initial random seed */
   SCIP_Bool             usepathmix;         /**< separate multi-row path-mixing cuts? */
   SCIP_Bool             tryflowcover;       /**< try flow cover cuts along paths? */
   SCIP_Bool             dynamiccuts;        /**< are generated cuts dynamic? */
};

/** variable bound metadata for continuous variables */
typedef struct PathVarInfo
{
   SCIP_Real             solval;             /**< current LP solution value */
   SCIP_Real             lb;                 /**< active lower bound (local or global) */
   SCIP_Real             ub;                 /**< active upper bound (local or global) */
   SCIP_Real             bounddist;          /**< distance to closest simple or variable bound */
   SCIP_Real             vlbcoef;            /**< best binary VLB coefficient (x >= coef*y + constant) */
   SCIP_Real             vlbconst;           /**< best binary VLB constant */
   SCIP_Real             vubcoef;            /**< best binary VUB coefficient (x <= coef*y + constant) */
   SCIP_Real             vubconst;           /**< best binary VUB constant */
   SCIP_Real             vlbdist;            /**< distance to best binary VLB */
   SCIP_Real             vubdist;            /**< distance to best binary VUB */
   int                   vlbvaridx;          /**< probindex of best binary VLB variable, or -1 */
   int                   vubvaridx;          /**< probindex of best binary VUB variable, or -1 */
   int                   substrow;           /**< single-cont equality row for substitution, or -1 */
   SCIP_Real             substcoef;          /**< coefficient in substrow */
   SCIP_Bool             iscont;             /**< TRUE iff !SCIPvarIsIntegral(var) */
   signed char           pmbndtype;          /**< preselected consistent bound type for path-mixing */
} PATHVARINFO;

/** snapshot of one cumulative valid <= inequality along a directed path */
typedef struct PathRowSnap
{
   int*                  inds;               /**< problem variable indices */
   SCIP_Real*            vals;               /**< coefficients */
   SCIP_Real             rhs;                /**< right-hand side */
   int                   nnz;                /**< number of nonzeros */
   int                   rank;               /**< cut rank */
   SCIP_Bool             islocal;            /**< is inequality local? */
} PATHROWSNAP;

/** directed row-arc entry */
typedef struct PathArc
{
   int                   row;                /**< LP row index */
   SCIP_Real             coef;               /**< variable coefficient in row */
} PATHARC;

/*
 * Local helper methods
 */

/** creates a SCIP_ROW from coefficients and adds it to the candidate array if efficacious */
static
SCIP_RETCODE storeCutCandidate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_VAR**            vars,               /**< active problem variables */
   const SCIP_Real*      cutcoefs,           /**< cut coefficients */
   const int*            cutinds,            /**< variable probindices */
   int                   cutnnz,             /**< number of nonzeros */
   SCIP_Real             cutrhs,             /**< cut rhs */
   int                   cutrank,            /**< cut rank */
   SCIP_Bool             cutislocal,         /**< is cut local? */
   const char*           prefix,             /**< cut name prefix */
   SCIP_ROW**            cuts,               /**< candidate cuts array */
   int*                  ncuts,              /**< number of stored candidates */
   int                   maxcuts,            /**< capacity of cuts array */
   SCIP_Bool*            added               /**< pointer to store whether candidate was kept */
   )
{
   SCIP_ROW* cut;
   char cutname[SCIP_MAXSTRLEN];
   SCIP_Real maxabs;
   SCIP_Real minabs;
   int i;

   assert(added != NULL);
   *added = FALSE;

   if( cutnnz <= 0 || *ncuts >= maxcuts || SCIPisInfinity(scip, REALABS(cutrhs)) )
      return SCIP_OKAY;

   maxabs = 0.0;
   minabs = SCIP_REAL_MAX;
   for( i = 0; i < cutnnz; ++i )
   {
      SCIP_Real a = REALABS(cutcoefs[i]);
      if( a > maxabs )
         maxabs = a;
      if( a > 0.0 && a < minabs )
         minabs = a;
   }
   if( minabs == SCIP_REAL_MAX || maxabs <= 0.0 || maxabs > 1e6 * minabs )
      return SCIP_OKAY;

   (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "%s%" SCIP_LONGINT_FORMAT "_%d", prefix, SCIPgetNLPs(scip), *ncuts);
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), cutrhs,
         cutislocal, FALSE, sepadata->dynamiccuts) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
   for( i = 0; i < cutnnz; ++i )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[cutinds[i]], cutcoefs[i]) );
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
   SCIProwChgRank(cut, cutrank);

   if( SCIProwGetNNonz(cut) > 0 && SCIPisCutEfficacious(scip, NULL, cut) && SCIPisCutNew(scip, cut) )
   {
      cuts[(*ncuts)++] = cut;
      *added = TRUE;
   }
   else
   {
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   }

   return SCIP_OKAY;
}

/** tries to generate a single-row cMIR or flow cover cut from aggrrow */
static
SCIP_RETCODE trySingleRowCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_VAR**            vars,               /**< active problem variables */
   SCIP_AGGRROW*         aggrrow,            /**< current aggregation row */
   SCIP_Bool             allowlocal,         /**< allow local cuts? */
   SCIP_Real*            cutcoefs,           /**< scratch coefs buffer (nvars) */
   int*                  cutinds,            /**< scratch inds buffer (nvars) */
   SCIP_Real*            bestcoefs,          /**< best coefs buffer (nvars) */
   int*                  bestinds,           /**< best inds buffer (nvars) */
   SCIP_ROW**            cuts,               /**< candidate cuts array */
   int*                  ncuts,              /**< candidate count */
   int                   maxcuts,            /**< candidate capacity */
   SCIP_Bool*            found               /**< set to TRUE if a cut candidate was added */
   )
{
   SCIP_Real cutefficacy;
   SCIP_Real besteff = 0.0;
   SCIP_Real cutrhs = 0.0;
   SCIP_Real bestrhs = 0.0;
   int cutnnz = 0;
   int bestnnz = 0;
   int cutrank = 0;
   int bestrank = 0;
   SCIP_Bool cutislocal = FALSE;
   SCIP_Bool bestlocal = FALSE;
   SCIP_Bool success = FALSE;
   const char* bestprefix = "pathcmir";

   assert(found != NULL);
   *found = FALSE;

   if( SCIPaggrRowGetNNz(aggrrow) == 0 || *ncuts >= maxcuts )
      return SCIP_OKAY;

   /* try cMIR scaling heuristic (includes VBD substitution on continuous vars and negative scaling) */
   cutefficacy = 0.0;
   SCIP_CALL( SCIPcutGenerationHeuristicCMIR(scip, NULL, POSTPROCESS, BOUNDSWITCH, VARTYPEUSEVBDS, allowlocal,
         MAXTESTDELTA, NULL, NULL, MINFRAC, MAXFRAC, aggrrow, cutcoefs, &cutrhs, cutinds, &cutnnz,
         &cutefficacy, &cutrank, &cutislocal, &success) );
   if( success && cutnnz > 0 && SCIPisEfficacious(scip, cutefficacy) )
   {
      besteff = cutefficacy;
      bestrhs = cutrhs;
      bestnnz = cutnnz;
      bestrank = cutrank;
      bestlocal = cutislocal;
      BMScopyMemoryArray(bestcoefs, cutcoefs, cutnnz);
      BMScopyMemoryArray(bestinds, cutinds, cutnnz);
      bestprefix = "pathcmir";
   }

   /* if cMIR did not find an efficacious cut, optionally try lifted flow cover cut */
   if( bestnnz == 0 && sepadata->tryflowcover )
   {
      SCIP_Real fceff = besteff;
      SCIP_Real fcrhs = 0.0;
      int fcnnz = 0;
      int fcrank = 0;
      SCIP_Bool fclocal = FALSE;
      SCIP_Bool fcsuccess = FALSE;

      SCIP_CALL( SCIPcalcFlowCover(scip, NULL, POSTPROCESS, BOUNDSWITCH, allowlocal, aggrrow,
            cutcoefs, &fcrhs, cutinds, &fcnnz, &fceff, &fcrank, &fclocal, &fcsuccess) );
      if( fcsuccess && fcnnz > 0 && fceff > besteff && SCIPisEfficacious(scip, fceff) )
      {
         besteff = fceff;
         bestrhs = fcrhs;
         bestnnz = fcnnz;
         bestrank = fcrank;
         bestlocal = fclocal;
         BMScopyMemoryArray(bestcoefs, cutcoefs, fcnnz);
         BMScopyMemoryArray(bestinds, cutinds, fcnnz);
         bestprefix = "pathfc";
      }
   }

   if( bestnnz > 0 )
   {
      SCIP_CALL( storeCutCandidate(scip, sepa, sepadata, vars, bestcoefs, bestinds, bestnnz, bestrhs,
            bestrank, bestlocal, bestprefix, cuts, ncuts, maxcuts, found) );
   }

   return SCIP_OKAY;
}

/** tries to separate a VBD-aware multi-row Günlük-Pochet path-mixing cut from cumulative valid <= path rows */
static
SCIP_RETCODE separatePathMixingCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_VAR**            vars,               /**< active problem variables */
   int                   nvars,              /**< number of active problem variables */
   const PATHVARINFO*    varinfo,            /**< precomputed variable metadata */
   SCIP_Bool             allowlocal,         /**< allow local bounds/cuts? */
   const PATHROWSNAP*    pathrows,           /**< valid <= cumulative path rows */
   int                   pathlen,            /**< number of path rows (>= 2) */
   SCIP_Real*            accvals,            /**< scratch dense accumulator (nvars, zeroed on entry/exit) */
   int*                  accinds,            /**< scratch active indices buffer (nvars) */
   int*                  accstamp,           /**< scratch stamp array (nvars) */
   int*                  curstamp,           /**< pointer to current stamp counter */
   SCIP_Real*            outcoefs,           /**< scratch output coefs buffer (nvars) */
   int*                  outinds,            /**< scratch output inds buffer (nvars) */
   SCIP_ROW**            cuts,               /**< candidate cuts array */
   int*                  ncuts,              /**< candidate count */
   int                   maxcuts,            /**< candidate capacity */
   SCIP_Bool*            found               /**< set to TRUE if a cut candidate was added */
   )
{
   SCIP_Real feastol;
   SCIP_Real epsilon;
   int* rowinds[MAXPATHLENHARD];
   SCIP_Real* rowvals[MAXPATHLENHARD];
   int rownnz[MAXPATHLENHARD];
   SCIP_Real trrhs[MAXPATHLENHARD];
   int rowperm[MAXPATHLENHARD];
   int nvalid = 0;
   int maxrank = 0;
   SCIP_Bool cutislocal = FALSE;
   int k;
   int j;

   assert(pathlen >= 2 && pathlen <= MAXPATHLENHARD);
   assert(found != NULL);
   *found = FALSE;

   if( *ncuts >= maxcuts )
      return SCIP_OKAY;

   feastol = SCIPfeastol(scip);
   epsilon = SCIPepsilon(scip);

   for( k = 0; k < pathlen; ++k )
   {
      rowinds[k] = NULL;
      rowvals[k] = NULL;
      rownnz[k] = 0;
   }

   /* Step 1: transform each cumulative path row into consistent nonnegative variables z_v >= 0,
    * substituting continuous variables using their preselected best VLB/VUB or simple bound
    */
   for( k = 0; k < pathlen; ++k )
   {
      long double r;
      SCIP_Bool rowok = TRUE;
      SCIP_Bool rowlocal = pathrows[k].islocal;
      int nacc = 0;
      int nnzkept = 0;
      int stamp;

      if( SCIPisInfinity(scip, REALABS(pathrows[k].rhs)) )
         continue;

      stamp = ++(*curstamp);
      r = (long double)pathrows[k].rhs;

      /* first pass: substitute continuous variables (may contribute to binary VBD variables) */
      for( j = 0; j < pathrows[k].nnz; ++j )
      {
         int v = pathrows[k].inds[j];
         SCIP_Real coef = pathrows[k].vals[j];
         const PATHVARINFO* vi = &varinfo[v];

         if( vi->pmbndtype == PMBND_NONE )
         {
            rowok = FALSE;
            break;
         }

         if( vi->iscont )
         {
            switch( vi->pmbndtype )
            {
               case PMBND_SIMPLE_LB:
                  r -= (long double)coef * (long double)vi->lb;
                  if( allowlocal && SCIPisGT(scip, vi->lb, SCIPvarGetLbGlobal(vars[v])) )
                     rowlocal = TRUE;
                  if( coef < 0.0 )
                  {
                     if( accstamp[v] != stamp )
                     {
                        accstamp[v] = stamp;
                        accinds[nacc++] = v;
                     }
                     accvals[v] += coef;
                  }
                  break;
               case PMBND_SIMPLE_UB:
                  r -= (long double)coef * (long double)vi->ub;
                  if( allowlocal && SCIPisLT(scip, vi->ub, SCIPvarGetUbGlobal(vars[v])) )
                     rowlocal = TRUE;
                  if( coef > 0.0 )
                  {
                     if( accstamp[v] != stamp )
                     {
                        accstamp[v] = stamp;
                        accinds[nacc++] = v;
                     }
                     accvals[v] -= coef;
                  }
                  break;
               case PMBND_VLB:
               {
                  int y = vi->vlbvaridx;
                  assert(y >= 0 && y < nvars);
                  r -= (long double)coef * (long double)vi->vlbconst;
                  if( accstamp[y] != stamp )
                  {
                     accstamp[y] = stamp;
                     accinds[nacc++] = y;
                  }
                  accvals[y] += coef * vi->vlbcoef;
                  if( coef < 0.0 )
                  {
                     if( accstamp[v] != stamp )
                     {
                        accstamp[v] = stamp;
                        accinds[nacc++] = v;
                     }
                     accvals[v] += coef;
                  }
                  break;
               }
               case PMBND_VUB:
               {
                  int y = vi->vubvaridx;
                  assert(y >= 0 && y < nvars);
                  r -= (long double)coef * (long double)vi->vubconst;
                  if( accstamp[y] != stamp )
                  {
                     accstamp[y] = stamp;
                     accinds[nacc++] = y;
                  }
                  accvals[y] += coef * vi->vubcoef;
                  if( coef > 0.0 )
                  {
                     if( accstamp[v] != stamp )
                     {
                        accstamp[v] = stamp;
                        accinds[nacc++] = v;
                     }
                     accvals[v] -= coef;
                  }
                  break;
               }
               default:
                  rowok = FALSE;
                  break;
            }
            if( !rowok )
               break;
         }
         else
         {
            /* accumulate integer variable in original space first (binary VBD contributions also land here) */
            if( accstamp[v] != stamp )
            {
               accstamp[v] = stamp;
               accinds[nacc++] = v;
            }
            accvals[v] += coef;
         }
      }

      if( !rowok )
      {
         for( j = 0; j < nacc; ++j )
            accvals[accinds[j]] = 0.0;
         continue;
      }

      /* second pass: shift/complement integer variables and extract nonzero transformed row */
      for( j = 0; j < nacc; ++j )
      {
         int v = accinds[j];
         SCIP_Real c = accvals[v];
         const PATHVARINFO* vi = &varinfo[v];

         if( c == 0.0 )
            continue;

         if( !vi->iscont )
         {
            if( vi->pmbndtype == PMBND_SIMPLE_LB )
            {
               r -= (long double)c * (long double)vi->lb;
               if( allowlocal && SCIPisGT(scip, vi->lb, SCIPvarGetLbGlobal(vars[v])) )
                  rowlocal = TRUE;
            }
            else if( vi->pmbndtype == PMBND_SIMPLE_UB )
            {
               r -= (long double)c * (long double)vi->ub;
               if( allowlocal && SCIPisLT(scip, vi->ub, SCIPvarGetUbGlobal(vars[v])) )
                  rowlocal = TRUE;
               c = -c;
               accvals[v] = c;
            }
            else
            {
               rowok = FALSE;
               break;
            }
         }

         /* safely drop tiny coefficients in nonnegative transformed space (0 <= x' <= u') */
         if( REALABS(c) <= epsilon )
         {
            if( c < 0.0 )
            {
               SCIP_Real glblb = SCIPvarGetLbGlobal(vars[v]);
               SCIP_Real glbub = SCIPvarGetUbGlobal(vars[v]);
               SCIP_Real loclb = allowlocal ? SCIPvarGetLbLocal(vars[v]) : glblb;
               SCIP_Real locub = allowlocal ? SCIPvarGetUbLocal(vars[v]) : glbub;

               /* relax c * x' >= c * (ub - lb) only when a finite simple domain range is known */
               if( !vi->iscont && !SCIPisInfinity(scip, -loclb) && !SCIPisInfinity(scip, locub) )
               {
                  r -= (long double)c * (long double)(locub - loclb);
                  if( allowlocal && (SCIPisGT(scip, loclb, glblb) || SCIPisLT(scip, locub, glbub)) )
                     rowlocal = TRUE;
               }
               else
               {
                  rowok = FALSE;
                  break;
               }
            }
            accvals[v] = 0.0;
            continue;
         }

         ++nnzkept;
      }

      if( !rowok || nnzkept == 0 || r > (long double)epsilon )
      {
         for( j = 0; j < nacc; ++j )
            accvals[accinds[j]] = 0.0;
         continue;
      }
      if( r >= -(long double)feastol )
         r = 0.0L;

      SCIP_CALL( SCIPallocBufferArray(scip, &rowinds[nvalid], nnzkept) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowvals[nvalid], nnzkept) );
      {
         int pos = 0;
         for( j = 0; j < nacc; ++j )
         {
            int v = accinds[j];
            SCIP_Real c = accvals[v];
            accvals[v] = 0.0;
            if( c != 0.0 )
            {
               rowinds[nvalid][pos] = v;
               rowvals[nvalid][pos] = c;
               ++pos;
            }
         }
         rownnz[nvalid] = pos;
      }

      if( rowlocal )
         cutislocal = TRUE;
      if( pathrows[k].rank > maxrank )
         maxrank = pathrows[k].rank;
      trrhs[nvalid] = (SCIP_Real)r;
      rowperm[nvalid] = nvalid;
      ++nvalid;
   }

   /* Step 2: sort transformed rows by decreasing RHS (0 >= rhs_0 > rhs_1 > ...) and run mixing */
   if( nvalid >= 2 )
   {
      int ordered[MAXPATHLENHARD];
      SCIP_Real orderedrhs[MAXPATHLENHARD];
      int nordered = 0;

      for( k = 0; k < nvalid; ++k )
      {
         int best = k;
         for( j = k + 1; j < nvalid; ++j )
         {
            if( trrhs[j] > trrhs[best] )
               best = j;
         }
         if( best != k )
         {
            SCIP_Real tr = trrhs[k];
            int tp = rowperm[k];
            trrhs[k] = trrhs[best];
            rowperm[k] = rowperm[best];
            trrhs[best] = tr;
            rowperm[best] = tp;
         }
      }

      for( k = 0; k < nvalid; ++k )
      {
         if( nordered == 0 || trrhs[k] < orderedrhs[nordered - 1] - feastol )
         {
            ordered[nordered] = rowperm[k];
            orderedrhs[nordered] = trrhs[k];
            ++nordered;
         }
      }

      if( nordered >= 2 )
      {
         int* varpos;
         int* uniquevars;
         SCIP_Real* trsolvals;
         int maxentries = 0;
         int nunique = 0;
         SCIP_Real delta = 1.0;

         for( k = 0; k < nordered; ++k )
            maxentries += rownnz[ordered[k]];

         SCIP_CALL( SCIPallocBufferArray(scip, &varpos, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &uniquevars, maxentries) );
         SCIP_CALL( SCIPallocBufferArray(scip, &trsolvals, maxentries) );

         for( k = 0; k < nordered; ++k )
         {
            int rk = ordered[k];
            for( j = 0; j < rownnz[rk]; ++j )
               varpos[rowinds[rk][j]] = -1;
         }

         for( k = 0; k < nordered; ++k )
         {
            int rk = ordered[k];
            SCIP_Real ar = REALABS(orderedrhs[k]);
            if( ar > delta )
               delta = ar;

            for( j = 0; j < rownnz[rk]; ++j )
            {
               int v = rowinds[rk][j];
               int pos = varpos[v];
               const PATHVARINFO* vi = &varinfo[v];

               if( pos < 0 )
               {
                  SCIP_Real sv = 0.0;
                  pos = nunique++;
                  varpos[v] = pos;
                  uniquevars[pos] = v;

                  switch( vi->pmbndtype )
                  {
                     case PMBND_SIMPLE_LB:
                        sv = MAX(0.0, vi->solval - vi->lb);
                        if( allowlocal && SCIPisGT(scip, vi->lb, SCIPvarGetLbGlobal(vars[v])) )
                           cutislocal = TRUE;
                        break;
                     case PMBND_SIMPLE_UB:
                        sv = MAX(0.0, vi->ub - vi->solval);
                        if( allowlocal && SCIPisLT(scip, vi->ub, SCIPvarGetUbGlobal(vars[v])) )
                           cutislocal = TRUE;
                        break;
                     case PMBND_VLB:
                        sv = MAX(0.0, vi->vlbdist);
                        break;
                     case PMBND_VUB:
                        sv = MAX(0.0, vi->vubdist);
                        break;
                     default:
                        sv = 0.0;
                        break;
                  }
                  trsolvals[pos] = sv;
               }

               if( !vi->iscont )
               {
                  SCIP_Real ac = REALABS(rowvals[rk][j]);
                  if( ac > delta )
                     delta = ac;
               }
            }
         }

         if( delta <= 1e12 && nunique > 0 )
         {
            SCIP_Real* cutvals;
            SCIP_Real* maxfrac;
            long double* downsum;
            long double* fsum;
            long double cutrhs = 0.0L;
            SCIP_Real flast = 0.0;
            SCIP_Real scale;

            delta = exp2(ceil(log2(delta + 1.0)));
            scale = -1.0 / delta;

            SCIP_CALL( SCIPallocClearBufferArray(scip, &cutvals, nunique) );
            SCIP_CALL( SCIPallocClearBufferArray(scip, &maxfrac, nunique) );
            SCIP_CALL( SCIPallocClearBufferArray(scip, &downsum, nunique) );
            SCIP_CALL( SCIPallocClearBufferArray(scip, &fsum, nunique) );

            for( k = 0; k < nordered; ++k )
            {
               int rk = ordered[k];
               SCIP_Real f = orderedrhs[k] * scale;
               long double fdiff = (long double)f - (long double)flast;
               cutrhs += fdiff;

               for( j = 0; j < rownnz[rk]; ++j )
               {
                  int v = rowinds[rk][j];
                  int pos = varpos[v];
                  SCIP_Real gj = rowvals[rk][j] * scale;

                  if( varinfo[v].iscont )
                  {
                     if( gj > cutvals[pos] )
                        cutvals[pos] = gj;
                  }
                  else
                  {
                     SCIP_Real gjdown = floor(gj);
                     SCIP_Real hj = gj - gjdown;
                     if( hj > maxfrac[pos] )
                        maxfrac[pos] = hj;
                     downsum[pos] += fdiff * (long double)gjdown;
                     fsum[pos] += fdiff;
                     if( fsum[pos] < (long double)maxfrac[pos] )
                        cutvals[pos] = (SCIP_Real)(downsum[pos] + fsum[pos]);
                     else
                        cutvals[pos] = (SCIP_Real)(downsum[pos] + (long double)maxfrac[pos]);
                  }
               }

               if( k > 0 )
               {
                  long double viol = cutrhs;
                  for( j = 0; j < nunique; ++j )
                     viol -= (long double)trsolvals[j] * (long double)cutvals[j];
                  viol *= (long double)delta;

                  if( viol > (long double)(10.0 * feastol) )
                  {
                     SCIP_Real unscale = -delta;
                     long double outrhs = cutrhs * (long double)unscale;
                     SCIP_Bool abortcut = FALSE;
                     int nacc = 0;
                     int outnnz = 0;
                     int stamp = ++(*curstamp);

                     /* untransform back to original problem variables */
                     for( j = 0; j < nunique; ++j )
                     {
                        int v = uniquevars[j];
                        SCIP_Real c = cutvals[j] * unscale;
                        const PATHVARINFO* vi = &varinfo[v];

                        if( c == 0.0 )
                           continue;

                        switch( vi->pmbndtype )
                        {
                           case PMBND_SIMPLE_LB:
                              outrhs += (long double)c * (long double)vi->lb;
                              if( accstamp[v] != stamp )
                              {
                                 accstamp[v] = stamp;
                                 accinds[nacc++] = v;
                              }
                              accvals[v] += c;
                              break;
                           case PMBND_SIMPLE_UB:
                              outrhs -= (long double)c * (long double)vi->ub;
                              if( accstamp[v] != stamp )
                              {
                                 accstamp[v] = stamp;
                                 accinds[nacc++] = v;
                              }
                              accvals[v] -= c;
                              break;
                           case PMBND_VLB:
                           {
                              int y = vi->vlbvaridx;
                              outrhs += (long double)c * (long double)vi->vlbconst;
                              if( accstamp[v] != stamp )
                              {
                                 accstamp[v] = stamp;
                                 accinds[nacc++] = v;
                              }
                              accvals[v] += c;
                              if( accstamp[y] != stamp )
                              {
                                 accstamp[y] = stamp;
                                 accinds[nacc++] = y;
                              }
                              accvals[y] -= c * vi->vlbcoef;
                              break;
                           }
                           case PMBND_VUB:
                           {
                              int y = vi->vubvaridx;
                              outrhs -= (long double)c * (long double)vi->vubconst;
                              if( accstamp[v] != stamp )
                              {
                                 accstamp[v] = stamp;
                                 accinds[nacc++] = v;
                              }
                              accvals[v] -= c;
                              if( accstamp[y] != stamp )
                              {
                                 accstamp[y] = stamp;
                                 accinds[nacc++] = y;
                              }
                              accvals[y] += c * vi->vubcoef;
                              break;
                           }
                           default:
                              abortcut = TRUE;
                              break;
                        }
                        if( abortcut )
                           break;
                     }

                     /* safely clean tiny coefficients using finite bounds (never drop unbounded vars) */
                     for( j = 0; j < nacc; ++j )
                     {
                        int v = accinds[j];
                        SCIP_Real c = accvals[v];
                        accvals[v] = 0.0;

                        if( abortcut || c == 0.0 )
                           continue;

                        if( REALABS(c) <= feastol )
                        {
                           SCIP_Real lb = allowlocal ? SCIPvarGetLbLocal(vars[v]) : SCIPvarGetLbGlobal(vars[v]);
                           SCIP_Real ub = allowlocal ? SCIPvarGetUbLocal(vars[v]) : SCIPvarGetUbGlobal(vars[v]);

                           if( c > 0.0 )
                           {
                              if( SCIPisInfinity(scip, -lb) )
                                 abortcut = TRUE;
                              else
                              {
                                 outrhs -= (long double)c * (long double)lb;
                                 if( allowlocal && SCIPisGT(scip, lb, SCIPvarGetLbGlobal(vars[v])) )
                                    cutislocal = TRUE;
                              }
                           }
                           else
                           {
                              if( SCIPisInfinity(scip, ub) )
                                 abortcut = TRUE;
                              else
                              {
                                 outrhs -= (long double)c * (long double)ub;
                                 if( allowlocal && SCIPisLT(scip, ub, SCIPvarGetUbGlobal(vars[v])) )
                                    cutislocal = TRUE;
                              }
                           }
                           continue;
                        }

                        outinds[outnnz] = v;
                        outcoefs[outnnz] = c;
                        ++outnnz;
                     }

                     if( !abortcut && outnnz > 0 && !SCIPisInfinity(scip, REALABS((SCIP_Real)outrhs)) )
                     {
                        SCIP_CALL( storeCutCandidate(scip, sepa, sepadata, vars, outcoefs, outinds, outnnz,
                              (SCIP_Real)outrhs, maxrank + 1, cutislocal, "pathmix", cuts, ncuts, maxcuts, found) );
                        if( *found )
                           break;
                     }
                  }
               }
               flast = f;
            }

            SCIPfreeBufferArray(scip, &fsum);
            SCIPfreeBufferArray(scip, &downsum);
            SCIPfreeBufferArray(scip, &maxfrac);
            SCIPfreeBufferArray(scip, &cutvals);
         }

         SCIPfreeBufferArray(scip, &trsolvals);
         SCIPfreeBufferArray(scip, &uniquevars);
         SCIPfreeBufferArray(scip, &varpos);
      }
   }

   for( k = nvalid - 1; k >= 0; --k )
   {
      SCIPfreeBufferArray(scip, &rowvals[k]);
      SCIPfreeBufferArray(scip, &rowinds[k]);
   }

   return SCIP_OKAY;
}

/** checks whether a row index is already on the current path */
static
SCIP_Bool isRowOnPath(
   const int*            currpath,           /**< current path row indices */
   int                   currpathlen,        /**< current path length */
   int                   row                 /**< row index to test */
   )
{
   int i;

   for( i = 0; i < currpathlen; ++i )
   {
      if( currpath[i] == row )
         return TRUE;
   }
   return FALSE;
}

/** checks whether a variable has any available arc row outside the current path */
static
SCIP_Bool hasAvailableArc(
   const PATHARC*        arcs,               /**< arc array */
   int                   start,              /**< start offset for variable */
   int                   end,                /**< end offset for variable */
   const int*            currpath,           /**< current path rows */
   int                   currpathlen         /**< current path length */
   )
{
   int k;

   if( start >= end )
      return FALSE;
   if( end - start > currpathlen )
      return TRUE;

   for( k = start; k < end; ++k )
   {
      if( !isRowOnPath(currpath, currpathlen, arcs[k].row) )
         return TRUE;
   }
   return FALSE;
}

/** selects a random unused arc row with numerically safe aggregation weight */
static
SCIP_Bool pickArcRow(
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   const PATHARC*        arcs,               /**< arc array */
   int                   start,              /**< start offset */
   int                   end,                /**< end offset */
   SCIP_Real             val,                /**< current coefficient of variable to cancel */
   SCIP_Real             minweight,          /**< minimal allowed |weight| */
   SCIP_Real             maxweight,          /**< maximal allowed |weight| */
   const int*            currpath,           /**< current path rows */
   int                   currpathlen,        /**< current path length */
   int*                  outrow,             /**< selected LP row index */
   SCIP_Real*            outweight           /**< selected aggregation weight */
   )
{
   int len = end - start;
   int pick;
   int k;

   if( len <= 0 )
      return FALSE;

   pick = start + SCIPrandomGetInt(randnumgen, 0, len - 1);

   for( k = pick; k < end; ++k )
   {
      int r = arcs[k].row;
      SCIP_Real w = -val / arcs[k].coef;
      SCIP_Real aw = REALABS(w);
      if( !isRowOnPath(currpath, currpathlen, r) && aw >= minweight && aw <= maxweight )
      {
         *outrow = r;
         *outweight = w;
         return TRUE;
      }
   }

   for( k = start; k < pick; ++k )
   {
      int r = arcs[k].row;
      SCIP_Real w = -val / arcs[k].coef;
      SCIP_Real aw = REALABS(w);
      if( !isRowOnPath(currpath, currpathlen, r) && aw >= minweight && aw <= maxweight )
      {
         *outrow = r;
         *outweight = w;
         return TRUE;
      }
   }

   return FALSE;
}

/** main separation routine */
static
SCIP_RETCODE separatePathCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_Bool             allowlocal,         /**< allow local cuts? */
   int                   depth,              /**< current tree depth */
   SCIP_RESULT*          result              /**< pointer to store result */
   )
{
   SCIP_VAR** vars;
   SCIP_ROW** rows;
   PATHVARINFO* varinfo = NULL;
   signed char* rowtype = NULL;
   int* numcont = NULL;
   int* inarcstart = NULL;
   int* inarcend = NULL;
   int* outarcstart = NULL;
   int* outarcend = NULL;
   PATHARC* inarcs = NULL;
   PATHARC* outarcs = NULL;
   SCIP_AGGRROW* aggrrow = NULL;
   SCIP_Real* cutcoefs = NULL;
   int* cutinds = NULL;
   SCIP_Real* bestcoefs = NULL;
   int* bestinds = NULL;
   SCIP_Real* accvals = NULL;
   int* accinds = NULL;
   int* accstamp = NULL;
   int curstamp = 0;
   SCIP_ROW** cuts = NULL;
   SCIP_Real feastol;
   SCIP_Real epsilon;
   SCIP_Real minweight;
   SCIP_Real maxweight;
   SCIP_Bool cutoff = FALSE;
   int nvars;
   int nrows;
   int maxsepacuts;
   int maxtries;
   int maxaggrnonzs;
   int maxpathlen;
   int nactivecont = 0;
   int ninarcs = 0;
   int noutarcs = 0;
   int ncuts = 0;
   int maxcuts;
   int ntries = 0;
   SCIP_Longint roundwork = 0;
   int r;
   int v;
   int j;

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   if( nrows == 0 || nvars == 0 )
      return SCIP_OKAY;

   maxsepacuts = (depth == 0) ? sepadata->maxsepacutsroot : sepadata->maxsepacuts;
   maxtries = (depth == 0) ? sepadata->maxtriesroot : sepadata->maxtries;
   maxcuts = sepadata->maxcutcands;
   maxpathlen = MIN(sepadata->maxpathlen, MAXPATHLENHARD);

   if( maxsepacuts <= 0 || maxcuts <= 0 || maxtries <= 0 || maxpathlen < 2 )
      return SCIP_OKAY;

   feastol = SCIPfeastol(scip);
   epsilon = SCIPepsilon(scip);
   minweight = feastol;
   maxweight = 1.0 / feastol;
   maxaggrnonzs = (int)(sepadata->maxaggdensity * (SCIP_Real)nvars) + sepadata->densityoffset;
   if( maxaggrnonzs < 50 )
      maxaggrnonzs = 50;
   if( maxaggrnonzs > MAXAGGRNNZHARD )
      maxaggrnonzs = MAXAGGRNNZHARD;

   SCIP_CALL( SCIPallocBufferArray(scip, &varinfo, nvars) );

   /* Step 1: compute bound distances and best binary VLB/VUB for all variables */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var = vars[v];
      PATHVARINFO* vi = &varinfo[v];
      SCIP_Real x = SCIPvarGetLPSol(var);
      SCIP_Real lb = allowlocal ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);
      SCIP_Real ub = allowlocal ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);
      SCIP_Bool finlb = !SCIPisInfinity(scip, -lb);
      SCIP_Bool finub = !SCIPisInfinity(scip, ub);
      SCIP_Real slbdist = finlb ? MAX(0.0, x - lb) : SCIP_REAL_MAX;
      SCIP_Real subdist = finub ? MAX(0.0, ub - x) : SCIP_REAL_MAX;

      vi->solval = x;
      vi->lb = lb;
      vi->ub = ub;
      vi->iscont = !SCIPvarIsIntegral(var);
      vi->vlbvaridx = -1;
      vi->vubvaridx = -1;
      vi->vlbcoef = 0.0;
      vi->vlbconst = 0.0;
      vi->vubcoef = 0.0;
      vi->vubconst = 0.0;
      vi->vlbdist = SCIP_REAL_MAX;
      vi->vubdist = SCIP_REAL_MAX;
      vi->substrow = -1;
      vi->substcoef = 0.0;
      vi->pmbndtype = PMBND_NONE;

      if( slbdist <= feastol )
         slbdist = 0.0;
      if( subdist <= feastol )
         subdist = 0.0;

      if( !vi->iscont )
      {
         vi->bounddist = MIN(slbdist, subdist);
         if( finlb && finub )
            vi->pmbndtype = (slbdist <= subdist) ? PMBND_SIMPLE_LB : PMBND_SIMPLE_UB;
         else if( finlb )
            vi->pmbndtype = PMBND_SIMPLE_LB;
         else if( finub )
            vi->pmbndtype = PMBND_SIMPLE_UB;
      }
      else
      {
         SCIP_Real bestlbdist = slbdist;
         SCIP_Real bestubdist = subdist;
         SCIP_Real closestvlb = -SCIPinfinity(scip);
         SCIP_Real closestvub = SCIPinfinity(scip);
         int closestvlbidx = -1;
         int closestvubidx = -1;

         SCIP_CALL( SCIPgetVarClosestVlb(scip, var, NULL, &closestvlb, &closestvlbidx) );
         if( closestvlbidx >= 0 )
         {
            SCIP_VAR* vlbvar = SCIPvarGetVlbVars(var)[closestvlbidx];
            SCIP_Real vlbcoef = SCIPvarGetVlbCoefs(var)[closestvlbidx];
            SCIP_Real vlbconst = SCIPvarGetVlbConstants(var)[closestvlbidx];
            if( SCIPvarIsActive(vlbvar) && SCIPvarIsBinary(vlbvar)
               && SCIPvarGetProbindex(vlbvar) >= 0
               && SCIPvarGetLbGlobal(vlbvar) == 0.0 && SCIPvarGetUbGlobal(vlbvar) == 1.0
               && !SCIPisInfinity(scip, REALABS(vlbcoef)) && !SCIPisInfinity(scip, REALABS(vlbconst)) )
            {
               SCIP_Real vdist = MAX(0.0, x - closestvlb);
               if( vdist <= feastol )
                  vdist = 0.0;
               if( vdist <= bestlbdist + feastol )
               {
                  vi->vlbvaridx = SCIPvarGetProbindex(vlbvar);
                  vi->vlbcoef = vlbcoef;
                  vi->vlbconst = vlbconst;
                  vi->vlbdist = vdist;
                  bestlbdist = vdist;
               }
            }
         }

         SCIP_CALL( SCIPgetVarClosestVub(scip, var, NULL, &closestvub, &closestvubidx) );
         if( closestvubidx >= 0 )
         {
            SCIP_VAR* vubvar = SCIPvarGetVubVars(var)[closestvubidx];
            SCIP_Real vubcoef = SCIPvarGetVubCoefs(var)[closestvubidx];
            SCIP_Real vubconst = SCIPvarGetVubConstants(var)[closestvubidx];
            if( SCIPvarIsActive(vubvar) && SCIPvarIsBinary(vubvar)
               && SCIPvarGetProbindex(vubvar) >= 0
               && SCIPvarGetLbGlobal(vubvar) == 0.0 && SCIPvarGetUbGlobal(vubvar) == 1.0
               && !SCIPisInfinity(scip, REALABS(vubcoef)) && !SCIPisInfinity(scip, REALABS(vubconst)) )
            {
               SCIP_Real vdist = MAX(0.0, closestvub - x);
               if( vdist <= feastol )
                  vdist = 0.0;
               if( vdist <= bestubdist + feastol )
               {
                  vi->vubvaridx = SCIPvarGetProbindex(vubvar);
                  vi->vubcoef = vubcoef;
                  vi->vubconst = vubconst;
                  vi->vubdist = vdist;
                  bestubdist = vdist;
               }
            }
         }

         vi->bounddist = MIN(bestlbdist, bestubdist);
         if( vi->bounddist <= feastol )
            vi->bounddist = 0.0;

         /* pick ONE consistent bound type for multi-row path-mixing */
         if( bestlbdist < SCIP_REAL_MAX || bestubdist < SCIP_REAL_MAX )
         {
            if( bestlbdist <= bestubdist )
               vi->pmbndtype = (vi->vlbvaridx >= 0) ? PMBND_VLB : PMBND_SIMPLE_LB;
            else
               vi->pmbndtype = (vi->vubvaridx >= 0) ? PMBND_VUB : PMBND_SIMPLE_UB;
         }

         if( vi->bounddist > 0.0 )
            ++nactivecont;
      }
   }

   /* fast exit on pure-integer problems or when all continuous vars are at bounds */
   if( nactivecont == 0 )
   {
      SCIPfreeBufferArray(scip, &varinfo);
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPallocBufferArray(scip, &rowtype, nrows) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &numcont, nrows) );

   /* Step 2: classify tight LP rows and count active continuous variables per row */
   for( r = 0; r < nrows; ++r )
   {
      SCIP_ROW* row = rows[r];
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real act;
      SCIP_Real lslack;
      SCIP_Real uslack;
      int nlpnonz;

      rowtype[r] = ROWTYPE_UNUSABLE;

      if( SCIProwIsModifiable(row) || (!allowlocal && SCIProwIsLocal(row)) )
         continue;

      nlpnonz = SCIProwGetNLPNonz(row);
      if( nlpnonz == 0 || SCIProwGetNNonz(row) != nlpnonz || nlpnonz > maxaggrnonzs )
         continue;

      lhs = SCIProwGetLhs(row);
      rhs = SCIProwGetRhs(row);
      act = SCIPgetRowLPActivity(scip, row);

      if( !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && SCIPisEQ(scip, lhs, rhs) )
      {
         rowtype[r] = ROWTYPE_EQ;
      }
      else
      {
         lslack = !SCIPisInfinity(scip, -lhs) ? (act - lhs) : SCIP_REAL_MAX;
         uslack = !SCIPisInfinity(scip, rhs) ? (rhs - act) : SCIP_REAL_MAX;

         if( lslack > feastol && uslack > feastol )
            continue;
         else if( lslack < uslack )
            rowtype[r] = ROWTYPE_GEQ;
         else
            rowtype[r] = ROWTYPE_LEQ;
      }

      {
         SCIP_COL** rcols = SCIProwGetCols(row);
         int ccnt = 0;
         for( j = 0; j < nlpnonz; ++j )
         {
            int pidx = SCIPvarGetProbindex(SCIPcolGetVar(rcols[j]));
            if( pidx >= 0 && varinfo[pidx].iscont && varinfo[pidx].bounddist > 0.0 )
               ++ccnt;
         }
         numcont[r] = ccnt;
      }
   }

   /* Step 3: register single-continuous-var equality rows as automatic substitutions */
   for( r = 0; r < nrows; ++r )
   {
      SCIP_COL** rcols;
      SCIP_Real* rvals;
      int nlpnonz;
      int contvar = -1;
      SCIP_Real contcoef = 0.0;

      if( rowtype[r] != ROWTYPE_EQ || numcont[r] != 1 )
         continue;

      rcols = SCIProwGetCols(rows[r]);
      rvals = SCIProwGetVals(rows[r]);
      nlpnonz = SCIProwGetNLPNonz(rows[r]);

      for( j = 0; j < nlpnonz; ++j )
      {
         int pidx = SCIPvarGetProbindex(SCIPcolGetVar(rcols[j]));
         if( pidx >= 0 && varinfo[pidx].iscont && varinfo[pidx].bounddist > 0.0 )
         {
            contvar = pidx;
            contcoef = rvals[j];
            break;
         }
      }

      if( contvar >= 0 && varinfo[contvar].substrow == -1 && REALABS(contcoef) > epsilon )
      {
         varinfo[contvar].substrow = r;
         varinfo[contvar].substcoef = contcoef;
         rowtype[r] = ROWTYPE_UNUSABLE;
      }
   }

   /* Step 4: build directed CSR inArcRows and outArcRows for non-substituted active continuous vars */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &inarcstart, nvars + 1) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &outarcstart, nvars + 1) );

   for( r = 0; r < nrows; ++r )
   {
      SCIP_COL** rcols;
      SCIP_Real* rvals;
      int nlpnonz;
      signed char rt = rowtype[r];

      if( rt == ROWTYPE_UNUSABLE || numcont[r] == 0 )
         continue;

      rcols = SCIProwGetCols(rows[r]);
      rvals = SCIProwGetVals(rows[r]);
      nlpnonz = SCIProwGetNLPNonz(rows[r]);

      for( j = 0; j < nlpnonz; ++j )
      {
         int pidx = SCIPvarGetProbindex(SCIPcolGetVar(rcols[j]));
         SCIP_Real coef = rvals[j];

         if( pidx < 0 || !varinfo[pidx].iscont || varinfo[pidx].bounddist <= 0.0 || varinfo[pidx].substrow >= 0 )
            continue;
         if( REALABS(coef) <= epsilon )
            continue;

         if( rt == ROWTYPE_LEQ )
         {
            if( coef < 0.0 )
               ++inarcstart[pidx + 1];
            else
               ++outarcstart[pidx + 1];
         }
         else if( rt == ROWTYPE_GEQ )
         {
            if( coef > 0.0 )
               ++inarcstart[pidx + 1];
            else
               ++outarcstart[pidx + 1];
         }
         else
         {
            ++inarcstart[pidx + 1];
            ++outarcstart[pidx + 1];
         }
      }
   }

   for( v = 0; v < nvars; ++v )
   {
      inarcstart[v + 1] += inarcstart[v];
      outarcstart[v + 1] += outarcstart[v];
   }
   ninarcs = inarcstart[nvars];
   noutarcs = outarcstart[nvars];

   SCIP_CALL( SCIPallocBufferArray(scip, &inarcend, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &outarcend, nvars) );
   BMScopyMemoryArray(inarcend, inarcstart, nvars);
   BMScopyMemoryArray(outarcend, outarcstart, nvars);

   SCIP_CALL( SCIPallocBufferArray(scip, &inarcs, MAX(ninarcs, 1)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &outarcs, MAX(noutarcs, 1)) );

   for( r = 0; r < nrows; ++r )
   {
      SCIP_COL** rcols;
      SCIP_Real* rvals;
      int nlpnonz;
      signed char rt = rowtype[r];

      if( rt == ROWTYPE_UNUSABLE || numcont[r] == 0 )
         continue;

      rcols = SCIProwGetCols(rows[r]);
      rvals = SCIProwGetVals(rows[r]);
      nlpnonz = SCIProwGetNLPNonz(rows[r]);

      for( j = 0; j < nlpnonz; ++j )
      {
         int pidx = SCIPvarGetProbindex(SCIPcolGetVar(rcols[j]));
         SCIP_Real coef = rvals[j];

         if( pidx < 0 || !varinfo[pidx].iscont || varinfo[pidx].bounddist <= 0.0 || varinfo[pidx].substrow >= 0 )
            continue;
         if( REALABS(coef) <= epsilon )
            continue;

         if( rt == ROWTYPE_LEQ )
         {
            if( coef < 0.0 )
            {
               int pos = inarcend[pidx]++;
               inarcs[pos].row = r;
               inarcs[pos].coef = coef;
            }
            else
            {
               int pos = outarcend[pidx]++;
               outarcs[pos].row = r;
               outarcs[pos].coef = coef;
            }
         }
         else if( rt == ROWTYPE_GEQ )
         {
            if( coef > 0.0 )
            {
               int pos = inarcend[pidx]++;
               inarcs[pos].row = r;
               inarcs[pos].coef = coef;
            }
            else
            {
               int pos = outarcend[pidx]++;
               outarcs[pos].row = r;
               outarcs[pos].coef = coef;
            }
         }
         else
         {
            int ipos = inarcend[pidx]++;
            int opos = outarcend[pidx]++;
            inarcs[ipos].row = r;
            inarcs[ipos].coef = coef;
            outarcs[opos].row = r;
            outarcs[opos].coef = coef;
         }
      }
   }

   /* Step 5: walk directed paths from every usable tight row with active continuous variables */
   SCIP_CALL( SCIPaggrRowCreate(scip, &aggrrow) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutinds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestinds, nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &accvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &accinds, nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &accstamp, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cuts, maxcuts) );

   for( r = 0; r < nrows && ntries < maxtries && ncuts < maxcuts && roundwork < MAXWORKPERROUND && !SCIPisStopped(scip); ++r )
   {
      SCIP_Real scales[2];
      int nscales;
      int s;
      signed char rt = rowtype[r];

      if( rt == ROWTYPE_UNUSABLE || numcont[r] == 0 )
         continue;

      ++ntries;

      if( rt == ROWTYPE_EQ )
      {
         SCIP_Real dualsol = SCIProwGetDualsol(rows[r]);
         nscales = 2;
         if( dualsol <= epsilon )
         {
            scales[0] = 1.0;
            scales[1] = -1.0;
         }
         else
         {
            scales[0] = -1.0;
            scales[1] = 1.0;
         }
      }
      else if( rt == ROWTYPE_LEQ )
      {
         nscales = 1;
         scales[0] = 1.0;
      }
      else
      {
         nscales = 1;
         scales[0] = -1.0;
      }

      for( s = 0; s < nscales && ncuts < maxcuts && !SCIPisStopped(scip); ++s )
      {
         int currpath[MAXPATHLENHARD];
         PATHROWSNAP pathsnaps[MAXPATHLENHARD];
         int currpathlen = 1;
         int nsnaps = 0;
         int subpasses = 0;
         SCIP_Bool foundcut = FALSE;
         SCIP_Bool didsub = FALSE;
         int p;

         SCIPaggrRowClear(aggrrow);
         SCIP_CALL( SCIPaggrRowAddRow(scip, aggrrow, rows[r], scales[s], rt == ROWTYPE_EQ ? 0 : rt) );
         currpath[0] = r;

         while( currpathlen <= maxpathlen && !SCIPisStopped(scip) )
         {
            int* aggrinds;
            int aggrnnz;
            SCIP_Bool addedsub = FALSE;
            int bestincol = -1;
            SCIP_Real bestinval = 0.0;
            SCIP_Real bestindist = 0.0;
            int bestoutcol = -1;
            SCIP_Real bestoutval = 0.0;
            SCIP_Real bestoutdist = 0.0;

         SCIP_Bool aggrvalid = TRUE;

         SCIPaggrRowRemoveZeros(scip, aggrrow, FALSE, &aggrvalid);
         if( !aggrvalid )
            break;
         aggrnnz = SCIPaggrRowGetNNz(aggrrow);
         if( aggrnnz == 0 || aggrnnz > maxaggrnonzs )
            break;

         aggrinds = SCIPaggrRowGetInds(aggrrow);

            /* first substitute any single-continuous-var equality columns */
            if( subpasses < MAXSUBPASSES )
            {
               for( j = 0; j < aggrnnz; ++j )
               {
                  int c = aggrinds[j];
                  if( varinfo[c].substrow >= 0 )
                  {
                     SCIP_Real cval = SCIPaggrRowGetProbvarValue(aggrrow, c);
                     SCIP_Real w = -cval / varinfo[c].substcoef;
                     SCIP_Real aw = REALABS(w);
                     if( aw >= minweight && aw <= maxweight )
                     {
                        SCIP_CALL( SCIPaggrRowAddRow(scip, aggrrow, rows[varinfo[c].substrow], w, 0) );
                        addedsub = TRUE;
                        didsub = TRUE;
                        ++subpasses;
                        break;
                     }
                  }
               }
            }
            if( addedsub )
               continue;

            /* select best continuous variables to cancel along directed arcs */
            for( j = 0; j < aggrnnz; ++j )
            {
               int c = aggrinds[j];
               SCIP_Real cval;

               if( !varinfo[c].iscont || varinfo[c].bounddist <= 0.0 || varinfo[c].substrow >= 0 )
                  continue;

               cval = SCIPaggrRowGetProbvarValue(aggrrow, c);
               if( REALABS(cval) <= epsilon )
                  continue;

               if( cval > 0.0 )
               {
                  if( inarcend[c] - inarcstart[c] <= MAXARCDEG
                     && hasAvailableArc(inarcs, inarcstart[c], inarcend[c], currpath, currpathlen) )
                  {
                     if( bestincol == -1 || varinfo[c].bounddist > bestindist )
                     {
                        bestincol = c;
                        bestinval = cval;
                        bestindist = varinfo[c].bounddist;
                     }
                  }
               }
               else
               {
                  if( outarcend[c] - outarcstart[c] <= MAXARCDEG
                     && hasAvailableArc(outarcs, outarcstart[c], outarcend[c], currpath, currpathlen) )
                  {
                     if( bestoutcol == -1 || varinfo[c].bounddist > bestoutdist )
                     {
                        bestoutcol = c;
                        bestoutval = cval;
                        bestoutdist = varinfo[c].bounddist;
                     }
                  }
               }
            }

            /* snapshot current valid <= row for multi-row path-mixing (if we will extend or already have snaps) */
            if( sepadata->usepathmix && nsnaps < maxpathlen && (nsnaps > 0 || bestincol >= 0 || bestoutcol >= 0) )
            {
               SCIP_CALL( SCIPallocBufferArray(scip, &pathsnaps[nsnaps].inds, aggrnnz) );
               SCIP_CALL( SCIPallocBufferArray(scip, &pathsnaps[nsnaps].vals, aggrnnz) );
               for( j = 0; j < aggrnnz; ++j )
               {
                  pathsnaps[nsnaps].inds[j] = aggrinds[j];
                  pathsnaps[nsnaps].vals[j] = SCIPaggrRowGetProbvarValue(aggrrow, aggrinds[j]);
               }
               pathsnaps[nsnaps].nnz = aggrnnz;
               pathsnaps[nsnaps].rhs = SCIPaggrRowGetRhs(aggrrow);
               pathsnaps[nsnaps].rank = SCIPaggrRowGetRank(aggrrow);
               pathsnaps[nsnaps].islocal = SCIPaggrRowIsLocal(aggrrow);
               ++nsnaps;
            }

            /* try single-row cuts once at least one aggregation/substitution occurred */
            if( currpathlen >= 2 || didsub )
            {
               roundwork += (SCIP_Longint)aggrnnz * (MAXTESTDELTA + 1);
               SCIP_CALL( trySingleRowCut(scip, sepa, sepadata, vars, aggrrow, allowlocal,
                     cutcoefs, cutinds, bestcoefs, bestinds, cuts, &ncuts, maxcuts, &foundcut) );
               if( foundcut || roundwork >= MAXWORKPERROUND )
                  break;
            }

            if( currpathlen == maxpathlen || (bestincol == -1 && bestoutcol == -1) )
               break;

            /* prefer out-arc direction on ties within feastol to keep path oriented and diverse */
            {
               int nextrow = -1;
               SCIP_Real nextw = 0.0;
               SCIP_Bool picked = FALSE;

               if( bestincol == -1 || (bestoutcol != -1 && bestoutdist >= bestindist - feastol) )
               {
                  picked = pickArcRow(sepadata->randnumgen, outarcs, outarcstart[bestoutcol], outarcend[bestoutcol],
                        bestoutval, minweight, maxweight, currpath, currpathlen, &nextrow, &nextw);
                  if( !picked && bestincol != -1 )
                  {
                     picked = pickArcRow(sepadata->randnumgen, inarcs, inarcstart[bestincol], inarcend[bestincol],
                           bestinval, minweight, maxweight, currpath, currpathlen, &nextrow, &nextw);
                  }
               }
               else
               {
                  picked = pickArcRow(sepadata->randnumgen, inarcs, inarcstart[bestincol], inarcend[bestincol],
                        bestinval, minweight, maxweight, currpath, currpathlen, &nextrow, &nextw);
                  if( !picked && bestoutcol != -1 )
                  {
                     picked = pickArcRow(sepadata->randnumgen, outarcs, outarcstart[bestoutcol], outarcend[bestoutcol],
                           bestoutval, minweight, maxweight, currpath, currpathlen, &nextrow, &nextw);
                  }
               }

               if( !picked )
                  break;

               assert(nextrow >= 0 && nextrow < nrows);
               SCIP_CALL( SCIPaggrRowAddRow(scip, aggrrow, rows[nextrow], nextw,
                     rowtype[nextrow] == ROWTYPE_EQ ? 0 : rowtype[nextrow]) );
               currpath[currpathlen++] = nextrow;
            }
         }

         /* if no single-row cut succeeded along the path, try VBD-aware multi-row path-mixing */
         if( !foundcut && sepadata->usepathmix && nsnaps >= 2 && ncuts < maxcuts )
         {
            SCIP_CALL( separatePathMixingCut(scip, sepa, sepadata, vars, nvars, varinfo, allowlocal,
                  pathsnaps, nsnaps, accvals, accinds, accstamp, &curstamp, bestcoefs, bestinds,
                  cuts, &ncuts, maxcuts, &foundcut) );
         }

         for( p = nsnaps - 1; p >= 0; --p )
         {
            SCIPfreeBufferArray(scip, &pathsnaps[p].vals);
            SCIPfreeBufferArray(scip, &pathsnaps[p].inds);
         }

         /* if first equality orientation succeeded, no need to repeat reversed orientation */
         if( foundcut )
            break;
      }
   }

   /* Step 6: select best cuts with hybrid cut selector */
   if( ncuts > 0 )
   {
      int nselected = 0;
      int i;

      SCIP_CALL( SCIPselectCutsHybrid(scip, cuts, NULL, sepadata->randnumgen, 1.0, 0.5, 0.1, 0.1, 0.0, 1.0,
            0.0, 0.0, ncuts, 0, maxsepacuts, &nselected) );

      for( i = 0; i < ncuts; ++i )
      {
         if( i < nselected && !cutoff )
         {
            if( SCIProwIsLocal(cuts[i]) )
            {
               SCIP_CALL( SCIPaddRow(scip, cuts[i], FALSE, &cutoff) );
            }
            else
            {
               SCIP_CALL( SCIPaddPoolCut(scip, cuts[i]) );
            }
         }
         SCIP_CALL( SCIPreleaseRow(scip, &cuts[i]) );
      }

      if( cutoff )
         *result = SCIP_CUTOFF;
      else if( nselected > 0 )
         *result = SCIP_SEPARATED;
   }

   SCIPfreeBufferArray(scip, &cuts);
   SCIPfreeBufferArray(scip, &accstamp);
   SCIPfreeBufferArray(scip, &accinds);
   SCIPfreeBufferArray(scip, &accvals);
   SCIPfreeBufferArray(scip, &bestinds);
   SCIPfreeBufferArray(scip, &bestcoefs);
   SCIPfreeBufferArray(scip, &cutinds);
   SCIPfreeBufferArray(scip, &cutcoefs);
   SCIPaggrRowFree(scip, &aggrrow);

   SCIPfreeBufferArray(scip, &outarcs);
   SCIPfreeBufferArray(scip, &inarcs);
   SCIPfreeBufferArray(scip, &outarcend);
   SCIPfreeBufferArray(scip, &inarcend);
   SCIPfreeBufferArray(scip, &outarcstart);
   SCIPfreeBufferArray(scip, &inarcstart);
   SCIPfreeBufferArray(scip, &numcont);
   SCIPfreeBufferArray(scip, &rowtype);
   SCIPfreeBufferArray(scip, &varinfo);

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

/** copy method for separator plugins */
static
SCIP_DECL_SEPACOPY(sepaCopyPath)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   SCIP_CALL( SCIPincludeSepaPath(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data */
static
SCIP_DECL_SEPAFREE(sepaFreePath)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeBlockMemory(scip, &sepadata);
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of separator */
static
SCIP_DECL_SEPAINITSOL(sepaInitsolPath)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIP_CALL( SCIPcreateRandom(scip, &sepadata->randnumgen, (unsigned int)sepadata->initseed, TRUE) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of separator */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolPath)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeRandom(scip, &sepadata->randnumgen);

   return SCIP_OKAY;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpPath)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   int ncalls;

   assert(result != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   if( SCIPisStopped(scip) || SCIPisExact(scip) )
      return SCIP_OKAY;

   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   ncalls = SCIPsepaGetNCallsAtNode(sepa);
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   SCIP_CALL( separatePathCuts(scip, sepa, sepadata, allowlocal, depth, result) );

   return SCIP_OKAY;
}

/*
 * separator specific interface methods
 */

/** creates the path cut separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaPath(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   BMSclearMemory(sepadata);

   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpPath, NULL, sepadata) );
   assert(sepa != NULL);

   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyPath) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreePath) );
   SCIP_CALL( SCIPsetSepaInitsol(scip, sepa, sepaInitsolPath) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolPath) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxrounds",
         "maximal number of path separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxroundsroot",
         "maximal number of path separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxsepacuts",
         "maximal number of path cuts separated per separation round",
         &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxsepacutsroot",
         "maximal number of path cuts separated per separation round in root node",
         &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxcutcands",
         "maximal number of path cut candidates collected per separation round",
         &sepadata->maxcutcands, FALSE, DEFAULT_MAXCUTCANDS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxtries",
         "maximal number of start rows tried per separation round in tree",
         &sepadata->maxtries, TRUE, DEFAULT_MAXTRIES, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxtriesroot",
         "maximal number of start rows tried per separation round in root node",
         &sepadata->maxtriesroot, TRUE, DEFAULT_MAXTRIESROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxpathlen",
         "maximal number of rows aggregated along a directed continuous path",
         &sepadata->maxpathlen, TRUE, DEFAULT_MAXPATHLEN, 2, MAXPATHLENHARD, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/" SEPA_NAME "/maxaggdensity",
         "maximal density of aggregated rows considered along a path",
         &sepadata->maxaggdensity, TRUE, DEFAULT_MAXAGGDENSITY, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/densityoffset",
         "additional nonzeros allowed on top of maxaggdensity",
         &sepadata->densityoffset, TRUE, DEFAULT_DENSITYOFFSET, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/usepathmix",
         "should VBD-aware multi-row path-mixing cuts be separated?",
         &sepadata->usepathmix, TRUE, DEFAULT_USEPATHMIX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/tryflowcover",
         "should flow cover cuts be tried along directed paths?",
         &sepadata->tryflowcover, TRUE, DEFAULT_TRYFLOWCOVER, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/initseed",
         "initial seed for random arc selection and cut selection",
         &sepadata->initseed, TRUE, DEFAULT_INITSEED, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
