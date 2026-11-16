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

/**@file   sepa_machinesched.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  separator for release-date cuts of single-machine scheduling substructures
 *
 * Detects LP rows with exactly three entries of the form M y + a s_i - a s_j >= d (or <=) with a binary variable y and
 * two variables s_i, s_j with finite lower bounds. Such a row implies s_i >= s_j + p for one value of y, where p > 0
 * is the larger of d and d - M (after scaling by 1/a). These implications form a directed graph on the job start times.
 * Greedily, starting from the job with most incoming arcs, a set N of jobs is built in which every pair of jobs is
 * ordered by one binary variable in both directions (a single machine). For each job j in N, the release-date cut
 *
 *    s_j >= min{ lb(s_i) : i in N } + sum_{i in N \\ j} p_i [i before j]
 *
 * is added to the global cut pool, where p_i is the minimal distance from job i to any other job in N that is ordered
 * after it, and [i before j] is the literal of the binary variable ordering i before j. The detection runs once per
 * solving run at the root node.
 *
 * Based on HiGHS (HighsMachineSchedSeparator, https://github.com/ERGO-Code/HiGHS/pull/2990) and
 *  G. L. Nemhauser and M. W. P. Savelsbergh, "A cutting plane algorithm for the single machine scheduling problem with
 *  release times", 1992; R. L. H. Wijfjes, "Recognition and Exploitation of Single-Machine Scheduling Subproblems in
 *  Mixed Integer Programs", MSc thesis, 2022.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_sepa.h"
#include "scip/pub_var.h"
#include "scip/scip_cut.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_pricer.h"
#include "scip/scip_prob.h"
#include "scip/scip_sepa.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include "scip/sepa_machinesched.h"
#include <string.h>

#define SEPA_NAME              "machinesched"
#define SEPA_DESC              "release-date cuts for single-machine scheduling substructures"
#define SEPA_PRIORITY             -1000
#define SEPA_FREQ                     0
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXMACHINES         100 /**< maximal number of machines (job sets) to detect */
#define DEFAULT_MINJOBS               3 /**< minimal number of jobs on a machine */
#define MAXORDCOLS               200000 /**< maximal number of columns for which order records are kept */
#define DEFAULT_RELEASECUTS        TRUE /**< should cuts for all release-date thresholds be generated? */
#define DEFAULT_ORDEREQUALITIES    TRUE /**< should binary variables determining the same order be identified? */

/** separator data */
struct SCIP_SepaData
{
   int                   maxmachines;        /**< maximal number of machines (job sets) to detect */
   int                   minjobs;            /**< minimal number of jobs on a machine */
   SCIP_Bool             orderequalities;    /**< should binary variables determining the same order be identified? */
   SCIP_Bool             releasecuts;        /**< should cuts for all release-date thresholds be generated? */
   SCIP_Bool             tried;              /**< was the detection already performed in this run? */
   int                   ncuts;              /**< number of cuts generated */
   int                   nmachines;          /**< number of machines detected */
};

/** directed arcs between jobs: arcto starts at least arcp after arcfrom if arcbin takes value arcval */
typedef struct Arcs
{
   int*                  from;               /**< LP position of the job that comes first */
   int*                  to;                 /**< LP position of the job that comes second */
   int*                  bin;                /**< LP position of the binary variable */
   SCIP_Real*            p;                  /**< minimal distance of the start times */
   SCIP_Bool*            val;                /**< value of the binary variable implying the order */
   int                   narcs;              /**< number of arcs */
   int                   size;               /**< size of the arrays */
   SCIP_HASHMAP*         map;                /**< map from (from, to) to the arc index */
   int                   ncols;              /**< number of LP columns (for the map key) */
   SCIP_Longint*         ordkeys;            /**< order records with u < v: ((u * ncols + v) * ncols + bin) * 16 + bits,
                                              *   bits: 1: u -> v if bin = 0, 2: u -> v if bin = 1, 4: v -> u if bin = 0,
                                              *   8: v -> u if bin = 1 */
   int                   nords;              /**< number of order records */
   int                   ordsize;            /**< size of the order record arrays */
} ARCS;

/*
 * Local methods
 */

/** returns the map key of an arc */
static
void* arcKey(
   ARCS*                 arcs,               /**< arcs */
   int                   from,               /**< first job */
   int                   to                  /**< second job */
   )
{
   return (void*)((size_t)from * (size_t)arcs->ncols + (size_t)to + 1); /*lint !e571*/
}

/** returns the index of the arc from -> to, or -1 */
static
int findArc(
   ARCS*                 arcs,               /**< arcs */
   int                   from,               /**< first job */
   int                   to                  /**< second job */
   )
{
   void* key = arcKey(arcs, from, to);

   if( !SCIPhashmapExists(arcs->map, key) )
      return -1;

   return SCIPhashmapGetImageInt(arcs->map, key);
}

/** adds the arc from -> to or strengthens an existing one */
static
SCIP_RETCODE addArc(
   SCIP*                 scip,               /**< SCIP data structure */
   ARCS*                 arcs,               /**< arcs */
   int                   from,               /**< first job */
   int                   to,                 /**< second job */
   int                   bin,                /**< binary variable */
   SCIP_Real             p,                  /**< distance */
   SCIP_Bool             val                 /**< value of the binary variable implying the order */
   )
{
   int idx;

   idx = findArc(arcs, from, to);

   if( idx >= 0 )
   {
      if( arcs->p[idx] < p )
      {
         arcs->p[idx] = p;
         arcs->bin[idx] = bin;
         arcs->val[idx] = val;
      }
      return SCIP_OKAY;
   }

   if( arcs->narcs == arcs->size )
   {
      int newsize = SCIPcalcMemGrowSize(scip, arcs->narcs + 1);

      SCIP_CALL( SCIPreallocBufferArray(scip, &arcs->from, newsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &arcs->to, newsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &arcs->bin, newsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &arcs->p, newsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &arcs->val, newsize) );
      arcs->size = newsize;
   }

   idx = arcs->narcs++;
   arcs->from[idx] = from;
   arcs->to[idx] = to;
   arcs->bin[idx] = bin;
   arcs->p[idx] = p;
   arcs->val[idx] = val;
   SCIP_CALL( SCIPhashmapInsertInt(arcs->map, arcKey(arcs, from, to), idx) );

   return SCIP_OKAY;
}

/** records the implication of a normalized row M y + s_i - s_j >= d, i.e., s_i >= s_j + d (y = 0), s_i >= s_j + d - M
 *  (y = 1), if one of the distances is positive
 */
static
SCIP_RETCODE tryAddArc(
   SCIP*                 scip,               /**< SCIP data structure */
   ARCS*                 arcs,               /**< arcs */
   int                   posjob,             /**< job i (coefficient +1) */
   int                   negjob,             /**< job j (coefficient -1) */
   int                   bin,                /**< binary variable y */
   SCIP_Real             d0,                 /**< distance if y = 0 */
   SCIP_Real             d1                  /**< distance if y = 1 */
   )
{
   if( !SCIPisPositive(scip, d0) && !SCIPisPositive(scip, d1) )
      return SCIP_OKAY;

   if( d0 > d1 )
   {
      SCIP_CALL( addArc(scip, arcs, negjob, posjob, bin, d0, FALSE) );
   }
   else
   {
      SCIP_CALL( addArc(scip, arcs, negjob, posjob, bin, d1, TRUE) );
   }

   /* record which value of the binary variable implies which order of the two jobs (if the key fits) */
   if( arcs->ncols > MAXORDCOLS )
      return SCIP_OKAY;

   if( arcs->nords == arcs->ordsize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, arcs->nords + 1);

      SCIP_CALL( SCIPreallocBufferArray(scip, &arcs->ordkeys, newsize) );
      arcs->ordsize = newsize;
   }
   {
      int u = MIN(negjob, posjob);
      int v = MAX(negjob, posjob);
      SCIP_Bool val = (d0 <= d1);

      int bits;

      /* arc negjob -> posjob */
      if( negjob == u )
         bits = val ? 2 : 1;
      else
         bits = val ? 8 : 4;
      arcs->ordkeys[arcs->nords] = (((SCIP_Longint)u * arcs->ncols + v) * arcs->ncols + bin) * 16 + bits;
      ++arcs->nords;
   }

   return SCIP_OKAY;
}

/** collects the arcs from the LP rows */
static
SCIP_RETCODE collectArcs(
   SCIP*                 scip,               /**< SCIP data structure */
   ARCS*                 arcs                /**< arcs */
   )
{
   SCIP_ROW** rows;
   int nrows;
   int r;

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   for( r = 0; r < nrows; ++r )
   {
      SCIP_ROW* row = rows[r];
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real scale;
      SCIP_Real bincoef = 0.0;
      int bin = -1;
      int pos = -1;
      int neg = -1;
      int k;

      /* modifiable rows may receive further columns by pricing, which would invalidate the inferred orderings */
      if( SCIProwIsLocal(row) || SCIProwIsModifiable(row) || SCIProwGetNNonz(row) != 3 || SCIProwGetNLPNonz(row) != 3 )
         continue;

      lhs = SCIProwGetLhs(row);
      rhs = SCIProwGetRhs(row);
      if( SCIPisEQ(scip, lhs, rhs) )
         continue;

      rowcols = SCIProwGetCols(row);
      rowvals = SCIProwGetVals(row);

      for( k = 0; k < 3; ++k )
      {
         SCIP_VAR* var = SCIPcolGetVar(rowcols[k]);

         if( SCIPvarIsBinary(var) )
         {
            if( bin >= 0 )
               break;
            bin = SCIPcolGetLPPos(rowcols[k]);
            bincoef = rowvals[k];
         }
         else if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
            break;
         else if( rowvals[k] > 0.0 )
         {
            if( pos >= 0 )
               break;
            pos = k;
         }
         else
         {
            if( neg >= 0 )
               break;
            neg = k;
         }
      }

      if( k < 3 || bin < 0 || pos < 0 || neg < 0 )
         continue;

      /* the two start time variables need coefficients of exactly the same absolute value; with a tolerance, the
       * discarded residual could make the inferred distances invalid for start times of large magnitude
       */
      scale = rowvals[pos];
      if( scale != -rowvals[neg] ) /*lint !e777*/
         continue;

      bincoef /= scale;
      pos = SCIPcolGetLPPos(rowcols[pos]);
      neg = SCIPcolGetLPPos(rowcols[neg]);
      assert(pos >= 0 && neg >= 0 && bin >= 0);

      /* M y + s_pos - s_neg >= lhs */
      if( !SCIPisInfinity(scip, -lhs) )
      {
         SCIP_Real d = (lhs - SCIProwGetConstant(row)) / scale;

         SCIP_CALL( tryAddArc(scip, arcs, pos, neg, bin, d, d - bincoef) );
      }

      /* -M y + s_neg - s_pos >= -rhs */
      if( !SCIPisInfinity(scip, rhs) )
      {
         SCIP_Real d = -(rhs - SCIProwGetConstant(row)) / scale;

         SCIP_CALL( tryAddArc(scip, arcs, neg, pos, bin, d, d + bincoef) );
      }
   }

   return SCIP_OKAY;
}

/** adds cliques forcing binary variables that determine the order of the same pair of jobs to be equal (or
 *  complementary): if y = 0 implies u before v and y = 1 implies v before u, then y is determined by the order
 */
static
SCIP_RETCODE addOrderEqualities(
   SCIP*                 scip,               /**< SCIP data structure */
   ARCS*                 arcs,               /**< arcs */
   SCIP_COL**            cols,               /**< LP columns */
   SCIP_Bool*            cutoff,             /**< pointer to store whether a cutoff was detected */
   int*                  nbdchgs             /**< pointer to count bound changes */
   )
{
   SCIP_Longint pairkey = -1;
   int basebin = -1;
   SCIP_Bool baseval = FALSE;
   int i;

   if( arcs->nords < 2 )
      return SCIP_OKAY;

   SCIPsortLong(arcs->ordkeys, arcs->nords);

   i = 0;
   while( i < arcs->nords && !(*cutoff) )
   {
      SCIP_Longint key = arcs->ordkeys[i] / 16;
      int bits = 0;
      int bin;

      /* combine the bits of the records of the same (u, v, bin) */
      while( i < arcs->nords && arcs->ordkeys[i] / 16 == key )
         bits |= (int)(arcs->ordkeys[i++] % 16);

      bin = (int)(key % arcs->ncols);

      if( key / arcs->ncols != pairkey )
      {
         pairkey = key / arcs->ncols;
         basebin = -1;
      }

      /* bin = 0 => u -> v and bin = 1 => v -> u, or bin = 1 => u -> v and bin = 0 => v -> u */
      if( (bits & 9) == 9 || (bits & 6) == 6 )
      {
         SCIP_Bool val = ((bits & 6) == 6);

         if( basebin < 0 )
         {
            basebin = bin;
            baseval = val;
         }
         else if( bin != basebin )
         {
            SCIP_VAR* vars[2];
            SCIP_Bool vals[2];
            SCIP_Bool infeasible;
            int nchg;

            /* literal (basebin == baseval) is equivalent to literal (bin == val) */
            vars[0] = SCIPcolGetVar(cols[basebin]);
            vars[1] = SCIPcolGetVar(cols[bin]);
            vals[0] = baseval;
            vals[1] = !val;
            SCIP_CALL( SCIPaddClique(scip, vars, vals, 2, FALSE, &infeasible, &nchg) );
            *nbdchgs += nchg;
            if( infeasible )
            {
               *cutoff = TRUE;
               break;
            }

            vals[0] = !baseval;
            vals[1] = val;
            SCIP_CALL( SCIPaddClique(scip, vars, vals, 2, FALSE, &infeasible, &nchg) );
            *nbdchgs += nchg;
            if( infeasible )
               *cutoff = TRUE;
         }
      }
   }

   return SCIP_OKAY;
}

/** greedily detects machines and adds release-date cuts to the global cut pool */
static
SCIP_RETCODE detectMachines(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   ARCS*                 arcs,               /**< arcs */
   SCIP_COL**            cols,               /**< LP columns */
   int                   ncols,              /**< number of LP columns */
   SCIP_Bool*            cutoff,             /**< pointer to store whether a cutoff was detected */
   int*                  ncuts               /**< pointer to count the number of cuts */
   )
{
   SCIP_Bool* used;
   SCIP_Real* proc;
   SCIP_Real* thresholds;
   SCIP_Real* keys;
   int* degree;
   int* cands;
   int* jobs;
   int m;

   SCIP_CALL( SCIPallocClearBufferArray(scip, &used, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &degree, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cands, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &keys, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &jobs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &proc, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &thresholds, ncols) );

   for( m = 0; m < sepadata->maxmachines && !(*cutoff); ++m )
   {
      SCIP_Real release;
      int largest;
      int ncands;
      int njobs;
      int a;
      int i;

      /* compute in-degrees among unused jobs */
      BMSclearMemoryArray(degree, ncols);
      for( a = 0; a < arcs->narcs; ++a )
      {
         if( !used[arcs->from[a]] && !used[arcs->to[a]] )
            ++degree[arcs->to[a]];
      }

      largest = -1;
      for( i = 0; i < ncols; ++i )
      {
         if( !used[i] && degree[i] > 0 && (largest < 0 || degree[i] > degree[largest]) )
            largest = i;
      }

      if( largest < 0 || degree[largest] < sepadata->minjobs - 1 )
         break;

      /* potential jobs: predecessors of the largest-degree job, sorted by decreasing degree */
      ncands = 0;
      for( a = 0; a < arcs->narcs; ++a )
      {
         if( arcs->to[a] == largest && !used[arcs->from[a]] )
         {
            cands[ncands] = arcs->from[a];
            keys[ncands] = (SCIP_Real)degree[arcs->from[a]] - 1e-7 * arcs->from[a];
            ++ncands;
         }
      }
      SCIPsortDownRealInt(keys, cands, ncands);

      jobs[0] = largest;
      proc[0] = SCIPinfinity(scip);
      njobs = 1;
      release = SCIPvarGetLbGlobal(SCIPcolGetVar(cols[largest]));

      for( i = 0; i < ncands; ++i )
      {
         SCIP_Real minp;
         SCIP_Bool valid = TRUE;
         int col = cands[i];
         int j;

         for( j = 0; j < njobs && valid; ++j )
         {
            int fromarc = findArc(arcs, col, jobs[j]);
            int toarc = findArc(arcs, jobs[j], col);

            valid = (fromarc >= 0 && toarc >= 0 && arcs->bin[fromarc] == arcs->bin[toarc]
               && arcs->val[fromarc] != arcs->val[toarc]);
         }

         if( !valid )
            continue;

         minp = SCIPinfinity(scip);
         for( j = 0; j < njobs; ++j )
         {
            int fromarc = findArc(arcs, col, jobs[j]);
            int toarc = findArc(arcs, jobs[j], col);

            minp = MIN(minp, arcs->p[fromarc]);
            proc[j] = MIN(proc[j], arcs->p[toarc]);
         }

         release = MIN(release, SCIPvarGetLbGlobal(SCIPcolGetVar(cols[col])));
         proc[njobs] = minp;
         jobs[njobs++] = col;
      }

      if( njobs < sepadata->minjobs )
         break;

      ++sepadata->nmachines;
      SCIPdebugMsg(scip, "machine %d with %d jobs, release date %g\n", m, njobs, release);

      /* release-date cuts for each job i and each threshold t <= r_i among the release dates of the machine:
       *    s_i >= t + sum_{j != i, r_j >= t} p_j [j before i],
       * i.e., -s_i + sum_{j != i, r_j >= t} p_j [j before i] <= -t; the smallest threshold yields the cut over all jobs
       */
      for( i = 0; i < njobs; ++i )
         thresholds[i] = SCIPvarGetLbGlobal(SCIPcolGetVar(cols[jobs[i]]));
      SCIPsortReal(thresholds, njobs);

      for( i = 0; i < njobs && !(*cutoff); ++i )
      {
         SCIP_Real ri = SCIPvarGetLbGlobal(SCIPcolGetVar(cols[jobs[i]]));
         int t;

         for( t = 0; t < njobs && !(*cutoff); ++t )
         {
            SCIP_ROW* cut;
            char name[SCIP_MAXSTRLEN];
            SCIP_Real thres;
            SCIP_Real rhs;
            int nterms;
            int j;

            /* skip duplicate thresholds and, unless enabled, all thresholds but the smallest */
            if( t > 0 && (thresholds[t] == thresholds[t - 1] || !sepadata->releasecuts) ) /*lint !e777*/
               continue;

            thres = (t == 0) ? release : thresholds[t];
            if( thres > ri )
               break;

            rhs = -thres;
            nterms = 0;

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "machinesched_%d_%d_%d", m, i, t);
            SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, name, -SCIPinfinity(scip), SCIPinfinity(scip),
                  FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

            for( j = 0; j < njobs; ++j )
            {
               int arc;

               if( j == i || SCIPvarGetLbGlobal(SCIPcolGetVar(cols[jobs[j]])) < thres )
                  continue;

               arc = findArc(arcs, jobs[j], jobs[i]);
               assert(arc >= 0);
               assert(!SCIPisInfinity(scip, proc[j]));

               if( arcs->val[arc] )
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, SCIPcolGetVar(cols[arcs->bin[arc]]), proc[j]) );
               }
               else
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, SCIPcolGetVar(cols[arcs->bin[arc]]), -proc[j]) );
                  rhs -= proc[j];
               }
               ++nterms;
            }

            if( nterms > 0 )
            {
               SCIP_CALL( SCIPaddVarToRow(scip, cut, SCIPcolGetVar(cols[jobs[i]]), -1.0) );
               SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
               SCIP_CALL( SCIPchgRowRhs(scip, cut, rhs) );

#ifdef SCIP_DEBUG
               SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif
               SCIP_CALL( SCIPaddPoolCut(scip, cut) );

               if( SCIPisCutEfficacious(scip, NULL, cut) )
               {
                  SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
               }
               ++(*ncuts);
            }

            SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         }
      }

      for( i = 0; i < njobs; ++i )
         used[jobs[i]] = TRUE;
   }

   SCIPfreeBufferArray(scip, &thresholds);
   SCIPfreeBufferArray(scip, &proc);
   SCIPfreeBufferArray(scip, &jobs);
   SCIPfreeBufferArray(scip, &keys);
   SCIPfreeBufferArray(scip, &cands);
   SCIPfreeBufferArray(scip, &degree);
   SCIPfreeBufferArray(scip, &used);

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyMachinesched)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);

   SCIP_CALL( SCIPincludeSepaMachinesched(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeMachinesched)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeBlockMemory(scip, &sepadata);
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of separator (called when branch and bound process is about to begin) */
static
SCIP_DECL_SEPAINITSOL(sepaInitsolMachinesched)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   sepadata->tried = FALSE;

   return SCIP_OKAY;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpMachinesched)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   SCIP_COL** cols;
   SCIP_Bool cutoff;
   ARCS arcs;
   int nbdchgs;
   int ncols;
   int ncuts;

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   if( sepadata->tried || SCIPgetDepth(scip) > 0 || sepadata->maxmachines <= 0 )
      return SCIP_OKAY;

   /* with active pricers, rows can change, so the detected structure may become invalid */
   if( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   sepadata->tried = TRUE;

   if( SCIPgetNBinVars(scip) == 0 || SCIPgetNBinVars(scip) == SCIPgetNVars(scip) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   if( ncols == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   BMSclearMemory(&arcs);
   arcs.ncols = ncols;
   SCIP_CALL( SCIPhashmapCreate(&arcs.map, SCIPblkmem(scip), 64) );

   SCIP_CALL( collectArcs(scip, &arcs) );

   cutoff = FALSE;
   ncuts = 0;
   nbdchgs = 0;
   if( sepadata->orderequalities )
   {
      SCIP_CALL( addOrderEqualities(scip, &arcs, cols, &cutoff, &nbdchgs) );
   }
   if( !cutoff && arcs.narcs >= 3 )
   {
      SCIP_CALL( detectMachines(scip, sepa, sepadata, &arcs, cols, ncols, &cutoff, &ncuts) );
   }

   SCIPhashmapFree(&arcs.map);
   SCIPfreeBufferArrayNull(scip, &arcs.ordkeys);
   SCIPfreeBufferArrayNull(scip, &arcs.val);
   SCIPfreeBufferArrayNull(scip, &arcs.p);
   SCIPfreeBufferArrayNull(scip, &arcs.bin);
   SCIPfreeBufferArrayNull(scip, &arcs.to);
   SCIPfreeBufferArrayNull(scip, &arcs.from);

   sepadata->ncuts += ncuts;

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nbdchgs > 0 )
      *result = SCIP_REDUCEDDOM;
   else if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}

/*
 * separator specific interface methods
 */

/** creates the machinesched separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaMachinesched(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   BMSclearMemory(sepadata);

   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpMachinesched, NULL, sepadata) );
   assert(sepa != NULL);

   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyMachinesched) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeMachinesched) );
   SCIP_CALL( SCIPsetSepaInitsol(scip, sepa, sepaInitsolMachinesched) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxmachines",
         "maximal number of machines (sets of pairwise ordered jobs) to detect",
         &sepadata->maxmachines, FALSE, DEFAULT_MAXMACHINES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/minjobs",
         "minimal number of jobs on a machine",
         &sepadata->minjobs, TRUE, DEFAULT_MINJOBS, 2, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/orderequalities",
         "should binary variables that determine the order of the same pair of jobs be identified via cliques?",
         &sepadata->orderequalities, TRUE, DEFAULT_ORDEREQUALITIES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/releasecuts",
         "should cuts be generated for all release-date thresholds (otherwise only for the smallest release date)?",
         &sepadata->releasecuts, TRUE, DEFAULT_RELEASECUTS, NULL, NULL) );

   return SCIP_OKAY;
}
