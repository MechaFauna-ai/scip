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

/**@file   sepa_modk.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  mod-k cuts separator for small primes k
 *
 * Separates mod-k Chvátal-Gomory cuts for the primes k = 3, 5, 7 (and optionally k = 2, which is usually covered by
 * the zerohalf separator) following the approach of the mod-k separator of HiGHS:
 *
 * -# Every LP row that is tight at the current LP solution and does not contain a continuous variable strictly between
 *    its bounds is written in "distance from bound" space: each column is complemented to its closest bound, so that
 *    all transformed variables are nonnegative. Only integral columns with a positive transformed value are kept, the
 *    row is scaled to integers and its (transformed, scaled) right hand side is appended as an extra column.
 * -# For each k, the linear system \f$ w^T \bar A \equiv 0 \pmod k \f$ on all kept columns and \f$ w^T \bar b \equiv 1
 *    \pmod k \f$ is solved over GF(k) by sparse Gaussian elimination. A bounded number of basic solutions \f$ w \f$ is
 *    enumerated by swapping nonbasic unknowns into the triangular factor one at a time.
 * -# For each solution, the original LP rows are aggregated with weights \f$ ((k-1) w_i \bmod k) / k \f$ (giving a
 *    right hand side with fractional part \f$ (k-1)/k \f$, i.e., a maximally violated mod-k cut) and, for k > 2, also
 *    with weights \f$ w_i / k \f$. From each aggregation, the MIR cut with scale 1 (the strengthened mod-k cut) is
 *    computed; optionally, lifted knapsack cover and flow cover cuts and the cMIR scaling heuristic are tried as well,
 *    and the most efficacious cut is kept.
 * -# The collected cuts are filtered with the hybrid cut selector; global cuts go to the cut pool.
 *
 * The validity of the generated cuts does not depend on the mod-k computations: every cut is computed by SCIP's cut
 * generation routines (SCIPcalcMIR(), SCIPcutGenerationHeuristicCMIR(), SCIPcalcKnapsackCover(), SCIPcalcFlowCover())
 * from a combination of LP rows in which each row is used with its valid side. The mod-k system only selects the
 * aggregation weights.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cuts.h"
#include "scip/cutsel_hybrid.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
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
#include "scip/sepa_modk.h"
#include <stdint.h>
#include <string.h>

#define SEPA_NAME              "modk"
#define SEPA_DESC              "mod-k cuts separator for small primes k"
#define SEPA_PRIORITY             -6100
#define SEPA_FREQ                     0
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE
#define SEPA_DELAY                FALSE

#define DEFAULT_MAXROUNDS             5 /**< maximal number of mod-k separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        -1 /**< maximal number of mod-k separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          50 /**< maximal number of mod-k cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT     200 /**< maximal number of mod-k cuts separated per separation round in root node */
#define DEFAULT_MAXCUTCANDS        2000 /**< maximal number of mod-k cut candidates collected per separation round */
#define DEFAULT_MAXSOLS             500 /**< maximal number of GF(k) solutions enumerated per k and separation round */
#define DEFAULT_USEMOD2           FALSE /**< should also mod-2 cuts be separated (usually done by zerohalf)? */
#define DEFAULT_STOPATFIRSTK      FALSE /**< should the separation stop at the first k that produced cuts? */
#define DEFAULT_CMIRHEUR           TRUE /**< should the cMIR scaling heuristic be tried on the aggregations? */
#define DEFAULT_MAXCMIRHEUR          20 /**< maximal number of aggregations per call on which the cMIR heuristic is tried */
#define DEFAULT_TRYCOVERS          TRUE /**< should lifted knapsack and flow cover cuts be tried on the aggregations? */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_INITSEED         0x40D7 /**< initial seed used for random tie-breaking in cut selection */

#define BOUNDSWITCH                 0.5 /**< threshold for bound switching in the MIR transformation */
#define POSTPROCESS                TRUE /**< apply postprocessing after MIR calculation */
#define USEVBDS                       0 /**< variable types for which variable bound substitution is allowed (none) */
#define FIXINTEGRALRHS            FALSE /**< should complementation be adjusted to make the rhs fractional? */
#define MINFRAC                    0.05 /**< minimal fractionality of rhs to produce MIR cut for */
#define MAXFRAC                   0.999 /**< maximal fractionality of rhs to produce MIR cut for */
#define MAXTESTDELTA                 10 /**< maximal number of deltas tested by the cMIR heuristic */
#define MAXINTSCALE               1e+06 /**< maximal scalar to make a row integral */
#define MAXDNOM                 1000000LL /**< maximal denominator for computing the integral scalar */
#define MAXFILLFACTOR                10 /**< maximal growth of the number of nonzeros during GF(k) elimination */
#define MAXFILLOFFSET            100000 /**< additional number of nonzeros allowed during GF(k) elimination */
#define MAXWORK            50000000LL /**< maximal number of elementary operations in one GF(k) elimination */
#define MAXAGGRLEN(ncols)        ((int)(0.1*(ncols)+1000)) /**< maximal number of integral entries of a used row */

/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator for tie-breaking in cut selection */
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int                   maxsepacutsroot;    /**< maximal number of cuts separated per separation round in root node */
   int                   maxcutcands;        /**< maximal number of cut candidates collected per separation round */
   int                   maxsols;            /**< maximal number of GF(k) solutions enumerated per k and round */
   int                   initseed;           /**< initial seed for the random number generator */
   SCIP_Bool             usemod2;            /**< should also mod-2 cuts be separated? */
   SCIP_Bool             stopatfirstk;       /**< should the separation stop at the first k that produced cuts? */
   SCIP_Bool             cmirheur;           /**< should the cMIR scaling heuristic be tried on the aggregations? */
   int                   maxcmirheur;        /**< maximal number of aggregations per call on which the cMIR heuristic is tried */
   SCIP_Bool             trycovers;          /**< should lifted knapsack and flow cover cuts be tried on the aggregations? */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
};

/** integral row system extracted from the tight LP rows */
typedef struct ModkRowSystem
{
   SCIP_ROW**            rows;               /**< LP rows used in the system */
   int*                  sign;               /**< +1 if the right hand side is used, -1 if the left hand side is used */
   SCIP_Real*            scale;              /**< integral scalar of each row */
   SCIP_Longint*         rhs;                /**< scaled, transformed integral right hand side of each row */
   int*                  beg;                /**< start of each row in the entry arrays (size nrows+1) */
   int*                  eqind;              /**< equation (compressed column) index of each entry */
   SCIP_Longint*         val;                /**< scaled, transformed integral coefficient of each entry */
   int                   nrows;              /**< number of rows in the system */
   int                   nentries;           /**< number of entries */
   int                   neqs;               /**< number of distinct columns in the system (without rhs equation) */
   int                   rowssize;           /**< allocated size of the row arrays */
   int                   entriessize;        /**< allocated size of the entry arrays */
} MODKROWSYSTEM;

/** linear system over GF(k): equations (one per column plus one for the right hand side) over the unknown row weights */
typedef struct GFkSystem
{
   int                   k;                  /**< the prime modulus */
   int                   neqs;               /**< number of equations */
   int                   nunks;              /**< number of unknowns */
   int**                 eqidx;              /**< sorted unknown indices of each equation */
   unsigned char**       eqval;              /**< values (1..k-1) of each equation */
   int*                  eqlen;              /**< length of each equation */
   int*                  eqsize;             /**< allocated size of each equation */
   unsigned char*        rhs;                /**< right hand side of each equation */
   SCIP_Bool*            eqused;             /**< has the equation been used as pivot row? */
   int**                 occ;                /**< equations containing each unknown (may contain stale entries) */
   int*                  occlen;             /**< length of occurrence lists */
   int*                  occsize;            /**< allocated size of occurrence lists */
   int*                  colcount;           /**< number of unused equations containing each unknown */
   signed char*          status;             /**< 1: basic, 0: nonbasic, -1: nonbasic scheduled for a basis swap */
   int*                  pivrow;             /**< pivot equation of each pivot step */
   int*                  pivcol;             /**< pivot unknown of each pivot step */
   int                   npiv;               /**< number of pivot steps */
   int*                  heapkey;            /**< keys (column counts) of the priority queue */
   int*                  heapcol;            /**< unknowns of the priority queue */
   int                   heaplen;            /**< number of elements in the priority queue */
   int                   heapsize;           /**< allocated size of the priority queue */
   int*                  tmpidx;             /**< temporary buffer for merging equations */
   unsigned char*        tmpval;             /**< temporary buffer for merging equations */
   int                   tmpsize;            /**< size of temporary buffers */
   SCIP_Longint          nnz;                /**< current number of nonzeros */
   SCIP_Longint          maxnnz;             /**< maximal number of nonzeros allowed */
   SCIP_Longint          work;               /**< performed work */
} GFKSYSTEM;

/** context for turning GF(k) solutions into cuts */
typedef struct ModkCutContext
{
   SCIP_SEPA*            sepa;               /**< the separator */
   SCIP_SEPADATA*        sepadata;           /**< separator data */
   MODKROWSYSTEM*        sys;                /**< the integral row system */
   SCIP_AGGRROW*         aggrrow;            /**< aggregation row */
   SCIP_HASHSET*         usedweights;        /**< hashes of weight vectors already processed for the current k */
   SCIP_Real*            cutcoefs;           /**< buffer for cut coefficients */
   int*                  cutinds;            /**< buffer for cut indices */
   SCIP_ROW**            cuts;               /**< collected cut candidates */
   int                   ncuts;              /**< number of collected cut candidates */
   int                   maxcuts;            /**< maximal number of cut candidates */
   SCIP_Bool             allowlocal;         /**< are local cuts allowed? */
   int                   k;                  /**< current modulus */
   int                   nsolsk;             /**< number of solutions processed for the current k */
   int                   nmirtried;          /**< number of MIR computations (statistics) */
   int                   nmirsuccess;        /**< number of successful MIR computations (statistics) */
   int                   nmireff;            /**< number of efficacious cuts (statistics) */
   int                   ncoverwins;         /**< number of times a cover cut was better (statistics) */
   int                   ncmirheur;          /**< number of calls of the cMIR heuristic in this separation call */
} MODKCUTCONTEXT;

/*
 * GF(k) arithmetic
 */

/** returns the multiplicative inverse of a in GF(k) for a prime k */
static
unsigned int gfkInverse(
   unsigned int          a,                  /**< nonzero element of GF(k) */
   int                   k                   /**< prime modulus */
   )
{
   unsigned int inv = 1;
   int e;

   assert(a > 0 && (int)a < k);

   /* Fermat: a^(k-2) is the inverse of a */
   for( e = 0; e < k - 2; ++e )
      inv = (inv * a) % (unsigned int)k;

   assert((inv * a) % (unsigned int)k == 1);

   return inv;
}

/** reduces an integer modulo k to a value in [0,k) */
static
unsigned char modK(
   SCIP_Longint          val,                /**< value */
   int                   k                   /**< modulus */
   )
{
   SCIP_Longint r = val % k;

   if( r < 0 )
      r += k;

   return (unsigned char) r;
}

/*
 * GF(k) system
 */

/** pushes an element to the priority queue */
static
SCIP_RETCODE heapPush(
   SCIP*                 scip,               /**< SCIP data structure */
   GFKSYSTEM*            gf,                 /**< GF(k) system */
   int                   key,                /**< key */
   int                   col                 /**< unknown */
   )
{
   int pos;

   if( gf->heaplen == gf->heapsize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, gf->heaplen + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &gf->heapkey, gf->heapsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &gf->heapcol, gf->heapsize, newsize) );
      gf->heapsize = newsize;
   }

   pos = gf->heaplen++;
   while( pos > 0 )
   {
      int parent = (pos - 1) / 2;
      if( gf->heapkey[parent] < key || (gf->heapkey[parent] == key && gf->heapcol[parent] <= col) )
         break;
      gf->heapkey[pos] = gf->heapkey[parent];
      gf->heapcol[pos] = gf->heapcol[parent];
      pos = parent;
   }
   gf->heapkey[pos] = key;
   gf->heapcol[pos] = col;

   return SCIP_OKAY;
}

/** pops the minimal element from the priority queue */
static
void heapPop(
   GFKSYSTEM*            gf,                 /**< GF(k) system */
   int*                  key,                /**< pointer to store the key */
   int*                  col                 /**< pointer to store the unknown */
   )
{
   int lastkey;
   int lastcol;
   int pos;

   assert(gf->heaplen > 0);

   *key = gf->heapkey[0];
   *col = gf->heapcol[0];

   --gf->heaplen;
   if( gf->heaplen == 0 )
      return;

   lastkey = gf->heapkey[gf->heaplen];
   lastcol = gf->heapcol[gf->heaplen];
   pos = 0;
   for( ;; )
   {
      int child = 2 * pos + 1;
      if( child >= gf->heaplen )
         break;
      if( child + 1 < gf->heaplen && (gf->heapkey[child + 1] < gf->heapkey[child]
            || (gf->heapkey[child + 1] == gf->heapkey[child] && gf->heapcol[child + 1] < gf->heapcol[child])) )
         ++child;
      if( lastkey < gf->heapkey[child] || (lastkey == gf->heapkey[child] && lastcol <= gf->heapcol[child]) )
         break;
      gf->heapkey[pos] = gf->heapkey[child];
      gf->heapcol[pos] = gf->heapcol[child];
      pos = child;
   }
   gf->heapkey[pos] = lastkey;
   gf->heapcol[pos] = lastcol;
}

/** appends an equation to the occurrence list of an unknown */
static
SCIP_RETCODE occAppend(
   SCIP*                 scip,               /**< SCIP data structure */
   GFKSYSTEM*            gf,                 /**< GF(k) system */
   int                   col,                /**< unknown */
   int                   eq                  /**< equation */
   )
{
   if( gf->occlen[col] == gf->occsize[col] )
   {
      int newsize = SCIPcalcMemGrowSize(scip, gf->occlen[col] + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &gf->occ[col], gf->occsize[col], newsize) );
      gf->occsize[col] = newsize;
   }
   gf->occ[col][gf->occlen[col]++] = eq;

   return SCIP_OKAY;
}

/** returns the position of an unknown in an equation, or -1 */
static
int eqFind(
   GFKSYSTEM*            gf,                 /**< GF(k) system */
   int                   eq,                 /**< equation */
   int                   col                 /**< unknown */
   )
{
   int pos;

   if( SCIPsortedvecFindInt(gf->eqidx[eq], col, gf->eqlen[eq], &pos) )
      return pos;

   return -1;
}

/** creates the GF(k) system from the integral row system; equation j < neqs corresponds to column j, equation neqs to
 *  the right hand side
 */
static
SCIP_RETCODE gfkCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   GFKSYSTEM*            gf,                 /**< GF(k) system to fill */
   MODKROWSYSTEM*        sys,                /**< integral row system */
   int                   k                   /**< prime modulus */
   )
{
   int* cnt;
   int i;
   int j;

   BMSclearMemory(gf);
   gf->k = k;
   gf->neqs = sys->neqs + 1;
   gf->nunks = sys->nrows;

   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &gf->eqidx, gf->neqs) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &gf->eqval, gf->neqs) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &gf->eqlen, gf->neqs) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &gf->eqsize, gf->neqs) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &gf->rhs, gf->neqs) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &gf->eqused, gf->neqs) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &gf->occ, gf->nunks) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &gf->occlen, gf->nunks) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &gf->occsize, gf->nunks) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &gf->colcount, gf->nunks) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &gf->status, gf->nunks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &gf->pivrow, MIN(gf->neqs, gf->nunks)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &gf->pivcol, MIN(gf->neqs, gf->nunks)) );

   /* count the nonzero entries per equation */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &cnt, gf->neqs) );
   for( i = 0; i < sys->nrows; ++i )
   {
      for( j = sys->beg[i]; j < sys->beg[i + 1]; ++j )
      {
         if( modK(sys->val[j], k) != 0 )
            ++cnt[sys->eqind[j]];
      }
      if( modK(sys->rhs[i], k) != 0 )
         ++cnt[sys->neqs];
   }

   for( j = 0; j < gf->neqs; ++j )
   {
      if( cnt[j] > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &gf->eqidx[j], cnt[j]) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &gf->eqval[j], cnt[j]) );
         gf->eqsize[j] = cnt[j];
      }
   }
   SCIPfreeBufferArray(scip, &cnt);

   /* fill equations; since the unknowns are processed in increasing order, the equations are sorted */
   for( i = 0; i < sys->nrows; ++i )
   {
      unsigned char v;

      for( j = sys->beg[i]; j < sys->beg[i + 1]; ++j )
      {
         int eq = sys->eqind[j];

         v = modK(sys->val[j], k);
         if( v == 0 )
            continue;

         assert(gf->eqlen[eq] == 0 || gf->eqidx[eq][gf->eqlen[eq] - 1] < i);
         gf->eqidx[eq][gf->eqlen[eq]] = i;
         gf->eqval[eq][gf->eqlen[eq]] = v;
         ++gf->eqlen[eq];
      }

      v = modK(sys->rhs[i], k);
      if( v != 0 )
      {
         gf->eqidx[sys->neqs][gf->eqlen[sys->neqs]] = i;
         gf->eqval[sys->neqs][gf->eqlen[sys->neqs]] = v;
         ++gf->eqlen[sys->neqs];
      }
   }

   /* the right hand side equation has right hand side 1, all others 0 */
   gf->rhs[sys->neqs] = 1;

   /* occurrence lists and column counts */
   for( j = 0; j < gf->neqs; ++j )
   {
      gf->nnz += gf->eqlen[j];
      for( i = 0; i < gf->eqlen[j]; ++i )
      {
         int col = gf->eqidx[j][i];
         SCIP_CALL( occAppend(scip, gf, col, j) );
         ++gf->colcount[col];
      }
   }

   gf->maxnnz = MAXFILLFACTOR * gf->nnz + MAXFILLOFFSET;
   gf->tmpsize = 0;
   gf->tmpidx = NULL;
   gf->tmpval = NULL;

   return SCIP_OKAY;
}

/** frees the GF(k) system */
static
void gfkFree(
   SCIP*                 scip,               /**< SCIP data structure */
   GFKSYSTEM*            gf                  /**< GF(k) system */
   )
{
   int j;

   SCIPfreeBlockMemoryArrayNull(scip, &gf->tmpval, gf->tmpsize);
   SCIPfreeBlockMemoryArrayNull(scip, &gf->tmpidx, gf->tmpsize);
   SCIPfreeBlockMemoryArrayNull(scip, &gf->heapcol, gf->heapsize);
   SCIPfreeBlockMemoryArrayNull(scip, &gf->heapkey, gf->heapsize);
   for( j = gf->nunks - 1; j >= 0; --j )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &gf->occ[j], gf->occsize[j]);
   }
   for( j = gf->neqs - 1; j >= 0; --j )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &gf->eqval[j], gf->eqsize[j]);
      SCIPfreeBlockMemoryArrayNull(scip, &gf->eqidx[j], gf->eqsize[j]);
   }
   SCIPfreeBlockMemoryArray(scip, &gf->pivcol, MIN(gf->neqs, gf->nunks));
   SCIPfreeBlockMemoryArray(scip, &gf->pivrow, MIN(gf->neqs, gf->nunks));
   SCIPfreeBlockMemoryArray(scip, &gf->status, gf->nunks);
   SCIPfreeBlockMemoryArray(scip, &gf->colcount, gf->nunks);
   SCIPfreeBlockMemoryArray(scip, &gf->occsize, gf->nunks);
   SCIPfreeBlockMemoryArray(scip, &gf->occlen, gf->nunks);
   SCIPfreeBlockMemoryArray(scip, &gf->occ, gf->nunks);
   SCIPfreeBlockMemoryArray(scip, &gf->eqused, gf->neqs);
   SCIPfreeBlockMemoryArray(scip, &gf->rhs, gf->neqs);
   SCIPfreeBlockMemoryArray(scip, &gf->eqsize, gf->neqs);
   SCIPfreeBlockMemoryArray(scip, &gf->eqlen, gf->neqs);
   SCIPfreeBlockMemoryArray(scip, &gf->eqval, gf->neqs);
   SCIPfreeBlockMemoryArray(scip, &gf->eqidx, gf->neqs);
}

/** replaces equation r by r + scale * p (all arithmetic in GF(k)) and updates column counts and occurrence lists */
static
SCIP_RETCODE gfkAddScaledEq(
   SCIP*                 scip,               /**< SCIP data structure */
   GFKSYSTEM*            gf,                 /**< GF(k) system */
   int                   r,                  /**< equation to modify */
   int                   p,                  /**< pivot equation */
   unsigned int          scale               /**< scale in 1..k-1 */
   )
{
   unsigned int k = (unsigned int) gf->k;
   int lenr = gf->eqlen[r];
   int lenp = gf->eqlen[p];
   int* ridx = gf->eqidx[r];
   unsigned char* rval = gf->eqval[r];
   int* pidx = gf->eqidx[p];
   unsigned char* pval = gf->eqval[p];
   int needed = lenr + lenp;
   int a = 0;
   int b = 0;
   int n = 0;

   if( needed > gf->tmpsize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, needed);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &gf->tmpidx, gf->tmpsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &gf->tmpval, gf->tmpsize, newsize) );
      gf->tmpsize = newsize;
   }

   gf->work += needed;

   while( a < lenr || b < lenp )
   {
      if( b >= lenp || (a < lenr && ridx[a] < pidx[b]) )
      {
         gf->tmpidx[n] = ridx[a];
         gf->tmpval[n] = rval[a];
         ++n;
         ++a;
      }
      else if( a >= lenr || pidx[b] < ridx[a] )
      {
         /* fill-in */
         unsigned char v = (unsigned char) ((scale * pval[b]) % k);
         assert(v != 0);
         gf->tmpidx[n] = pidx[b];
         gf->tmpval[n] = v;
         ++n;
         ++gf->colcount[pidx[b]];
         SCIP_CALL( occAppend(scip, gf, pidx[b], r) );
         ++b;
      }
      else
      {
         unsigned char v = (unsigned char) ((rval[a] + scale * pval[b]) % k);
         if( v != 0 )
         {
            gf->tmpidx[n] = ridx[a];
            gf->tmpval[n] = v;
            ++n;
         }
         else
            --gf->colcount[ridx[a]];
         ++a;
         ++b;
      }
   }

   gf->nnz += n - lenr;

   if( n > gf->eqsize[r] )
   {
      int newsize = SCIPcalcMemGrowSize(scip, n);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &gf->eqidx[r], gf->eqsize[r], newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &gf->eqval[r], gf->eqsize[r], newsize) );
      gf->eqsize[r] = newsize;
   }
   if( n > 0 )
   {
      BMScopyMemoryArray(gf->eqidx[r], gf->tmpidx, n);
      BMScopyMemoryArray(gf->eqval[r], gf->tmpval, n);
   }
   gf->eqlen[r] = n;

   gf->rhs[r] = (unsigned char) ((gf->rhs[r] + scale * gf->rhs[p]) % k);

   return SCIP_OKAY;
}

/** performs the sparse Gaussian elimination; returns whether the system is solvable (and the limits were respected) */
static
SCIP_RETCODE gfkEliminate(
   SCIP*                 scip,               /**< SCIP data structure */
   GFKSYSTEM*            gf,                 /**< GF(k) system */
   SCIP_Bool*            solvable            /**< pointer to store whether a solution exists */
   )
{
   unsigned int k = (unsigned int) gf->k;
   int maxpiv = MIN(gf->neqs, gf->nunks);
   int j;

   *solvable = FALSE;

   for( j = 0; j < gf->nunks; ++j )
   {
      if( gf->colcount[j] > 0 )
      {
         SCIP_CALL( heapPush(scip, gf, gf->colcount[j], j) );
      }
   }

   while( gf->heaplen > 0 && gf->npiv < maxpiv )
   {
      int key;
      int col;
      int pivot = -1;
      int pivotlen = INT_MAX;
      int pivotpos = -1;
      unsigned int pivinv;
      int i;
      int noccs;

      heapPop(gf, &key, &col);

      if( gf->colcount[col] == 0 )
         continue;

      assert(gf->status[col] == 0);

      if( gf->colcount[col] != key )
      {
         SCIP_CALL( heapPush(scip, gf, gf->colcount[col], col) );
         continue;
      }

      /* choose the shortest unused equation containing the column as pivot row */
      for( i = 0; i < gf->occlen[col]; ++i )
      {
         int eq = gf->occ[col][i];
         int pos;

         if( gf->eqused[eq] || gf->eqlen[eq] >= pivotlen )
            continue;

         pos = eqFind(gf, eq, col);
         if( pos < 0 )
            continue;

         pivot = eq;
         pivotlen = gf->eqlen[eq];
         pivotpos = pos;
      }
      assert(pivot >= 0);

      pivinv = gfkInverse(gf->eqval[pivot][pivotpos], gf->k);

      /* eliminate the column from all other unused equations; the occurrence list of col does not grow here, since
       * col is canceled and the pivot equation only contains unknowns that are not yet basic
       */
      noccs = gf->occlen[col];
      for( i = 0; i < noccs; ++i )
      {
         int eq = gf->occ[col][i];
         int pos;
         unsigned int scale;

         if( eq == pivot || gf->eqused[eq] )
            continue;

         pos = eqFind(gf, eq, col);
         if( pos < 0 )
            continue;

         scale = (pivinv * (k - gf->eqval[eq][pos])) % k;
         SCIP_CALL( gfkAddScaledEq(scip, gf, eq, pivot, scale) );
         assert(eqFind(gf, eq, col) < 0);

         if( gf->nnz > gf->maxnnz || gf->work > MAXWORK )
            return SCIP_OKAY;
      }
      assert(noccs == gf->occlen[col]);

      /* the occurrence list of col is not needed anymore */
      gf->occlen[col] = 0;

      gf->pivrow[gf->npiv] = pivot;
      gf->pivcol[gf->npiv] = col;
      ++gf->npiv;
      gf->status[col] = 1;
      gf->eqused[pivot] = TRUE;

      /* update column counts of the unknowns in the pivot row */
      for( i = 0; i < gf->eqlen[pivot]; ++i )
      {
         int c = gf->eqidx[pivot][i];

         --gf->colcount[c];
         assert(gf->colcount[c] >= 0);
         assert(c != col || gf->colcount[c] == 0);

         if( gf->colcount[c] > 0 )
         {
            SCIP_CALL( heapPush(scip, gf, gf->colcount[c], c) );
         }
      }
   }

   /* the system is solvable iff all remaining (now empty) equations have zero right hand side */
   for( j = 0; j < gf->neqs; ++j )
   {
      if( gf->eqused[j] )
         continue;

      assert(gf->eqlen[j] == 0);
      if( gf->rhs[j] != 0 )
         return SCIP_OKAY;
   }

   *solvable = TRUE;

   return SCIP_OKAY;
}

/*
 * Cut generation
 */

/** computes a hash value of a sorted weight vector */
static
void* hashWeights(
   int*                  idx,                /**< sorted unknown indices */
   int*                  w,                  /**< weights */
   int                   n                   /**< length */
   )
{
   uint64_t h = 14695981039346656037ULL;
   int i;

   for( i = 0; i < n; ++i )
   {
      h ^= (uint64_t)(unsigned int)idx[i] * 8ULL + (uint64_t)(unsigned int)w[i];
      h *= 1099511628211ULL;
   }

   /* never return NULL */
   h |= 1ULL;

   return (void*)(size_t)h;
}

/** computes an MIR cut from the current aggregation row and stores it as a candidate if it is efficacious */
static
SCIP_RETCODE computeCut(
   SCIP*                 scip,               /**< SCIP data structure */
   MODKCUTCONTEXT*       ctx,                /**< cut context */
   int                   variant             /**< variant (only used for naming) */
   )
{
   SCIP_Real cutrhs;
   SCIP_Real cutefficacy;
   int cutnnz;
   int cutrank;
   SCIP_Bool cutislocal;
   SCIP_Bool success;
   SCIP_Bool valid;

   SCIPaggrRowRemoveZeros(scip, ctx->aggrrow, FALSE, &valid);
   if( !valid || SCIPaggrRowGetNNz(ctx->aggrrow) == 0 )
      return SCIP_OKAY;

   ++ctx->nmirtried;

   /* each method only returns a cut if it is efficacious and strictly better than the current best efficacy; on
    * failure, the cut arrays are left as they were
    */
   success = FALSE;
   cutnnz = 0;
   cutrhs = 0.0;
   cutrank = 0;
   cutislocal = FALSE;
   cutefficacy = -SCIPinfinity(scip);

   /* MIR cut with scale 1, i.e., the (strengthened) mod-k cut */
   {
      SCIP_Real eff = cutefficacy;
      SCIP_Real rhs;
      int nnz;
      int rank;
      SCIP_Bool islocal;
      SCIP_Bool found;

      SCIP_CALL( SCIPcalcMIR(scip, NULL, POSTPROCESS, BOUNDSWITCH, USEVBDS, ctx->allowlocal, FIXINTEGRALRHS, NULL, NULL,
            MINFRAC, MAXFRAC, 1.0, ctx->aggrrow, ctx->cutcoefs, &rhs, ctx->cutinds, &nnz, &eff, &rank, &islocal,
            &found) );
      if( found )
      {
         ++ctx->nmirsuccess;
         success = TRUE;
         cutnnz = nnz;
         cutrhs = rhs;
         cutrank = rank;
         cutislocal = islocal;
         cutefficacy = eff;
      }
   }

   /* lifted knapsack cover cut */
   if( ctx->sepadata->trycovers )
   {
      SCIP_Real eff = cutefficacy;
      SCIP_Real rhs;
      int nnz;
      int rank;
      SCIP_Bool islocal;
      SCIP_Bool found;

      SCIP_CALL( SCIPcalcKnapsackCover(scip, NULL, ctx->allowlocal, ctx->aggrrow, ctx->cutcoefs, &rhs, ctx->cutinds,
            &nnz, &eff, &rank, &islocal, &found) );
      if( found )
      {
         ++ctx->ncoverwins;
         success = TRUE;
         cutnnz = nnz;
         cutrhs = rhs;
         cutrank = rank;
         cutislocal = islocal;
         cutefficacy = eff;
      }
   }

   /* lifted flow cover cut */
   if( ctx->sepadata->trycovers )
   {
      SCIP_Real eff = cutefficacy;
      SCIP_Real rhs;
      int nnz;
      int rank;
      SCIP_Bool islocal;
      SCIP_Bool found;

      SCIP_CALL( SCIPcalcFlowCover(scip, NULL, POSTPROCESS, BOUNDSWITCH, ctx->allowlocal, ctx->aggrrow, ctx->cutcoefs,
            &rhs, ctx->cutinds, &nnz, &eff, &rank, &islocal, &found) );
      if( found )
      {
         ++ctx->ncoverwins;
         success = TRUE;
         cutnnz = nnz;
         cutrhs = rhs;
         cutrank = rank;
         cutislocal = islocal;
         cutefficacy = eff;
      }
   }

   /* cMIR heuristic over different scalings and complementations; since it is expensive, it is only applied to a
    * limited number of aggregations per separation call
    */
   if( ctx->sepadata->cmirheur && ctx->ncmirheur < ctx->sepadata->maxcmirheur )
   {
      SCIP_Real eff = cutefficacy;
      SCIP_Real rhs;
      int nnz;
      int rank;
      SCIP_Bool islocal;
      SCIP_Bool found;

      ++ctx->ncmirheur;
      SCIP_CALL( SCIPcutGenerationHeuristicCMIR(scip, NULL, POSTPROCESS, BOUNDSWITCH, USEVBDS, ctx->allowlocal,
            MAXTESTDELTA, NULL, NULL, MINFRAC, MAXFRAC, ctx->aggrrow, ctx->cutcoefs, &rhs, ctx->cutinds, &nnz, &eff,
            &rank, &islocal, &found) );
      if( found )
      {
         success = TRUE;
         cutnnz = nnz;
         cutrhs = rhs;
         cutrank = rank;
         cutislocal = islocal;
         cutefficacy = eff;
      }
   }

   SCIPdebugMsg(scip, "  variant %d: aggr rhs %.10g nnz %d nrows %d -> success %u eff %g\n", variant,
      SCIPaggrRowGetRhs(ctx->aggrrow), SCIPaggrRowGetNNz(ctx->aggrrow), SCIPaggrRowGetNRows(ctx->aggrrow), success,
      cutefficacy);

   if( success && cutnnz > 0 )
   {
      SCIP_VAR** vars = SCIPgetVars(scip);
      SCIP_ROW* cut;
      char cutname[SCIP_MAXSTRLEN];
      int i;

      ++ctx->nmireff;

      (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "mod%d_%" SCIP_LONGINT_FORMAT "_%d_%d", ctx->k, SCIPgetNLPs(scip),
         ctx->ncuts, variant);

      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, ctx->sepa, cutname, -SCIPinfinity(scip), cutrhs, cutislocal,
            FALSE, ctx->sepadata->dynamiccuts) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
      for( i = 0; i < cutnnz; ++i )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cut, vars[ctx->cutinds[i]], ctx->cutcoefs[i]) );
      }
      SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
      SCIProwChgRank(cut, cutrank);

      if( SCIPisCutEfficacious(scip, NULL, cut) && SCIPisCutNew(scip, cut) )
      {
         assert(ctx->ncuts < ctx->maxcuts);
         ctx->cuts[ctx->ncuts++] = cut;
      }
      else
      {
         SCIP_CALL( SCIPreleaseRow(scip, &cut) );
      }
   }

   return SCIP_OKAY;
}

/** processes a GF(k) solution: aggregates the rows with the derived weights and computes cuts */
static
SCIP_RETCODE processSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   MODKCUTCONTEXT*       ctx,                /**< cut context */
   int*                  solidx,             /**< unknowns with nonzero weight (will be sorted) */
   int*                  solw,               /**< weights in 1..k-1 */
   int                   nsol                /**< number of nonzero weights */
   )
{
   MODKROWSYSTEM* sys = ctx->sys;
   int k = ctx->k;
   void* hash;
   int i;

   if( nsol == 0 )
      return SCIP_OKAY;

   SCIPsortIntInt(solidx, solw, nsol);
   hash = hashWeights(solidx, solw, nsol);
   if( SCIPhashsetExists(ctx->usedweights, hash) )
      return SCIP_OKAY;
   SCIP_CALL( SCIPhashsetInsert(ctx->usedweights, SCIPblkmem(scip), hash) );

   ++ctx->nsolsk;

#ifndef NDEBUG
   /* check the solution against the integral system */
   {
      SCIP_Longint* chk;
      SCIP_Longint rhschk = 0;

      SCIP_CALL( SCIPallocClearBufferArray(scip, &chk, sys->neqs) );
      for( i = 0; i < nsol; ++i )
      {
         int r = solidx[i];
         int j;

         for( j = sys->beg[r]; j < sys->beg[r + 1]; ++j )
            chk[sys->eqind[j]] = (chk[sys->eqind[j]] + (SCIP_Longint)modK(sys->val[j], k) * solw[i]) % k;
         rhschk = (rhschk + (SCIP_Longint)modK(sys->rhs[r], k) * solw[i]) % k;
      }
      for( i = 0; i < sys->neqs; ++i )
         assert(chk[i] == 0);
      assert(rhschk == 1);
      SCIPfreeBufferArray(scip, &chk);
   }
#endif

   /* weights ((k-1) w mod k) / k: the aggregated right hand side has fractional part (k-1)/k */
   SCIPaggrRowClear(ctx->aggrrow);
   for( i = 0; i < nsol; ++i )
   {
      int r = solidx[i];
      int lambda = ((k - 1) * solw[i]) % k;

      if( lambda == 0 )
         continue;

      SCIP_CALL( SCIPaggrRowAddRow(scip, ctx->aggrrow, sys->rows[r], sys->sign[r] * sys->scale[r] * lambda / (SCIP_Real)k,
            sys->sign[r]) );
   }
   SCIP_CALL( computeCut(scip, ctx, 0) );

   /* weights w / k: the aggregated right hand side has fractional part 1/k */
   if( k != 2 && ctx->ncuts < ctx->maxcuts )
   {
      SCIPaggrRowClear(ctx->aggrrow);
      for( i = 0; i < nsol; ++i )
      {
         int r = solidx[i];

         SCIP_CALL( SCIPaggrRowAddRow(scip, ctx->aggrrow, sys->rows[r], sys->sign[r] * sys->scale[r] * solw[i] / (SCIP_Real)k,
               sys->sign[r]) );
      }
      SCIP_CALL( computeCut(scip, ctx, 1) );
   }
   SCIPaggrRowClear(ctx->aggrrow);

   return SCIP_OKAY;
}

/** enumerates basic solutions of the eliminated system and processes them */
static
SCIP_RETCODE gfkEnumerate(
   SCIP*                 scip,               /**< SCIP data structure */
   GFKSYSTEM*            gf,                 /**< eliminated GF(k) system */
   MODKCUTCONTEXT*       ctx,                /**< cut context */
   int                   maxsols             /**< maximal number of solutions to enumerate */
   )
{
   unsigned int k = (unsigned int) gf->k;
   int* swapidx;
   int* swapcol;
   int nswaps = 0;
   int swappos = 0;
   unsigned char* solval;
   int* solidx;
   int* solw;
   int nsol;
   int nenum = 0;
   int i;

   SCIP_CALL( SCIPallocBufferArray(scip, &swapidx, gf->nunks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &swapcol, gf->nunks) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &solval, gf->nunks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solidx, gf->nunks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solw, gf->nunks) );

   /* for each nonbasic unknown, schedule one basis swap into the last pivot row (in back substitution order, i.e.,
    * the pivot row with the largest index) that contains it; this keeps the factor triangular
    */
   for( i = gf->npiv - 1; i >= 0; --i )
   {
      int eq = gf->pivrow[i];
      int j;

      for( j = 0; j < gf->eqlen[eq]; ++j )
      {
         int col = gf->eqidx[eq][j];

         if( gf->status[col] != 0 )
            continue;

         gf->status[col] = -1;
         swapidx[nswaps] = i;
         swapcol[nswaps] = col;
         ++nswaps;
      }
   }

   for( ;; )
   {
      /* back substitution with all non-pivot unknowns set to zero */
      nsol = 0;
      for( i = gf->npiv - 1; i >= 0; --i )
      {
         int eq = gf->pivrow[i];
         int col = gf->pivcol[i];
         unsigned int sum = 0;
         unsigned int pivval = 0;
         unsigned int val;
         int j;

         for( j = 0; j < gf->eqlen[eq]; ++j )
         {
            int c = gf->eqidx[eq][j];

            if( c == col )
               pivval = gf->eqval[eq][j];
            else if( solval[c] != 0 )
               sum = (sum + (unsigned int)gf->eqval[eq][j] * solval[c]) % k;
         }
         assert(pivval != 0);

         val = ((gf->rhs[eq] + k - sum) * gfkInverse(pivval, gf->k)) % k;
         if( val != 0 )
         {
            solval[col] = (unsigned char) val;
            solidx[nsol] = col;
            solw[nsol] = (int) val;
            ++nsol;
         }
      }

      /* reset dense solution values before processing (which sorts the arrays) */
      for( i = 0; i < nsol; ++i )
         solval[solidx[i]] = 0;

      SCIP_CALL( processSolution(scip, ctx, solidx, solw, nsol) );
      ++nenum;

      if( nenum >= maxsols || ctx->ncuts >= ctx->maxcuts || swappos >= nswaps || SCIPisStopped(scip) )
         break;

      /* swap the next nonbasic unknown into the basis */
      {
         int pi = swapidx[swappos];
         int entering = swapcol[swappos];
         int leaving = gf->pivcol[pi];

         assert(gf->status[leaving] == 1);
         assert(eqFind(gf, gf->pivrow[pi], entering) >= 0);
         gf->pivcol[pi] = entering;
         gf->status[entering] = 1;
         gf->status[leaving] = 0;
         ++swappos;
      }
   }

   SCIPfreeBufferArray(scip, &solw);
   SCIPfreeBufferArray(scip, &solidx);
   SCIPfreeBufferArray(scip, &solval);
   SCIPfreeBufferArray(scip, &swapcol);
   SCIPfreeBufferArray(scip, &swapidx);

   return SCIP_OKAY;
}

/*
 * Row system construction
 */

/** builds the integral row system from the tight LP rows */
static
SCIP_RETCODE buildRowSystem(
   SCIP*                 scip,               /**< SCIP data structure */
   MODKROWSYSTEM*        sys,                /**< row system to fill */
   SCIP_Bool             allowlocal          /**< are local rows allowed? */
   )
{
   SCIP_COL** cols;
   SCIP_ROW** rows;
   SCIP_Real* tval;      /* transformed solution value of each LP column */
   SCIP_Real* shift;     /* bound used for complementing each LP column */
   signed char* csign;   /* +1: lower bound complement, -1: upper bound complement, 0: unusable */
   SCIP_Bool* cisint;
   int* eqofcol;
   SCIP_Real* scalevals;
   SCIP_Real feastol;
   SCIP_Real epsilon;
   int ncols;
   int nrows;
   int maxintrowlen;
   int entriessize;
   int i;

   BMSclearMemory(sys);

   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   if( ncols == 0 || nrows == 0 )
      return SCIP_OKAY;

   feastol = SCIPfeastol(scip);
   epsilon = SCIPepsilon(scip);
   maxintrowlen = MAXAGGRLEN(ncols);

   SCIP_CALL( SCIPallocBufferArray(scip, &tval, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &shift, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &csign, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cisint, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &eqofcol, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scalevals, ncols + 1) );

   /* compute the bound transformation of all LP columns */
   for( i = 0; i < ncols; ++i )
   {
      SCIP_VAR* var = SCIPcolGetVar(cols[i]);
      SCIP_Real x = SCIPcolGetPrimsol(cols[i]);
      SCIP_Real lb = allowlocal ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);
      SCIP_Real ub = allowlocal ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);
      SCIP_Bool lbinf = SCIPisInfinity(scip, -lb);
      SCIP_Bool ubinf = SCIPisInfinity(scip, ub);

      eqofcol[i] = -1;
      cisint[i] = SCIPcolIsIntegral(cols[i]);

      if( lbinf && ubinf )
      {
         csign[i] = 0;
         tval[i] = SCIPinfinity(scip);
         shift[i] = 0.0;
      }
      else if( !lbinf && (ubinf || x - lb <= ub - x) )
      {
         csign[i] = 1;
         tval[i] = MAX(x - lb, 0.0);
         shift[i] = lb;
      }
      else
      {
         csign[i] = -1;
         tval[i] = MAX(ub - x, 0.0);
         shift[i] = ub;
      }
   }

   entriessize = 0;
   for( i = 0; i < nrows; ++i )
      entriessize += SCIProwGetNLPNonz(rows[i]);

   entriessize = MAX(entriessize, 1);
   sys->rowssize = nrows;
   sys->entriessize = entriessize;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sys->rows, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sys->sign, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sys->scale, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sys->rhs, nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sys->beg, nrows + 1) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sys->eqind, entriessize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sys->val, entriessize) );
   sys->beg[0] = 0;

   for( i = 0; i < nrows; ++i )
   {
      SCIP_ROW* row = rows[i];
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      SCIP_Real activity;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real side;
      SCIP_Real transrhs;
      SCIP_Real intscale;
      SCIP_Longint intrhs;
      SCIP_Bool integral;
      SCIP_Bool skip;
      int rowlen;
      int nintentries;
      int nscalevals;
      int sign;
      int start;
      int j;

      if( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !allowlocal) )
         continue;

      rowlen = SCIProwGetNNonz(row);
      if( rowlen == 0 || SCIProwGetNLPNonz(row) != rowlen )
         continue;

      lhs = SCIProwGetLhs(row);
      rhs = SCIProwGetRhs(row);
      activity = SCIPgetRowLPActivity(scip, row);

      /* only use rows that are tight at the LP solution */
      if( !SCIPisInfinity(scip, rhs) && rhs - activity <= feastol )
         sign = 1;
      else if( !SCIPisInfinity(scip, -lhs) && activity - lhs <= feastol )
         sign = -1;
      else
         continue;

      integral = SCIProwIsIntegral(row);

      /* side as used by SCIPaggrRowAddRow() */
      if( sign == 1 )
      {
         side = rhs - SCIProwGetConstant(row);
         if( integral )
            side = SCIPfloor(scip, side);
      }
      else
      {
         side = lhs - SCIProwGetConstant(row);
         if( integral )
            side = SCIPceil(scip, side);
      }

      rowcols = SCIProwGetCols(row);
      rowvals = SCIProwGetVals(row);

      /* transform the row into the space of nonnegative bound distances; skip rows with continuous variables away from
       * their bounds and rows with free variables
       */
      skip = FALSE;
      transrhs = sign * side;
      nintentries = 0;
      for( j = 0; j < rowlen; ++j )
      {
         int c = SCIPcolGetLPPos(rowcols[j]);

         assert(c >= 0 && c < ncols);

         if( csign[c] == 0 || (!cisint[c] && tval[c] > feastol) )
         {
            skip = TRUE;
            break;
         }

         transrhs -= sign * rowvals[j] * shift[c];

         if( cisint[c] && tval[c] > feastol )
            ++nintentries;
      }

      if( skip || nintentries > maxintrowlen )
         continue;

      if( nintentries == 0 && REALABS(transrhs) <= epsilon )
         continue;

      /* compute the integral scalar */
      if( integral )
         intscale = 1.0;
      else
      {
         SCIP_Bool success;

         nscalevals = 0;
         for( j = 0; j < rowlen; ++j )
         {
            int c = SCIPcolGetLPPos(rowcols[j]);
            if( cisint[c] && tval[c] > feastol )
               scalevals[nscalevals++] = rowvals[j];
         }
         if( REALABS(transrhs) > epsilon )
            scalevals[nscalevals++] = transrhs;

         SCIP_CALL( SCIPcalcIntegralScalar(scalevals, nscalevals, -feastol, feastol, MAXDNOM, MAXINTSCALE, &intscale,
               &success) );
         if( !success || intscale <= 0.0 || intscale > MAXINTSCALE )
            continue;
      }

      /* the integral rhs has to fit into a long integer */
      if( REALABS(intscale * transrhs) > 1e+15 )
         continue;

      intrhs = (SCIP_Longint) SCIPround(scip, intscale * transrhs);

      /* store entries */
      start = sys->beg[sys->nrows];
      for( j = 0; j < rowlen; ++j )
      {
         int c = SCIPcolGetLPPos(rowcols[j]);
         SCIP_Real v;

         if( !cisint[c] || tval[c] <= feastol )
            continue;

         v = intscale * sign * csign[c] * rowvals[j];
         if( REALABS(v) > 1e+15 )
         {
            skip = TRUE;
            break;
         }

         if( eqofcol[c] < 0 )
            eqofcol[c] = sys->neqs++;

         sys->eqind[sys->nentries] = eqofcol[c];
         sys->val[sys->nentries] = (SCIP_Longint) SCIPround(scip, v);
         ++sys->nentries;
      }

      if( skip )
      {
         sys->nentries = start;
         continue;
      }

      sys->rows[sys->nrows] = row;
      sys->sign[sys->nrows] = sign;
      sys->scale[sys->nrows] = intscale;
      sys->rhs[sys->nrows] = intrhs;
      ++sys->nrows;
      sys->beg[sys->nrows] = sys->nentries;
   }

   SCIPfreeBufferArray(scip, &scalevals);
   SCIPfreeBufferArray(scip, &eqofcol);
   SCIPfreeBufferArray(scip, &cisint);
   SCIPfreeBufferArray(scip, &csign);
   SCIPfreeBufferArray(scip, &shift);
   SCIPfreeBufferArray(scip, &tval);

   return SCIP_OKAY;
}

/** frees the integral row system */
static
void freeRowSystem(
   SCIP*                 scip,               /**< SCIP data structure */
   MODKROWSYSTEM*        sys                 /**< row system */
   )
{
   if( sys->rows == NULL )
      return;

   SCIPfreeBlockMemoryArray(scip, &sys->val, sys->entriessize);
   SCIPfreeBlockMemoryArray(scip, &sys->eqind, sys->entriessize);
   SCIPfreeBlockMemoryArray(scip, &sys->beg, sys->rowssize + 1);
   SCIPfreeBlockMemoryArray(scip, &sys->rhs, sys->rowssize);
   SCIPfreeBlockMemoryArray(scip, &sys->scale, sys->rowssize);
   SCIPfreeBlockMemoryArray(scip, &sys->sign, sys->rowssize);
   SCIPfreeBlockMemoryArray(scip, &sys->rows, sys->rowssize);
}

/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyModk)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   SCIP_CALL( SCIPincludeSepaModk(scip) );

   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeModk)
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
SCIP_DECL_SEPAINITSOL(sepaInitsolModk)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIP_CALL( SCIPcreateRandom(scip, &sepadata->randnumgen, (unsigned int)sepadata->initseed, TRUE) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of separator */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolModk)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeRandom(scip, &sepadata->randnumgen);

   return SCIP_OKAY;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpModk)
{  /*lint --e{715}*/
   static const int primes[] = {2, 3, 5, 7};
   SCIP_SEPADATA* sepadata;
   MODKROWSYSTEM sys;
   MODKCUTCONTEXT ctx;
   SCIP_Bool cutoff = FALSE;
   int ncalls;
   int maxsepacuts;
   int nselected;
   int p;
   int i;

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

   maxsepacuts = depth == 0 ? sepadata->maxsepacutsroot : sepadata->maxsepacuts;
   if( maxsepacuts == 0 || sepadata->maxcutcands == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( buildRowSystem(scip, &sys, allowlocal) );

   SCIPdebugMsg(scip, "mod-k system: %d rows, %d columns, %d entries\n", sys.nrows, sys.neqs, sys.nentries);

   if( sys.nrows > 0 )
   {
      int nnonzerorhs = 0;

      for( i = 0; i < sys.nrows; ++i )
         nnonzerorhs += (sys.rhs[i] != 0);

      if( nnonzerorhs > 0 )
      {
         ctx.sepa = sepa;
         ctx.sepadata = sepadata;
         ctx.sys = &sys;
         ctx.allowlocal = allowlocal;
         ctx.ncuts = 0;
         ctx.nmirtried = 0;
         ctx.nmirsuccess = 0;
         ctx.nmireff = 0;
         ctx.ncoverwins = 0;
         ctx.ncmirheur = 0;
         ctx.maxcuts = sepadata->maxcutcands;
         SCIP_CALL( SCIPaggrRowCreate(scip, &ctx.aggrrow) );
         SCIP_CALL( SCIPallocBufferArray(scip, &ctx.cutcoefs, SCIPgetNVars(scip)) );
         SCIP_CALL( SCIPallocBufferArray(scip, &ctx.cutinds, SCIPgetNVars(scip)) );
         SCIP_CALL( SCIPallocBufferArray(scip, &ctx.cuts, ctx.maxcuts) );
         SCIP_CALL( SCIPhashsetCreate(&ctx.usedweights, SCIPblkmem(scip), MIN(sepadata->maxsols, 10000) + 1) );

         for( p = 0; p < 4; ++p )
         {
            GFKSYSTEM gf;
            SCIP_Bool solvable;
            int ncutsbefore = ctx.ncuts;

            if( primes[p] == 2 && !sepadata->usemod2 )
               continue;

            if( ctx.ncuts >= ctx.maxcuts || SCIPisStopped(scip) )
               break;

            ctx.k = primes[p];
            ctx.nsolsk = 0;
            SCIPhashsetRemoveAll(ctx.usedweights);

            SCIP_CALL( gfkCreate(scip, &gf, &sys, ctx.k) );
            SCIP_CALL( gfkEliminate(scip, &gf, &solvable) );

            if( solvable )
            {
               SCIP_CALL( gfkEnumerate(scip, &gf, &ctx, sepadata->maxsols) );
            }

            SCIPdebugMsg(scip, "k=%d: system %d x %d, %d pivots, solvable=%u, %d solutions, %d MIR tried, %d success, "
               "%d efficacious, %d cover wins, %d cuts\n", ctx.k, gf.neqs, gf.nunks, gf.npiv, solvable, ctx.nsolsk,
               ctx.nmirtried, ctx.nmirsuccess, ctx.nmireff, ctx.ncoverwins, ctx.ncuts - ncutsbefore);

            gfkFree(scip, &gf);

            if( sepadata->stopatfirstk && ctx.ncuts > ncutsbefore )
               break;
         }

         SCIPhashsetFree(&ctx.usedweights, SCIPblkmem(scip));

         /* select and add cuts: global cuts go to the pool, local cuts to the separation storage */
         if( ctx.ncuts > 0 )
         {
            SCIP_CALL( SCIPselectCutsHybrid(scip, ctx.cuts, NULL, sepadata->randnumgen, 1.0, 0.5, 0.1, 0.1, 0.0, 1.0,
                  0.0, 0.0, ctx.ncuts, 0, maxsepacuts, &nselected) );

            for( i = 0; i < ctx.ncuts; ++i )
            {
               if( i < nselected && !cutoff )
               {
                  if( SCIProwIsLocal(ctx.cuts[i]) )
                  {
                     SCIP_CALL( SCIPaddRow(scip, ctx.cuts[i], FALSE, &cutoff) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPaddPoolCut(scip, ctx.cuts[i]) );
                  }
               }
               SCIP_CALL( SCIPreleaseRow(scip, &ctx.cuts[i]) );
            }

            if( cutoff )
               *result = SCIP_CUTOFF;
            else if( nselected > 0 )
               *result = SCIP_SEPARATED;
         }

         SCIPfreeBufferArray(scip, &ctx.cuts);
         SCIPfreeBufferArray(scip, &ctx.cutinds);
         SCIPfreeBufferArray(scip, &ctx.cutcoefs);
         SCIPaggrRowFree(scip, &ctx.aggrrow);
      }
   }

   freeRowSystem(scip, &sys);

   return SCIP_OKAY;
}

/*
 * separator specific interface methods
 */

/** creates the mod-k cut separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaModk(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   BMSclearMemory(sepadata);

   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpModk, NULL, sepadata) );
   assert(sepa != NULL);

   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyModk) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeModk) );
   SCIP_CALL( SCIPsetSepaInitsol(scip, sepa, sepaInitsolModk) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolModk) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxrounds",
         "maximal number of mod-k separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxroundsroot",
         "maximal number of mod-k separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxsepacuts",
         "maximal number of mod-k cuts separated per separation round",
         &sepadata->maxsepacuts, FALSE, DEFAULT_MAXSEPACUTS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxsepacutsroot",
         "maximal number of mod-k cuts separated per separation round in the root node",
         &sepadata->maxsepacutsroot, FALSE, DEFAULT_MAXSEPACUTSROOT, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxcutcands",
         "maximal number of mod-k cut candidates collected per separation round",
         &sepadata->maxcutcands, FALSE, DEFAULT_MAXCUTCANDS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxsols",
         "maximal number of GF(k) solutions (row weight vectors) enumerated per modulus and separation round",
         &sepadata->maxsols, FALSE, DEFAULT_MAXSOLS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/initseed",
         "initial seed used for random tie-breaking in cut selection",
         &sepadata->initseed, FALSE, DEFAULT_INITSEED, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/usemod2",
         "should also mod-2 cuts be separated (usually done by the zerohalf separator)?",
         &sepadata->usemod2, TRUE, DEFAULT_USEMOD2, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/stopatfirstk",
         "should the separation stop at the first modulus that produced cuts?",
         &sepadata->stopatfirstk, TRUE, DEFAULT_STOPATFIRSTK, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/cmirheur",
         "should the cMIR scaling heuristic (different scalings and complementations) be tried on the aggregations?",
         &sepadata->cmirheur, TRUE, DEFAULT_CMIRHEUR, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/" SEPA_NAME "/maxcmirheur",
         "maximal number of aggregations per separation call on which the cMIR heuristic is tried",
         &sepadata->maxcmirheur, TRUE, DEFAULT_MAXCMIRHEUR, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/trycovers",
         "should lifted knapsack cover and flow cover cuts be tried on the mod-k aggregations as well?",
         &sepadata->trycovers, TRUE, DEFAULT_TRYCOVERS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/" SEPA_NAME "/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
