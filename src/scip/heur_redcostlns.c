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

/**@file   heur_redcostlns.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  reduced cost large neighborhood search: fixes variables by lurking reduced cost fixings and solves a sub-MIP
 * @author Max Buckley
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <limits.h>
#include <string.h>

#include "scip/heur_redcostlns.h"
#include "scip/pub_heur.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_copy.h"
#include "scip/scip_general.h"
#include "scip/scip_heur.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nodesel.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"

#define HEUR_NAME             "redcostlns"
#define HEUR_DESC             "LNS heuristic fixing variables by lurking reduced cost fixings at the root"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_LNS
#define HEUR_PRIORITY         -1090000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         0
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_TARGETFIXINGRATE 0.5   /**< fraction of unfixed integral variables at which fixing stops */
#define DEFAULT_MINFIXINGRATE    0.3   /**< minimal fraction of fixed integral variables to solve the sub-MIP */
#define DEFAULT_MINCANDRATE      0.1   /**< minimal fraction of unfixed integral variables with a lurking fixing */
#define DEFAULT_MAXPROPROUNDS    2     /**< maximal number of propagation rounds after each fixing (-1: unlimited) */
#define DEFAULT_MAXBACKTRACKS    100   /**< maximal number of fixings that are undone because of infeasibility */
#define DEFAULT_MAXNODES         500LL /**< maximal number of nodes of the sub-MIP */
#define DEFAULT_MINNODES         50LL  /**< minimal number of nodes required to start the sub-MIP */
#define DEFAULT_NODESOFS         500LL /**< number of nodes added to the contingent of the total nodes */
#define DEFAULT_NODESQUOT        0.1   /**< sub-MIP nodes in relation to the number of nodes of the original problem */
#define DEFAULT_MINIMPROVE       0.01  /**< factor by which the sub-MIP should improve the incumbent */
#define DEFAULT_COPYCUTS         TRUE  /**< should all active cuts of the cut pool be copied to the sub-MIP? */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Real             targetfixingrate;   /**< fraction of unfixed integral variables at which fixing stops */
   SCIP_Real             minfixingrate;      /**< minimal fraction of fixed integral variables to solve the sub-MIP */
   SCIP_Real             mincandrate;        /**< minimal fraction of unfixed integral variables with a lurking fixing */
   int                   maxproprounds;      /**< maximal number of propagation rounds after each fixing */
   int                   maxbacktracks;      /**< maximal number of fixings undone because of infeasibility */
   SCIP_Longint          maxnodes;           /**< maximal number of nodes of the sub-MIP */
   SCIP_Longint          minnodes;           /**< minimal number of nodes required to start the sub-MIP */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes */
   SCIP_Real             nodesquot;          /**< sub-MIP nodes in relation to the number of nodes of the original problem */
   SCIP_Real             minimprove;         /**< factor by which the sub-MIP should improve the incumbent */
   SCIP_Bool             copycuts;           /**< should all active cuts of the cut pool be copied to the sub-MIP? */
   SCIP_Longint          usednodes;          /**< nodes already used by the sub-MIPs */
};

/** collects the lurking reduced cost fixings of the current (root) LP solution
 *
 *  For a nonbasic integral column at its lower (upper) bound with reduced cost r, every solution with objective value
 *  below lpobj + |r| has the variable at that bound; the fixing is thus implied by reduced cost fixing once the cutoff
 *  bound drops below this threshold.
 */
static
SCIP_RETCODE collectLurkingFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            candvars,           /**< array to store candidate variables (size at least number of LP columns) */
   SCIP_Real*            candvals,           /**< array to store fixing values */
   SCIP_Real*            candthresholds,     /**< array to store required cutoff bounds */
   int*                  ncands              /**< pointer to store the number of candidates */
   )
{
   SCIP_COL** cols;
   SCIP_Real lpobj;
   int ncols;
   int c;

   *ncands = 0;
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   lpobj = SCIPgetLPObjval(scip);

   for( c = 0; c < ncols; ++c )
   {
      SCIP_VAR* var;
      SCIP_Real redcost;
      SCIP_BASESTAT basestat;

      var = SCIPcolGetVar(cols[c]);
      if( !SCIPvarIsIntegral(var) )
         continue;
      if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         continue;

      basestat = SCIPcolGetBasisStatus(cols[c]);
      if( basestat != SCIP_BASESTAT_LOWER && basestat != SCIP_BASESTAT_UPPER )
         continue;

      redcost = SCIPgetColRedcost(scip, cols[c]);

      if( basestat == SCIP_BASESTAT_LOWER && SCIPisDualfeasPositive(scip, redcost) )
      {
         candvals[*ncands] = SCIPvarGetLbLocal(var);
         candthresholds[*ncands] = lpobj + redcost;
      }
      else if( basestat == SCIP_BASESTAT_UPPER && SCIPisDualfeasNegative(scip, redcost) )
      {
         candvals[*ncands] = SCIPvarGetUbLocal(var);
         candthresholds[*ncands] = lpobj - redcost;
      }
      else
         continue;

      if( SCIPisInfinity(scip, REALABS(candvals[*ncands])) )
         continue;

      candvars[*ncands] = var;
      ++(*ncands);
   }

   return SCIP_OKAY;
}

/** sets up and solves the sub-MIP on the current probing node; copies the local bounds, i.e., the fixings */
static
SCIP_RETCODE solveSubMIP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< the heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Longint          nstallnodes,        /**< stall node limit for the sub-MIP */
   SCIP_Real             lowerbound,         /**< lower bound used for the minimal improvement */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP* subscip;
   SCIP_VAR** vars;
   SCIP_VAR** subvars;
   SCIP_HASHMAP* varmap;
   SCIP_Bool valid;
   SCIP_Bool success;
   SCIP_RETCODE retstat;
   int nvars;
   int i;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(subscip), nvars) );

   /* copy with local bounds, so that the probing fixings are part of the sub-MIP */
   SCIP_CALL( SCIPcopy(scip, subscip, varmap, NULL, "_redcostlns", FALSE, FALSE, FALSE, FALSE, TRUE, &valid) );

   if( heurdata->copycuts )
   {
      SCIP_CALL( SCIPcopyCuts(scip, subscip, varmap, NULL, FALSE, NULL) );
   }

   for( i = 0; i < nvars; ++i )
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmap, vars[i]);
   SCIPhashmapFree(&varmap);

   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );
#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#else
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

   SCIP_CALL( SCIPcopyLimits(scip, subscip) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", heurdata->maxnodes) );

   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_FAST, TRUE) );
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }
   SCIP_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

   /* require an improvement over the incumbent */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real upperbound;
      SCIP_Real cutoffbound;

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
      if( !SCIPisInfinity(scip, -lowerbound) )
         cutoffbound = (1.0 - heurdata->minimprove) * SCIPgetUpperbound(scip) + heurdata->minimprove * lowerbound;
      else if( SCIPgetUpperbound(scip) >= 0.0 )
         cutoffbound = (1.0 - heurdata->minimprove) * SCIPgetUpperbound(scip);
      else
         cutoffbound = (1.0 + heurdata->minimprove) * SCIPgetUpperbound(scip);
      cutoffbound = MIN(upperbound, cutoffbound);
      SCIP_CALL( SCIPsetObjlimit(subscip, cutoffbound) );
   }

   /* errors in the sub-MIP should not abort the main solve in optimized mode */
   retstat = SCIPsolve(subscip);
   if( retstat != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retstat );
#endif
      SCIPwarningMessage(scip, "Error while solving sub-MIP in " HEUR_NAME " heuristic; sub-SCIP terminated with code <%d>\n",
         retstat);
   }
   else
   {
      SCIP_CALL( SCIPtranslateSubSols(scip, subscip, heur, subvars, &success, NULL) );
      if( success )
         *result = SCIP_FOUNDSOL;
   }

   heurdata->usednodes += SCIPgetNNodes(subscip);

   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRedcostlns)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   SCIP_CALL( SCIPincludeHeurRedcostlns(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRedcostlns)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitRedcostlns)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   heurdata->usednodes = 0;

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRedcostlns)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** candvars;
   SCIP_Real* candvals;
   SCIP_Real* candthresholds;
   SCIP_Real lowerbound;
   SCIP_Real fixingrate;
   SCIP_Longint nstallnodes;
   SCIP_Bool cutoff;
   SCIP_Bool valid;
   int nunfixedint;
   int ncands;
   int nbacktracks;
   int ncols;
   int c;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* we need an optimal basic LP solution at the root */
   if( SCIPgetDepth(scip) != 0 || !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL
      || !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   if( SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) == 0 )
      return SCIP_OKAY;

   /* node budget, as in other LNS heuristics */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur) + 1.0) / (SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);
   nstallnodes += heurdata->nodesofs;
   nstallnodes -= heurdata->usednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);
   if( nstallnodes < heurdata->minnodes )
      return SCIP_OKAY;

   SCIP_CALL( SCIPcheckCopyLimits(scip, &valid) );
   if( !valid )
      return SCIP_OKAY;

   nunfixedint = SCIPgetNPseudoBranchCands(scip);
   if( nunfixedint == 0 )
      return SCIP_OKAY;

   ncols = SCIPgetNLPCols(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &candvars, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candvals, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candthresholds, ncols) );

   SCIP_CALL( collectLurkingFixings(scip, candvars, candvals, candthresholds, &ncands) );

   SCIPdebugMsg(scip, "%d lurking reduced cost fixings for %d unfixed integral variables\n", ncands, nunfixedint);

   if( ncands < heurdata->mincandrate * nunfixedint )
      goto FREEARRAYS;

   *result = SCIP_DIDNOTFIND;

   /* most robust fixings (largest required cutoff, i.e., largest reduced cost) first */
   SCIPsortDownRealRealPtr(candthresholds, candvals, (void**)candvars, ncands);

   lowerbound = SCIPgetLowerbound(scip);

   SCIP_CALL( SCIPstartProbing(scip) );

   cutoff = FALSE;
   nbacktracks = 0;
   fixingrate = 0.0;

   for( c = 0; c < ncands && fixingrate < heurdata->targetfixingrate && !SCIPisStopped(scip); ++c )
   {
      SCIP_VAR* var = candvars[c];

      /* fixings that are only implied by cutoffs below the dual bound cannot be part of an improving solution */
      if( candthresholds[c] <= lowerbound + SCIPfeastol(scip) )
         break;

      /* skip variables fixed or tightened away from the candidate value by propagation */
      if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var))
         || SCIPisFeasLT(scip, candvals[c], SCIPvarGetLbLocal(var))
         || SCIPisFeasGT(scip, candvals[c], SCIPvarGetUbLocal(var)) )
         continue;

      SCIP_CALL( SCIPnewProbingNode(scip) );
      SCIP_CALL( SCIPfixVarProbing(scip, var, candvals[c]) );
      SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, &cutoff, NULL) );

      if( cutoff )
      {
         /* undo this fixing and continue with the next candidate */
         SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip) - 1) );
         cutoff = FALSE;
         if( ++nbacktracks > heurdata->maxbacktracks )
            break;
         continue;
      }

      fixingrate = 1.0 - SCIPgetNPseudoBranchCands(scip) / (SCIP_Real)nunfixedint;
   }

   fixingrate = 1.0 - SCIPgetNPseudoBranchCands(scip) / (SCIP_Real)nunfixedint;

   SCIPdebugMsg(scip, "fixing rate %.3f after %d candidates, %d backtracks\n", fixingrate, c, nbacktracks);

   if( fixingrate >= heurdata->minfixingrate && !SCIPisStopped(scip) )
   {
      SCIP_CALL( solveSubMIP(scip, heur, heurdata, nstallnodes, lowerbound, result) );
   }

   SCIP_CALL( SCIPendProbing(scip) );

FREEARRAYS:
   SCIPfreeBufferArray(scip, &candthresholds);
   SCIPfreeBufferArray(scip, &candvals);
   SCIPfreeBufferArray(scip, &candvars);

   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the reduced cost LNS primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRedcostlns(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   heurdata->usednodes = 0;

   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ,
         HEUR_FREQOFS, HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRedcostlns, heurdata) );
   assert(heur != NULL);

   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRedcostlns) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRedcostlns) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRedcostlns) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/targetfixingrate",
         "fraction of unfixed integral variables at which fixing stops",
         &heurdata->targetfixingrate, TRUE, DEFAULT_TARGETFIXINGRATE, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minfixingrate",
         "minimal fraction of fixed integral variables to solve the sub-MIP",
         &heurdata->minfixingrate, TRUE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/mincandrate",
         "minimal fraction of unfixed integral variables with a lurking reduced cost fixing",
         &heurdata->mincandrate, TRUE, DEFAULT_MINCANDRATE, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxproprounds",
         "maximal number of propagation rounds after each fixing (-1: unlimited)",
         &heurdata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxbacktracks",
         "maximal number of fixings that are undone because of infeasibility",
         &heurdata->maxbacktracks, TRUE, DEFAULT_MAXBACKTRACKS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/maxnodes",
         "maximal number of nodes of the sub-MIP",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/minnodes",
         "minimal number of nodes required to start the sub-MIP",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/nodesquot",
         "sub-MIP nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/minimprove",
         "factor by which the sub-MIP should improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/copycuts",
         "should all active cuts of the cut pool be copied to the sub-MIP?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
