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
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_redcostlns.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  root reduced-cost probing and LNS heuristic
 *
 * Collects integer variables with strong root LP reduced costs, orders their
 * lurking fixings by descending implied cutoff bound, propagates them in
 * probing mode (lifting the root dual bound whenever propagation hits
 * infeasibility), and solves a tightly budgeted LNS sub-MIP if the propagated
 * fixing rate reaches the target threshold.
 */

#include <assert.h>
#include <string.h>

#include "blockmemshell/memory.h"
#include "scip/heur_redcostlns.h"
#include "scip/heuristics.h"
#include "scip/pub_event.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/scip_branch.h"
#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_event.h"
#include "scip/scip_exact.h"
#include "scip/scip_general.h"
#include "scip/scip_heur.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nodesel.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_pricer.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"

#define HEUR_NAME             "redcostlns"
#define HEUR_DESC             "root reduced-cost probing and LNS heuristic"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_LNS
#define HEUR_PRIORITY         -1050000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         0
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      TRUE

#define DEFAULT_MINFIXINGRATE 0.30
#define DEFAULT_TARGETFIXRATE 0.50
#define DEFAULT_MINLURKRATIO  0.10
#define DEFAULT_MINIMPROVE    0.01
#define DEFAULT_MAXNODES      500LL
#define DEFAULT_STALLNODES    200LL
#define DEFAULT_LPLIMFAC      2.0
#define DEFAULT_COPYCUTS      TRUE
#define DEFAULT_LIFTDUALBOUND TRUE

#define EVENTHDLR_NAME        "Redcostlns"
#define EVENTHDLR_DESC        "LP event handler for " HEUR_NAME " heuristic"

#define REDCOSTLNS_CALL(x)    do { if( (retcode = (x)) != SCIP_OKAY ) goto TERMINATE; } while( FALSE )

struct SCIP_HeurData
{
   SCIP_Real             minfixingrate;
   SCIP_Real             targetfixrate;
   SCIP_Real             minlurkratio;
   SCIP_Real             minimprove;
   SCIP_Longint          maxnodes;
   SCIP_Longint          stallnodes;
   SCIP_Real             lplimfac;
   SCIP_Real             nodelimit;
   SCIP_Bool             copycuts;
   SCIP_Bool             liftdualbound;
};

typedef struct
{
   SCIP_Real             lurkcutoff;
   SCIP_VAR*             var;
   SCIP_Real             fixval;
   SCIP_Bool             fixupper;
} LURK_BOUND;

static
SCIP_DECL_SORTPTRCOMP(compLurkDesc)
{
   const LURK_BOUND* a = (const LURK_BOUND*)elem1;
   const LURK_BOUND* b = (const LURK_BOUND*)elem2;

   if( a->lurkcutoff > b->lurkcutoff )
      return -1;
   if( a->lurkcutoff < b->lurkcutoff )
      return 1;
   if( SCIPvarGetProbindex(a->var) < SCIPvarGetProbindex(b->var) )
      return -1;
   if( SCIPvarGetProbindex(a->var) > SCIPvarGetProbindex(b->var) )
      return 1;
   return 0;
}

static
SCIP_DECL_EVENTEXEC(eventExecRedcostlns)
{
   SCIP_HEURDATA* heurdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_LPSOLVED);
   SCIP_STRINGEQ( SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME, SCIP_INVALIDCALL );

   heurdata = (SCIP_HEURDATA*)eventdata;

   if( SCIPgetNLPs(scip) > heurdata->lplimfac * heurdata->nodelimit )
   {
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE solveSubscip(
   SCIP*                 scip,
   SCIP*                 subscip,
   SCIP_HEUR*            heur,
   SCIP_HEURDATA*        heurdata,
   SCIP_VAR**            fixedvars,
   SCIP_Real*            fixedvals,
   int                   nfixedvars,
   SCIP_RESULT*          result
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** subvars = NULL;
   SCIP_HASHMAP* varmapfw = NULL;
   SCIP_EVENTHDLR* eventhdlr = NULL;
   SCIP_Real allfixingrate;
   SCIP_Bool success;
   SCIP_Bool caughtevent = FALSE;
   SCIP_RETCODE retcode = SCIP_OKAY;
   SCIP_RETCODE subretcode;
   int nvars;
   int i;

   REDCOSTLNS_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   REDCOSTLNS_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );

   REDCOSTLNS_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapfw, "redcostlns", fixedvars, fixedvals,
         nfixedvars, FALSE, heurdata->copycuts, &success, NULL) );

   REDCOSTLNS_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecRedcostlns, NULL) );
   if( eventhdlr == NULL )
   {
      retcode = SCIP_PLUGINNOTFOUND;
      goto TERMINATE;
   }

   REDCOSTLNS_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   for( i = 0; i < nvars; ++i )
      subvars[i] = (SCIP_VAR*)SCIPhashmapGetImage(varmapfw, vars[i]);
   SCIPhashmapFree(&varmapfw);

   REDCOSTLNS_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );
   REDCOSTLNS_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   REDCOSTLNS_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );

   REDCOSTLNS_CALL( SCIPcopyLimits(scip, subscip) );
   heurdata->nodelimit = (SCIP_Real)heurdata->maxnodes;
   REDCOSTLNS_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", heurdata->stallnodes) );
   REDCOSTLNS_CALL( SCIPsetLongintParam(subscip, "limits/nodes", heurdata->maxnodes) );
   REDCOSTLNS_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   REDCOSTLNS_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );
   REDCOSTLNS_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   if( SCIPfindNodesel(subscip, "estimate") != NULL && !SCIPisParamFixed(subscip, "nodeselection/estimate/stdpriority") )
   {
      REDCOSTLNS_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX / 4) );
   }
   if( SCIPfindBranchrule(subscip, "inference") != NULL && !SCIPisParamFixed(subscip, "branching/inference/priority") )
   {
      REDCOSTLNS_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX / 4) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/enable") )
   {
      REDCOSTLNS_CALL( SCIPsetBoolParam(subscip, "conflict/enable", TRUE) );
   }
   if( !SCIPisParamFixed(subscip, "conflict/useboundlp") )
   {
      REDCOSTLNS_CALL( SCIPsetCharParam(subscip, "conflict/useboundlp", 'o') );
   }
   if( !SCIPisParamFixed(subscip, "conflict/maxstoresize") )
   {
      REDCOSTLNS_CALL( SCIPsetIntParam(subscip, "conflict/maxstoresize", 100) );
   }
   REDCOSTLNS_CALL( SCIPsetBoolParam(subscip, "lp/checkdualfeas", FALSE) );

   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);
      SCIP_Real cutoff;

      if( !SCIPisInfinity(scip, -SCIPgetLowerbound(scip)) )
         cutoff = (1.0 - heurdata->minimprove) * SCIPgetUpperbound(scip)
            + heurdata->minimprove * SCIPgetLowerbound(scip);
      else
         cutoff = SCIPgetUpperbound(scip) >= 0.0
            ? (1.0 - heurdata->minimprove) * SCIPgetUpperbound(scip)
            : (1.0 + heurdata->minimprove) * SCIPgetUpperbound(scip);

      cutoff = MIN(upperbound, cutoff);
      REDCOSTLNS_CALL( SCIPsetObjlimit(subscip, cutoff) );
   }

   subretcode = SCIPpresolve(subscip);
   if( subretcode != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Error while presolving subproblem in %s heuristic; sub-SCIP terminated with code <%d>\n",
         HEUR_NAME, subretcode);
      retcode = subretcode;
      goto TERMINATE;
   }

   allfixingrate = (SCIPgetNOrigVars(subscip) > 0)
      ? (SCIPgetNOrigVars(subscip) - SCIPgetNVars(subscip)) / (SCIP_Real)SCIPgetNOrigVars(subscip)
      : 0.0;
   allfixingrate = MAX(allfixingrate, 0.0);

   if( allfixingrate >= heurdata->minfixingrate / 2.0 )
   {
      SCIP_SOL** subsols;
      int nsubsols;

      REDCOSTLNS_CALL( SCIPtransformProb(subscip) );
      REDCOSTLNS_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*)heurdata, NULL) );
      caughtevent = TRUE;
      subretcode = SCIPsolve(subscip);
      if( subretcode != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Error while solving subproblem in %s heuristic; sub-SCIP terminated with code <%d>\n",
            HEUR_NAME, subretcode);
         retcode = subretcode;
         goto TERMINATE;
      }
      REDCOSTLNS_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*)heurdata, -1) );
      caughtevent = FALSE;

      REDCOSTLNS_CALL( SCIPmergeVariableStatistics(subscip, scip, subvars, vars, nvars) );

      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      success = FALSE;
      for( i = 0; i < nsubsols && !success; ++i )
      {
         SCIP_SOL* newsol;

         REDCOSTLNS_CALL( SCIPtranslateSubSol(scip, subscip, subsols[i], heur, subvars, &newsol) );
         REDCOSTLNS_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );
         if( success )
            *result = SCIP_FOUNDSOL;
      }
   }

TERMINATE:
   if( caughtevent )
   {
      SCIP_RETCODE dropret = SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*)heurdata, -1);
      if( retcode == SCIP_OKAY )
         retcode = dropret;
   }
   SCIPfreeBufferArrayNull(scip, &subvars);
   if( varmapfw != NULL )
      SCIPhashmapFree(&varmapfw);
   return retcode;
}

static
SCIP_DECL_HEURCOPY(heurCopyRedcostlns)
{
   assert(scip != NULL);
   assert(heur != NULL);
   SCIP_STRINGEQ( SCIPheurGetName(heur), HEUR_NAME, SCIP_INVALIDCALL );

   SCIP_CALL( SCIPincludeHeurRedcostlns(scip) );
   return SCIP_OKAY;
}

static
SCIP_DECL_HEURFREE(heurFreeRedcostlns)
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);
   return SCIP_OKAY;
}

static
SCIP_DECL_HEUREXEC(heurExecRedcostlns)
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** vars;
   LURK_BOUND* lurk = NULL;
   LURK_BOUND** lurkptrs = NULL;
   SCIP_VAR** fixedvars = NULL;
   SCIP_Real* fixedvals = NULL;
   SCIP_Real cutoffbound;
   SCIP_Real globallb;
   SCIP_Real liftedlb;
   SCIP_Bool canlift;
   SCIP_Bool success;
   SCIP_RETCODE retcode = SCIP_OKAY;
   int nbinvars;
   int nintvars;
   int ntotalint;
   int nlurk;
   int nfixedvars;
   int nfixedint;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( nodeinfeasible || SCIPgetDepth(scip) != 0 || SCIPgetSubscipDepth(scip) > 0 )
      return SCIP_OKAY;

   if( !SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL
      || !SCIPisLPDualReliable(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   cutoffbound = SCIPgetCutoffbound(scip);
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), cutoffbound) )
      return SCIP_OKAY;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );
   ntotalint = nbinvars + nintvars;
   if( ntotalint < 4 )
      return SCIP_OKAY;

   globallb = SCIPgetLocalLowerbound(scip);
   canlift = heurdata->liftdualbound && !SCIPisExact(scip) && SCIPallowStrongDualReds(scip);

   REDCOSTLNS_CALL( SCIPallocBufferArray(scip, &lurk, ntotalint) );
   REDCOSTLNS_CALL( SCIPallocBufferArray(scip, &lurkptrs, ntotalint) );
   nlurk = 0;

   for( i = 0; i < ntotalint; ++i )
   {
      SCIP_VAR* var = vars[i];
      SCIP_Real lb = SCIPvarGetLbLocal(var);
      SCIP_Real ub = SCIPvarGetUbLocal(var);
      SCIP_Real rc;
      SCIP_Real rsol;
      SCIP_Real rlpobj;
      SCIP_Real delta;
      SCIP_Real lcut;
      SCIP_Bool fixupper;
      SCIP_Real fixval;

      if( ub - lb < 0.5 )
         continue;

      rlpobj = SCIPvarGetBestRootLPObjval(var);
      if( rlpobj == SCIP_INVALID || SCIPisInfinity(scip, -rlpobj) || SCIPisInfinity(scip, rlpobj) )
         continue;

      rc = SCIPvarGetBestRootRedcost(var);
      if( SCIPisDualfeasZero(scip, rc) )
         continue;

      rsol = SCIPvarGetBestRootSol(var);
      if( rsol == SCIP_INVALID )
         continue;

      if( SCIPisDualfeasPositive(scip, rc) && !SCIPisInfinity(scip, -lb) && SCIPisFeasEQ(scip, rsol, lb) )
      {
         delta = MIN(1.0, (lb + 1.0) - rsol);
         if( delta <= 0.5 )
            continue;
         lcut = rlpobj + rc * delta;
         fixupper = TRUE;
         fixval = lb;
      }
      else if( SCIPisDualfeasNegative(scip, rc) && !SCIPisInfinity(scip, ub) && SCIPisFeasEQ(scip, rsol, ub) )
      {
         delta = MIN(1.0, rsol - (ub - 1.0));
         if( delta <= 0.5 )
            continue;
         lcut = rlpobj - rc * delta;
         fixupper = FALSE;
         fixval = ub;
      }
      else
         continue;

      if( lcut <= globallb + SCIPfeastol(scip) || lcut >= cutoffbound )
         continue;

      lurk[nlurk].lurkcutoff = lcut;
      lurk[nlurk].var = var;
      lurk[nlurk].fixval = fixval;
      lurk[nlurk].fixupper = fixupper;
      lurkptrs[nlurk] = &lurk[nlurk];
      ++nlurk;
   }

   if( (SCIP_Real)nlurk < heurdata->minlurkratio * (SCIP_Real)ntotalint )
      goto TERMINATE;

   *result = SCIP_DIDNOTFIND;
   SCIPsortPtr((void**)lurkptrs, compLurkDesc, nlurk);

   REDCOSTLNS_CALL( SCIPallocBufferArray(scip, &fixedvars, ntotalint) );
   REDCOSTLNS_CALL( SCIPallocBufferArray(scip, &fixedvals, ntotalint) );
   nfixedvars = 0;
   nfixedint = 0;
   liftedlb = globallb;

   /* propagate sorted lurking fixings in probing mode */
   REDCOSTLNS_CALL( SCIPstartProbing(scip) );

   for( i = 0; i < nlurk && !SCIPisStopped(scip); ++i )
   {
      LURK_BOUND* lbnd = lurkptrs[i];
      SCIP_VAR* var = lbnd->var;
      SCIP_Bool cutoff = FALSE;
      SCIP_Real currlb = SCIPvarGetLbLocal(var);
      SCIP_Real currub = SCIPvarGetUbLocal(var);
      int pdepth;

      if( lbnd->lurkcutoff <= liftedlb + SCIPfeastol(scip) )
         break;

      if( lbnd->fixupper )
      {
         if( currub <= lbnd->fixval + 0.5 )
            continue;
         if( currlb > lbnd->fixval + 0.5 )
         {
            /* direct contradiction with previously propagated bounds */
            if( canlift )
            {
               SCIP_Real safelb = lbnd->lurkcutoff - 1e-5 * MAX(1.0, REALABS(lbnd->lurkcutoff));
               if( safelb > liftedlb && safelb < cutoffbound )
                  liftedlb = safelb;
            }
            break;
         }
      }
      else
      {
         if( currlb >= lbnd->fixval - 0.5 )
            continue;
         if( currub < lbnd->fixval - 0.5 )
         {
            if( canlift )
            {
               SCIP_Real safelb = lbnd->lurkcutoff - 1e-5 * MAX(1.0, REALABS(lbnd->lurkcutoff));
               if( safelb > liftedlb && safelb < cutoffbound )
                  liftedlb = safelb;
            }
            break;
         }
      }

      if( SCIPgetDepth(scip) >= SCIP_MAXTREEDEPTH )
         break;

      pdepth = SCIPgetProbingDepth(scip);
      REDCOSTLNS_CALL( SCIPnewProbingNode(scip) );

      if( lbnd->fixupper )
      {
         REDCOSTLNS_CALL( SCIPchgVarUbProbing(scip, var, lbnd->fixval) );
      }
      else
      {
         REDCOSTLNS_CALL( SCIPchgVarLbProbing(scip, var, lbnd->fixval) );
      }

      REDCOSTLNS_CALL( SCIPpropagateProbing(scip, 2, &cutoff, NULL) );
      if( cutoff )
      {
         if( canlift )
         {
            SCIP_Real safelb = lbnd->lurkcutoff - 1e-5 * MAX(1.0, REALABS(lbnd->lurkcutoff));
            if( safelb > liftedlb && safelb < cutoffbound )
               liftedlb = safelb;
         }
         REDCOSTLNS_CALL( SCIPbacktrackProbing(scip, pdepth) );
         break;
      }

      /* check fixing rate every 16 applied fixings */
      if( (i & 15) == 15 )
      {
         int k;
         int cnt = 0;

         for( k = 0; k < ntotalint; ++k )
         {
            if( SCIPvarGetUbLocal(vars[k]) - SCIPvarGetLbLocal(vars[k]) < 0.5 )
               ++cnt;
         }
         if( (SCIP_Real)cnt >= heurdata->targetfixrate * (SCIP_Real)ntotalint )
            break;
      }
   }

   /* extract propagated integer fixings for LNS neighborhood */
   for( i = 0; i < ntotalint; ++i )
   {
      SCIP_Real llb = SCIPvarGetLbLocal(vars[i]);
      SCIP_Real lub = SCIPvarGetUbLocal(vars[i]);

      if( lub - llb < 0.5 )
      {
         ++nfixedint;
         if( SCIPvarGetUbGlobal(vars[i]) - SCIPvarGetLbGlobal(vars[i]) >= 0.5 )
         {
            fixedvars[nfixedvars] = vars[i];
            fixedvals[nfixedvars] = SCIPfeasFloor(scip, llb + 0.5);
            ++nfixedvars;
         }
      }
   }

   REDCOSTLNS_CALL( SCIPendProbing(scip) );

   /* apply any proven global root dual bound lift (sound at depth 0 where local == global) */
   if( canlift && liftedlb > globallb + SCIPfeastol(scip) && liftedlb < cutoffbound )
   {
      assert(SCIPgetDepth(scip) == 0);
      REDCOSTLNS_CALL( SCIPupdateLocalLowerbound(scip, liftedlb) );
   }

   /* solve fast LNS sub-MIP if fixing rate is sufficient */
   if( nfixedvars > 0 && (SCIP_Real)nfixedint >= heurdata->minfixingrate * (SCIP_Real)ntotalint
      && !SCIPisStopped(scip) && SCIPisLT(scip, SCIPgetLocalLowerbound(scip), SCIPgetCutoffbound(scip)) )
   {
      REDCOSTLNS_CALL( SCIPcheckCopyLimits(scip, &success) );
      if( success )
      {
         SCIP* subscip;
         SCIP_RETCODE subret;

         REDCOSTLNS_CALL( SCIPcreate(&subscip) );
         subret = solveSubscip(scip, subscip, heur, heurdata, fixedvars, fixedvals, nfixedvars, result);
         REDCOSTLNS_CALL( SCIPfree(&subscip) );
         REDCOSTLNS_CALL( subret );
      }
   }

TERMINATE:
   if( SCIPinProbing(scip) )
   {
      SCIP_RETCODE endret = SCIPendProbing(scip);
      if( retcode == SCIP_OKAY )
         retcode = endret;
   }
   SCIPfreeBufferArrayNull(scip, &fixedvals);
   SCIPfreeBufferArrayNull(scip, &fixedvars);
   SCIPfreeBufferArrayNull(scip, &lurkptrs);
   SCIPfreeBufferArrayNull(scip, &lurk);

   return retcode;
}

/** creates the root reduced-cost LNS primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRedcostlns(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR,
         HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS, HEUR_MAXDEPTH, HEUR_TIMING,
         HEUR_USESSUBSCIP, heurExecRedcostlns, heurdata) );

   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRedcostlns) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRedcostlns) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/minfixingrate",
         "minimum fraction of integer variables fixed to trigger sub-MIP",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/targetfixrate",
         "target fraction of integer variables fixed during reduced-cost probing",
         &heurdata->targetfixrate, FALSE, DEFAULT_TARGETFIXRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/minlurkratio",
         "minimum fraction of integer variables with nonzero root reduced cost",
         &heurdata->minlurkratio, FALSE, DEFAULT_MINLURKRATIO, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/minimprove",
         "factor by which the heuristic should at least improve the incumbent",
         &heurdata->minimprove, FALSE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip,
         "heuristics/" HEUR_NAME "/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip,
         "heuristics/" HEUR_NAME "/stallnodes",
         "number of stalling nodes for the subproblem",
         &heurdata->stallnodes, TRUE, DEFAULT_STALLNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/" HEUR_NAME "/lplimfac",
         "factor by which the LP limit depends on the node limit",
         &heurdata->lplimfac, TRUE, DEFAULT_LPLIMFAC, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/copycuts",
         "should all active cuts from cutpool be copied to sub-SCIP?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/" HEUR_NAME "/liftdualbound",
         "should probing infeasibility lift the root dual bound?",
         &heurdata->liftdualbound, TRUE, DEFAULT_LIFTDUALBOUND, NULL, NULL) );

   return SCIP_OKAY;
}
