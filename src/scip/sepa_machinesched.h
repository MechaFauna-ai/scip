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

/**@file   sepa_machinesched.h
 * @ingroup SEPARATORS
 * @brief  separator for release-date cuts of single-machine scheduling substructures
 *
 * Detects big-M disjunctions M y + s_i - s_j >= p between start times s_i, s_j of jobs (with finite lower bounds,
 * the release dates) controlled by a binary variable y, and greedily searches for a set of jobs in which every pair
 * of jobs is ordered by such a disjunction (a single machine). For each job j of such a set N,
 *
 *    s_j >= min{ r_i : i in N } + sum_{i in N \ j} p_i [i before j]
 *
 * is valid, where p_i is the minimal distance of job i to any successor in N and [i before j] is the literal of the
 * binary variable that orders i before j.
 *
 * Based on HiGHS (HighsMachineSchedSeparator, https://github.com/ERGO-Code/HiGHS/pull/2990) and
 *  G. L. Nemhauser and M. W. P. Savelsbergh, "A cutting plane algorithm for the single machine scheduling problem with
 *  release times", 1992; R. L. H. Wijfjes, "Recognition and Exploitation of Single-Machine Scheduling Subproblems in
 *  Mixed Integer Programs", MSc thesis, 2022.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_MACHINESCHED_H__
#define __SCIP_SEPA_MACHINESCHED_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the machinesched separator and includes it in SCIP
 *
 *  @ingroup SeparatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSepaMachinesched(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
