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

/**@file   sort.c
 * @brief  unit tests for some sort implementations
 * @author Christoph Schubert
 */

#include<stdio.h>

#include "scip/pub_misc.h"
#include "scip/scip.h"

#include "include/scip_test.h"

#include <math.h>

/* helper methods */

/** compares two entries of an int array, returns 0 if both are equal, a negative number if the entry at position ind2
 *  is larger than the one at position ind1, and a positive number otherwise; used for sorting an array in ascending
 *  order
 */
static
SCIP_DECL_SORTINDCOMP(intComparatorAscending)
{  /*lint --e{715}*/
   int* tosort = (int*)dataptr;

   assert(tosort != NULL);

   return tosort[ind1] - tosort[ind2];
}

/** compares two entries of an int array, returns 0 if both are equal, a positive number if the entry at position ind2
 *  is larger than the one at position ind1, and a negative number otherwise; used for sorting an array in descending
 *  order
 */
static
SCIP_DECL_SORTINDCOMP(intComparatorDescending)
{  /*lint --e{715}*/
   int* tosort = (int*)dataptr;

   assert(tosort != NULL);

   return tosort[ind2] - tosort[ind1];
}

/** lexicographic comparison of two key pairs */
static
SCIP_DECL_SORTPTRCOMP(compareKeypairs)
{
   const SCIP_SORTKEYPAIR* first = (const SCIP_SORTKEYPAIR*)elem1;
   const SCIP_SORTKEYPAIR* second = (const SCIP_SORTKEYPAIR*)elem2;

   if( first->key1 < second->key1 )
      return -1;
   if( first->key1 > second->key1 )
      return 1;
   if( first->key2 < second->key2 )
      return -1;
   return first->key2 > second->key2;
}

/** GLOBAL VARIABLES **/
static SCIP* scip;
static int* tosort;
static int ntosort;

/* TEST SUITE */
static
void setup(void)
{
   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   ntosort = 8;
   SCIP_CALL( SCIPallocBufferArray(scip, &tosort, ntosort) );
   tosort[0] = 2;
   tosort[1] = 4;
   tosort[2] = 3;
   tosort[3] = 7;
   tosort[4] = 1;
   tosort[5] = 5;
   tosort[6] = 6;
   tosort[7] = 0;

   TEST_ASSERT_NOT_NULL(scip);
   TEST_ASSERT_NOT_NULL(tosort);
}

static
void teardown(void)
{

   SCIPfreeBufferArray(scip, &tosort);

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   /* check for memory leaks */
   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

/* Unity setUp/tearDown */
void setUp(void)
{
   setup();
}

void tearDown(void)
{
   teardown();
}

/* TESTS */

void test_sort_create_and_free(void)
{
   /* calls setup and teardown */
}

/** @brief tests SCIPsort by checking that the permutation is sorted ascending */
void test_sort_asc_sort_up(void)
{
   int* perm;
   int i;

   TEST_ASSERT_NOT_NULL(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &perm, ntosort) );

   SCIPsort(perm, intComparatorAscending, (void*)tosort, ntosort);

   for( i = 0; i < ntosort-1; i++ )
   {
      TEST_ASSERT_LESS_THAN(tosort[perm[i]], tosort[perm[i+1]]);
   }

   SCIPfreeBufferArray(scip, &perm);
}

/** @brief tests SCIPsortDown by checking that the permutation is sorted descending */
void test_sort_asc_sort_down(void)
{
   int* perm;
   int i;

   TEST_ASSERT_NOT_NULL(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &perm, ntosort) );

   SCIPsortDown(perm, intComparatorAscending, (void*)tosort, ntosort);

   for( i = 0; i < ntosort-1; i++ )
   {
      TEST_ASSERT_GREATER_THAN(tosort[perm[i]], tosort[perm[i+1]]);
   }

   SCIPfreeBufferArray(scip, &perm);
}

/** @brief tests SCIPsort by checking that the permutation is sorted descending */
void test_sort_desc_sort_up(void)
{
   int* perm;
   int i;

   TEST_ASSERT_NOT_NULL(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &perm, ntosort) );

   SCIPsort(perm, intComparatorDescending, (void*)tosort, ntosort);

   for( i = 0; i < ntosort-1; i++ )
   {
      TEST_ASSERT_GREATER_THAN(tosort[perm[i]], tosort[perm[i+1]]);
   }

   SCIPfreeBufferArray(scip, &perm);
}

/** @brief tests SCIPsortDown by checking that the permutation is sorted ascending */
void test_sort_desc_sort_down(void)
{
   int* perm;
   int i;

   TEST_ASSERT_NOT_NULL(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &perm, ntosort) );

   SCIPsortDown(perm, intComparatorDescending, (void*)tosort, ntosort);

   for( i = 0; i < ntosort-1; i++ )
   {
      TEST_ASSERT_LESS_THAN(tosort[perm[i]], tosort[perm[i+1]]);
   }

   SCIPfreeBufferArray(scip, &perm);
}

/** checks key-pair sorts against the generic API, including shell-sort/ninther boundaries and tied keys */
void test_sort_keypair_permutations(void)
{
   const int lengths[] = {0, 1, 2, 24, 25, 26, 728, 729, 730, 800};
   SCIP_SORTKEYPAIR data[800];
   SCIP_SORTKEYPAIR* localkeys[800];
   void* generickeys[800];
   int localidxs[800];
   int genericidxs[800];
   int n;
   int pattern;
   int down;
   int i;

   for( n = 0; n < (int)(sizeof(lengths) / sizeof(lengths[0])); ++n )
   {
      int len = lengths[n];

      for( pattern = 0; pattern < 6; ++pattern )
      {
         for( down = 0; down < 2; ++down )
         {
            for( i = 0; i < len; ++i )
            {
               /* signed-zero ties, repeated pairs, sorted, reverse sorted, mixed finite keys, and infinities */
               data[i].key1 = pattern == 0 ? (i % 2 == 0 ? 0.0 : -0.0)
                  : pattern == 1 ? (SCIP_Real)(i % 7) : pattern == 2 ? (SCIP_Real)i
                  : pattern == 3 ? (SCIP_Real)(len - i) : pattern == 4 ? (SCIP_Real)((i * 37) % 101 - 50)
                  : i % 3 == 0 ? -INFINITY : i % 3 == 1 ? 0.0 : INFINITY;
               data[i].key2 = pattern == 1 ? (SCIP_Real)((i / 7) % 3)
                  : pattern == 4 ? (SCIP_Real)((i * 53) % 29 - 14) * 0.125
                  : pattern == 5 ? (SCIP_Real)((i * 17) % 7) : 0.0;
               localkeys[i] = &data[i];
               generickeys[i] = &data[i];
               localidxs[i] = i;
               genericidxs[i] = i;
            }

            if( down )
            {
               SCIPsortDownPtrInt(generickeys, genericidxs, compareKeypairs, len);
               SCIPsortDownKeypairInt(localkeys, localidxs, len);
            }
            else
            {
               SCIPsortPtrInt(generickeys, genericidxs, compareKeypairs, len);
               SCIPsortKeypairInt(localkeys, localidxs, len);
            }

            for( i = 0; i < len; ++i )
            {
               TEST_ASSERT_EQUAL_INT(genericidxs[i], localidxs[i]);
               TEST_ASSERT_TRUE(generickeys[i] == localkeys[i]);
            }
         }
      }
   }
}

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_sort_create_and_free);
   RUN_TEST(test_sort_asc_sort_up);
   RUN_TEST(test_sort_asc_sort_down);
   RUN_TEST(test_sort_desc_sort_up);
   RUN_TEST(test_sort_desc_sort_down);
   RUN_TEST(test_sort_keypair_permutations);
   return UNITY_END();
}
