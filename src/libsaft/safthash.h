/* safthash.h
 * Copyright (C) 2013  Sylvain FORET
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *                                                                       
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *                                                                       
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Based on Glib's ghash table
 */

#ifndef __SAFT_KMERHASH_H__
#define __SAFT_KMERHASH_H__

#include "saftsequence.h"


#ifdef __cplusplus
extern "C"
{
#endif

/***************************/
/* Generic saft hash table */
/***************************/

/* Linear probing */
#define saft_hash_table_probe(step) (step)
/* For qudratic probing probing (linear seems marginally faster here) */
/* #define saft_hash_table_probe(step) ((step) * (step)) */


#define KMER_VAL_BYTES (sizeof (long))
#define KMER_VAL_NUCS  (KMER_VAL_BYTES * NUCS_PER_BYTE)

typedef union _SaftHashKmer SaftHashKmer;

union _SaftHashKmer
{
  unsigned char *kmer_ptr;
  unsigned char  kmer_vala[KMER_VAL_BYTES];
  unsigned long  kmer_vall;
};

typedef union _SaftHashValue SaftHashValue;

union _SaftHashValue
{
  long  count;
  void *ptr;
};

typedef struct _SaftHashNode SaftHashNode;

struct _SaftHashNode
{
  SaftHashKmer      kmer;

  /* If key_hash == 0, node is not in use
   * If key_hash == 1, node is a tombstone
   * If key_hash >= 2, node contains data
   */
  unsigned long key_hash;

  SaftHashValue value;
};

typedef unsigned long (*SaftHashFunc)  (const SaftHashKmer  *kmer,
                                        size_t               size);
typedef int           (*SaftEqualFunc) (const SaftHashKmer  *saft1,
                                        const SaftHashKmer  *saft2,
                                        size_t               size);

typedef struct _SaftHashTable  SaftHashTable;

struct _SaftHashTable
{
  SaftHashNode *nodes;

  long          size;
  long          mod;
  unsigned long mask;
  long          nnodes;
  long          noccupied;  /* nnodes + tombstones */

  unsigned int  k;
  unsigned int  kmer_bytes;

  SaftHashFunc  hash_func;
  SaftEqualFunc key_equal_func;
};

typedef struct _SaftHashTableIter SaftHashTableIter;

struct _SaftHashTableIter
{
  SaftHashTable *hash_table;
  long           position;
};

unsigned long  saft_hash_generic                   (const SaftHashKmer  *kmer,
                                                    size_t               size);

int            saft_equal_generic                  (const SaftHashKmer  *kmer1,
                                                    const SaftHashKmer  *kmer2,
                                                    size_t               size);

unsigned long  saft_hash_long                      (const SaftHashKmer  *kmer,
                                                    size_t               size);

int            saft_equal_long                     (const SaftHashKmer  *kmer1,
                                                    const SaftHashKmer  *kmer2,
                                                    size_t               size);

SaftHashTable* saft_hash_table_new                 (size_t               k);

SaftHashTable* saft_hash_table_new_full            (SaftHashFunc         hash_func,
                                                    SaftEqualFunc        key_equal_func,
                                                    size_t               k);

void           saft_hash_table_destroy             (SaftHashTable       *hash_table);

void           saft_hash_table_increment           (SaftHashTable       *hash_table,
                                                    const SaftHashKmer  *kmer);

void           saft_hash_table_add_count           (SaftHashTable       *hash_table,
                                                    const SaftHashKmer  *kmer,
                                                    long                 count);

SaftHashNode*  saft_hash_table_lookup_or_create    (SaftHashTable       *hash_table,
                                                    const SaftHashKmer  *kmer);

SaftHashNode*  saft_hash_table_lookup              (SaftHashTable       *hash_table,
                                                    const SaftHashKmer  *kmer);

void           saft_hash_table_iter_init           (SaftHashTableIter   *iter,
                                                    SaftHashTable       *hash_table);

SaftHashNode*  saft_hash_table_iter_next           (SaftHashTableIter   *iter);

#ifdef __cplusplus
}
#endif

#endif /* __SAFT_KMERHASH_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
