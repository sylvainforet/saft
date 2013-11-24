/* safthash.c
 * Copyright (C) 2013  Sylvain FORET
 *
 * This file is part of libngs.
 *
 * libngs is free software: you can redistribute it and/or modify
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

#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "safthash.h"
#include "saftutils.h"

#define HASH_TABLE_MIN_SHIFT 3  /* 1 << 3 == 8 buckets */

/* Each table size has an associated prime modulo (the first prime
 * lower than the table size) used to find the initial bucket. Probing
 * then works modulo 2^n. The prime modulo is necessary to get a
 * good distribution with poor hash functions. */
static const long prime_mod[] =
{
  1ul,          /* For 1 << 0 */
  2ul,
  3ul,
  7ul,
  13ul,
  31ul,
  61ul,
  127ul,
  251ul,
  509ul,
  1021ul,
  2039ul,
  4093ul,
  8191ul,
  16381ul,
  32749ul,
  65521ul,      /* For 1 << 16 */
  131071ul,
  262139ul,
  524287ul,
  1048573ul,
  2097143ul,
  4194301ul,
  8388593ul,
  16777213ul,
  33554393ul,
  67108859ul,
  134217689ul,
  268435399ul,
  536870909ul,
  1073741789ul,
  2147483647ul,  /* For 1 << 31 */
  4294967291ul,
  8589934583ul,
  17179869143ul,
  34359738337ul,
  68719476731ul,
  137438953447ul,
  274877906899ul,
  549755813881ul,
  1099511627689ul,
  2199023255531ul,
  4398046511093ul,
  8796093022151ul,
  17592186044399ul,
  35184372088777ul,
  70368744177643ul,
  140737488355213ul,
  281474976710597ul,
  562949953421231ul,
  1125899906842597ul,
  2251799813685119ul,
  4503599627370449ul,
  9007199254740881ul,
  18014398509481951ul,
  36028797018963913ul,
  72057594037927931ul,
  144115188075855859ul,
  288230376151711717ul,
  576460752303423433ul,
  1152921504606846883ul,
  2305843009213693951ul,
  4611686018427387847ul,
  9223372036854775783ul,
  18446744073709551557ul
};


static void                 saft_hash_table_set_shift                 (SaftHashTable       *hash_table,
                                                                       int                  shift);

static int                  saft_hash_table_find_closest_shift        (long n);

static void                 saft_hash_table_set_shift_from_size       (SaftHashTable       *hash_table,
                                                                       long                 size);

static inline void          saft_hash_table_copy_kmer                 (SaftHashTable       *hash_table,
                                                                       const unsigned char *kmer1,
                                                                       SaftHashKmer        *kmer2);

static inline unsigned long saft_hash_table_lookup_node_for_insertion (SaftHashTable       *hash_table,
                                                                       const unsigned char *kmer,
                                                                       unsigned long       *hash_return);


static void                 saft_hash_table_resize                    (SaftHashTable       *hash_table);

static inline int           saft_hash_table_maybe_resize              (SaftHashTable       *hash_table);


unsigned long
saft_hash_generic (const unsigned char *kmer,
                   size_t               size)
{
  unsigned long int hash = 0;
  unsigned int      i    = 0;

  do
    {
      hash <<= 8;
      hash  |= kmer[i];
    }
  while (++i < size);

  return ((long)(hash * 2654435769U));
}

int
saft_equal_generic (const unsigned char *kmer1,
                    SaftHashKmer        *kmer2,
                    size_t               size)
{
  unsigned int i;

  for (i = 0; i < size; i++)
    if (kmer1[i] != kmer2->kmer_ptr[i])
      return 0;

  return 1;
}

unsigned long
saft_hash_32bp (const unsigned char *kmer,
                size_t               size)
{
  unsigned long key = *((unsigned long*)kmer);

  key = (~key) + (key << 21);
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8);
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4);
  key = key ^ (key >> 28);
  key = key + (key << 31);

  return key;
}

int
saft_equal_32bp (const unsigned char *kmer1,
                 SaftHashKmer        *kmer2,
                 size_t               size)
{
  /* Completely unrolled */
  const int b1 = (kmer1[0] == kmer2->kmer_val[0]);
  const int b2 = (kmer1[1] == kmer2->kmer_val[1]);
  const int b3 = (kmer1[2] == kmer2->kmer_val[2]);
  const int b4 = (kmer1[3] == kmer2->kmer_val[3]);
  const int b5 = (kmer1[4] == kmer2->kmer_val[4]);
  const int b6 = (kmer1[5] == kmer2->kmer_val[5]);
  const int b7 = (kmer1[6] == kmer2->kmer_val[6]);
  const int b8 = (kmer1[7] == kmer2->kmer_val[7]);

  const int b12 = b1 && b2;
  const int b34 = b3 && b4;
  const int b56 = b5 && b6;
  const int b78 = b7 && b8;

  const int b1234 = b12 && b34;
  const int b5678 = b56 && b78;

  return b1234 && b5678;
}

static void
saft_hash_table_set_shift (SaftHashTable *hash_table,
                           int            shift)
{
  unsigned long mask = 0;
  int           i;

  hash_table->size = 1L << shift;
  hash_table->mod  = prime_mod[shift];

  for (i = 0; i < shift; i++)
    {
      mask <<= 1L;
      mask |= 1;
    }

  hash_table->mask = mask;
}

static int
saft_hash_table_find_closest_shift (long n)
{
  int i;

  for (i = 0; n; i++)
    n >>= 1L;

  return i;
}

static void
saft_hash_table_set_shift_from_size (SaftHashTable *hash_table,
                                     long           size)
{
  int shift;

  shift = saft_hash_table_find_closest_shift (size);

  if (shift < HASH_TABLE_MIN_SHIFT)
    shift = HASH_TABLE_MIN_SHIFT;

  saft_hash_table_set_shift (hash_table, shift);
}

static inline void
saft_hash_table_copy_kmer (SaftHashTable       *hash_table,
                           const unsigned char *kmer1,
                           SaftHashKmer        *kmer2)
{
  if (hash_table->kmer_bytes > KMER_VAL_BYTES)
    /* FIXME FIXME FIXME do we really want a terminating '\0' */
    kmer2->kmer_ptr = (unsigned char*)strndup ((const char*)kmer1, hash_table->kmer_bytes);
  else
    {
      int i;

      for (i = 4; i < hash_table->kmer_bytes; i += 4)
        {
          const int idx1 = i - 4;
          const int idx2 = i - 3;
          const int idx3 = i - 2;
          const int idx4 = i - 1;

          kmer2->kmer_val[idx1] = kmer1[idx1];
          kmer2->kmer_val[idx2] = kmer1[idx2];
          kmer2->kmer_val[idx3] = kmer1[idx3];
          kmer2->kmer_val[idx4] = kmer1[idx4];
        }
      for (i -= 4; i < hash_table->kmer_bytes; i++)
        kmer2->kmer_val[i] = kmer1[i];
    }
}

SaftHashNode*
saft_hash_table_lookup (SaftHashTable       *hash_table,
                        const unsigned char *kmer)
{
  unsigned long hash_value;

  /* Empty buckets have hash_value set to 0, and for tombstones, it's 1.
   * We need to make sure our hash value is not one of these.
   */
  hash_value = (* hash_table->hash_func) (kmer, hash_table->kmer_bytes);
  if (SAFT_UNLIKELY (hash_value <= 1))
    hash_value = 2;

  return saft_hash_table_lookup_with_key (hash_table, kmer, hash_value);
}

SaftHashNode*
saft_hash_table_lookup_with_key (SaftHashTable       *hash_table,
                                 const unsigned char *kmer,
                                 unsigned long        key)
{
  SaftHashNode       *node;
  unsigned long       node_index;
  const unsigned long hash_value = key;
  unsigned long       step       = 0;

  /* Empty buckets have hash_value set to 0, and for tombstones, it's 1.
   * We need to make sure our hash value is not one of these.
   */
  node_index = hash_value % hash_table->mod;
  node       = &hash_table->nodes[node_index];

  while (node->key_hash)
    {
      /* We first check if our full hash values
       * are equal so we can avoid calling the full-blown
       * key equality function in most cases.
       */
      if (node->key_hash == hash_value)
        if (hash_table->key_equal_func (kmer, &node->kmer, hash_table->kmer_bytes))
          break;
      step++;
      /* Linear probing
       * TODO compare this with quadratic probing:
       * node_index += (step * step);
       * */
      node_index += step;
      node_index &= hash_table->mask;
      node        = &hash_table->nodes[node_index];
    }
  return node->key_hash ? node : NULL;
}

static inline unsigned long
saft_hash_table_lookup_node_for_insertion (SaftHashTable       *hash_table,
                                           const unsigned char *kmer,
                                           unsigned long       *hash_return)
{
  SaftHashNode *node;
  unsigned long node_index;
  unsigned long hash_value;
  unsigned long first_tombstone = 0;
  unsigned long step            = 0;
  int           have_tombstone  = 0;

  /* Empty buckets have hash_value set to 0, and for tombstones, it's 1.
   * We need to make sure our hash value is not one of these.
   */
  hash_value = (* hash_table->hash_func) (kmer, hash_table->kmer_bytes);
  if (SAFT_UNLIKELY (hash_value <= 1))
    hash_value = 2;

  *hash_return = hash_value;
  node_index   = hash_value % hash_table->mod;
  node         = &hash_table->nodes[node_index];

  while (node->key_hash)
    {
      /* We first check if our full hash values
       * are equal so we can avoid calling the full-blown
       * key equality function in most cases.
       */
      if (node->key_hash == hash_value)
        {
          if (hash_table->key_equal_func (kmer, &node->kmer, hash_table->kmer_bytes))
            return node_index;
        }
      else if (node->key_hash == 1 && !have_tombstone)
        {
          first_tombstone = node_index;
          have_tombstone  = 1;
        }
      step++;
      node_index += step;
      node_index &= hash_table->mask;
      node        = &hash_table->nodes[node_index];
    }
  if (have_tombstone)
    return first_tombstone;

  return node_index;
}

static void
saft_hash_table_resize (SaftHashTable *hash_table)
{
  SaftHashNode *new_nodes;
  long old_size;
  long i;

  old_size = hash_table->size;
  saft_hash_table_set_shift_from_size (hash_table, hash_table->nnodes * 2);

  new_nodes = calloc (hash_table->size, sizeof (*new_nodes));

  for (i = 0; i < old_size; i++)
    {
      SaftHashNode  *node = &hash_table->nodes[i];
      SaftHashNode  *new_node;
      unsigned long  hash_val;
      unsigned long  step = 0;

      if (node->key_hash <= 1)
        continue;

      hash_val = node->key_hash % hash_table->mod;
      new_node = &new_nodes[hash_val];

      while (new_node->key_hash)
        {
          step++;
          hash_val += step;
          hash_val &= hash_table->mask;
          new_node  = &new_nodes[hash_val];
        }
      new_node->kmer     = node->kmer;
      new_node->key_hash = node->key_hash;
      new_node->value    = node->value;
    }

  free (hash_table->nodes);
  hash_table->nodes = new_nodes;
  hash_table->noccupied = hash_table->nnodes;

#if 0
/* TODO check how to tombstones are handled */

#define IS_TOUCHED(i) ((touched[(i) / 8] >> ((i) & 7)) & 1)
#define SET_TOUCHED(i) (touched[(i) / 8] |= (1 << ((i) & 7)))

  const unsigned long old_size = hash_table->size;
  unsigned char      *touched;
  long                i;

  saft_hash_table_set_shift_from_size (hash_table, hash_table->nnodes * 2);
  hash_table->nodes = g_realloc (hash_table->nodes,
                                 hash_table->size * sizeof (*hash_table->nodes));
  memset (hash_table->nodes + old_size,
          0,
          (hash_table->size - old_size) * sizeof (*hash_table->nodes));
  touched = g_malloc0 ((old_size / 8 + 1) * sizeof (*touched));

  for (i = 0; i < old_size; i++)
    {
      SaftHashNode  *node;
      SaftHashNode  *next_node;
      SaftHashNode  tmp_node;
      unsigned long hash_val;
      long          step;

      if (IS_TOUCHED (i))
        continue;
      node = &hash_table->nodes[i];
      if (node->key_hash == 0)
        continue;
      hash_val = node->key_hash % hash_table->mod;
      if (hash_val == i)
        {
          SET_TOUCHED (i);
          continue;
        }
      tmp_node = *node;
      next_node = &hash_table->nodes[hash_val];
      node->key_hash = 0;
      step = 0;
      while (next_node->key_hash)
        {
          /* Hit an untouched element */
          if (hash_val < old_size && !IS_TOUCHED (hash_val))
            {
              SaftHashNode tmpp;

              tmpp = tmp_node;
              tmp_node = *next_node;
              *next_node = tmpp;
              SET_TOUCHED (hash_val);
              hash_val = tmp_node.key_hash % hash_table->mod;
              step = 0;
            }
          else
            {
              step++;
              hash_val += step;
              hash_val &= hash_table->mask;
            }
          next_node = &hash_table->nodes[hash_val];
        }
      *next_node = tmp_node;
    }
  g_free (touched);
#undef IS_TOUCHED
#undef SET_TOUCHED
#endif
}

static inline int
saft_hash_table_maybe_resize (SaftHashTable *hash_table)
{
  long noccupied = hash_table->noccupied;
  long size = hash_table->size;

  if ((size > hash_table->nnodes * 4 && size > 1 << HASH_TABLE_MIN_SHIFT) ||
      (size <= noccupied + (noccupied / 16)))
    {
      saft_hash_table_resize (hash_table);
      return 1;
    }
  return 0;
}

SaftHashTable*
saft_hash_table_new (size_t k)
{
  SaftHashTable *hash_table;

  if (k <= KMER_VAL_NUCS)
    hash_table = saft_hash_table_new_full (saft_hash_32bp,
                                           saft_equal_32bp,
                                           k);
  else
    hash_table = saft_hash_table_new_full (saft_hash_generic,
                                           saft_equal_generic,
                                           k);

  return hash_table;
}

SaftHashTable*
saft_hash_table_new_full (SaftHashFunc  hash_func,
                          SaftEqualFunc key_equal_func,
                          size_t        k)
{
  SaftHashTable *hash_table;

  hash_table                     = malloc (sizeof (*hash_table));
  saft_hash_table_set_shift (hash_table, HASH_TABLE_MIN_SHIFT);
  hash_table->nnodes             = 0;
  hash_table->noccupied          = 0;
  hash_table->hash_func          = hash_func;
  hash_table->key_equal_func     = key_equal_func;
  hash_table->nodes              = calloc (hash_table->size, sizeof (*hash_table->nodes));
  hash_table->k                  = k;
  hash_table->kmer_bytes         = (k + NUCS_PER_BYTE - 1) / NUCS_PER_BYTE;

  return hash_table;
}

void
saft_hash_table_iter_init (SaftHashTableIter *iter,
                           SaftHashTable     *hash_table)
{
  iter->hash_table = hash_table;
  iter->position   = -1;
}

SaftHashNode*
saft_hash_table_iter_next (SaftHashTableIter *iter)
{
  SaftHashNode *node      = NULL;
  long          position  = iter->position;

  do
    {
      position++;
      if (position >= iter->hash_table->size)
        {
          iter->position = position;
          return NULL;
        }
      node = &iter->hash_table->nodes[position];
    }
  while (node->key_hash <= 1);
  iter->position = position;

  return node;
}

void
saft_hash_table_destroy (SaftHashTable *hash_table)
{
  if (hash_table)
    {
      /* FIXME FIXME FIXME free the copied kmers if any */
      free (hash_table);
    }
}

void
saft_hash_table_increment (SaftHashTable       *hash_table,
                           const unsigned char *kmer)
{
  saft_hash_table_add_count (hash_table, kmer, 1);
}

void
saft_hash_table_add_count (SaftHashTable       *hash_table,
                           const unsigned char *kmer,
                           long                 count)
{
  SaftHashNode *node;
  unsigned long node_index;
  unsigned long key_hash;
  unsigned long old_hash;

  node_index = saft_hash_table_lookup_node_for_insertion (hash_table, kmer, &key_hash);
  node       = &hash_table->nodes[node_index];
  old_hash   = node->key_hash;

  if (old_hash > 1)
    node->value.count += count;
  else
    {
      saft_hash_table_copy_kmer (hash_table, kmer, &node->kmer);
      node->value.count = count;
      node->key_hash    = key_hash;
      hash_table->nnodes++;
      if (old_hash == 0)
        {
          /* We replaced an empty node, and not a tombstone */
          hash_table->noccupied++;
          saft_hash_table_maybe_resize (hash_table);
        }
    }
}

SaftHashNode*
saft_hash_table_lookup_or_create (SaftHashTable       *hash_table,
                                  const unsigned char *kmer)
{
  SaftHashNode *node;
  unsigned long node_index;
  unsigned long key_hash;
  unsigned long old_hash;

  node_index = saft_hash_table_lookup_node_for_insertion (hash_table, kmer, &key_hash);
  node       = &hash_table->nodes[node_index];
  old_hash   = node->key_hash;

  if (old_hash <= 1)
    {
      saft_hash_table_copy_kmer (hash_table, kmer, &node->kmer);
      node->key_hash    = key_hash;
      node->value.count = 0;
      node->value.ptr   = NULL;
      hash_table->nnodes++;
      if (old_hash == 0)
        {
          /* We replaced an empty node, and not a tombstone */
          hash_table->noccupied++;
          /* Addresses changed */
          if (saft_hash_table_maybe_resize (hash_table))
            {
              node_index = saft_hash_table_lookup_node_for_insertion (hash_table, kmer, &key_hash);
              node       = &hash_table->nodes[node_index];
            }
        }
    }
  return node;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
