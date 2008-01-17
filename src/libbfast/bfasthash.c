/* bfasthash.c
 * Copyright (C) 2008  Sylvain FORET
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
 *
 */

#include <stdlib.h>

#include "bfasthash.h"


#define BFAST_HTABLE_SIZE 1024

static unsigned int bfast_get_high_bit (unsigned int x);


BfastHNode*
bfast_hnode_new ()
{
  BfastHNode *node;

  node        = malloc (sizeof (*node));
  node->seq   = NULL;
  node->next  = NULL;
  node->count = 0;

  return node;
}

void
bfast_hnode_free (BfastHNode *node)
{
  if (node)
    {
      if (node->next)
        bfast_hnode_free (node->next);
      free (node);
    }
}

BfastHTable*
bfast_htable_new (BfastAlphabet *alphabet,
                  unsigned int   word_size)
{
  BfastHTable *table;
  unsigned int i;

  table            = malloc (sizeof (*table));
  table->size      = BFAST_HTABLE_SIZE;
  table->shift     = bfast_get_high_bit (alphabet->size);
  table->table     = malloc (sizeof (*table->table) * table->size);
  table->word_size = word_size;

  for (i = 0; i < table->size; i++)
    table->table[i] = NULL;

  return table;
}

void
bfast_htable_clear (BfastHTable *table)
{
  unsigned int i;

  if (table)
    for (i = 0; i < table->size; i++)
      {
        bfast_hnode_free (table->table[i]);
        table->table[i] = NULL;
      }
}

void
bfast_htable_free (BfastHTable *table)
{
  unsigned int i;

  if (table)
    {
      for (i = 0; i < table->size; i++)
        bfast_hnode_free (table->table[i]);
      free (table->table);
      free (table);
    }
}

void
bfast_htable_add_seq (BfastHTable   *table,
                      BfastSequence *seq)
{
  unsigned int i;

  for (i = seq->size - table->word_size + 1; i-- > 0; )
    bfast_htable_add (table, &seq->seq[i]);
}

void
bfast_htable_add (BfastHTable *table,
                  BfastLetter *start)
{
  const unsigned int idx = bfast_htable_hash (table, start);
  BfastHNode        *tmp;

  for (tmp = table->table[idx]; tmp; tmp = tmp->next)
    if (bfast_htable_cmp (table, tmp, start))
      {
        tmp->count++;
        return;
      }

  tmp               = bfast_hnode_new ();
  tmp->seq          = start;
  tmp->count        = 1;
  tmp->next         = table->table[idx];
  table->table[idx] = tmp;
}

unsigned int
bfast_htable_hash (BfastHTable *table,
                   BfastLetter *start)
{
  unsigned int ret = *start;
  unsigned int i;

  for (i = 1; i < table->word_size; i++)
    ret = (ret << table->shift) | *++start;

  return ret % BFAST_HTABLE_SIZE;
}

int
bfast_htable_cmp (BfastHTable *table,
                  BfastHNode  *node,
                  BfastLetter *start)
{
  unsigned int i;
  BfastLetter *tmp = node->seq - 1;

  --start;
  for (i = 0; i < table->word_size; i++)
    if (*++start != *++tmp)
      return 0;

  return 1;
}

unsigned int
bfast_htable_d2 (BfastHTable *tab1,
                 BfastHTable *tab2)
{
  BfastHNode   *tmp1;
  BfastHNode   *tmp2;
  unsigned int  i;
  unsigned int  d2 = 0;

  for (i = 0; i < tab1->size; i++)
    for (tmp1 = tab1->table[i]; tmp1; tmp1 = tmp1->next)
      {
        tmp2 = bfast_htable_lookup (tab2, tmp1->seq);
        if (tmp2)
          d2 += tmp1->count * tmp2->count;
      }

  return d2;
}

BfastHNode*
bfast_htable_lookup (BfastHTable *table,
                     BfastLetter *start)
{
  const unsigned int idx = bfast_htable_hash (table, start);
  BfastHNode        *tmp;

  for (tmp = table->table[idx]; tmp; tmp = tmp->next)
    if (bfast_htable_cmp (table, tmp, start))
      return tmp;

  return NULL;
}

/**
 * From JJ (Joerg Arndt)
 */
static unsigned int
bfast_get_high_bit (unsigned int x)
{
  unsigned int r = 0;

  if (x & 0xffff0000)
    {
      x >>= 16;
      r  += 16;
    }
  if (x & 0x0000ff00)
    {
      x >>= 8;
      r  += 8;
    }
  if (x & 0x000000f0)
    {
      x >>= 4;
      r  += 4;
    }
  if (x & 0x0000000c)
    {
      x >>= 2;
      r  += 2;
    }
  if (x & 0x00000002)
    r += 1;

  return r;
}
