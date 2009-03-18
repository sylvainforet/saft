/* wafthash.c
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

#include "wafthash.h"


#define WAFT_HTABLE_SIZE 1024

static unsigned int waft_get_high_bit (unsigned int x);


WaftHNode*
waft_hnode_new ()
{
  WaftHNode *node;

  node        = malloc (sizeof (*node));
  node->seq   = NULL;
  node->next  = NULL;
  node->count = 0;

  return node;
}

void
waft_hnode_free (WaftHNode *node)
{
  if (node)
    {
      if (node->next)
        waft_hnode_free (node->next);
      free (node);
    }
}

WaftHTable*
waft_htable_new (WaftAlphabet *alphabet,
                 unsigned int  word_size)
{
  WaftHTable  *table;
  unsigned int i;

  table            = malloc (sizeof (*table));
  table->size      = WAFT_HTABLE_SIZE;
  table->shift     = waft_get_high_bit (alphabet->size);
  table->table     = malloc (sizeof (*table->table) * table->size);
  table->word_size = word_size;

  for (i = 0; i < table->size; i++)
    table->table[i] = NULL;

  return table;
}

void
waft_htable_clear (WaftHTable *table)
{
  unsigned int i;

  if (table)
    for (i = 0; i < table->size; i++)
      {
        waft_hnode_free (table->table[i]);
        table->table[i] = NULL;
      }
}

void
waft_htable_free (WaftHTable *table)
{
  unsigned int i;

  if (table)
    {
      for (i = 0; i < table->size; i++)
        waft_hnode_free (table->table[i]);
      free (table->table);
      free (table);
    }
}

void
waft_htable_add_seq (WaftHTable   *table,
                     WaftSequence *seq)
{
  unsigned int i;

  for (i = seq->size - table->word_size + 1; i-- > 0; )
    waft_htable_add (table, &seq->seq[i]);
}

void
waft_htable_add (WaftHTable *table,
                 WaftLetter *start)
{
  const unsigned int idx = waft_htable_hash (table, start);
  WaftHNode         *tmp;

  for (tmp = table->table[idx]; tmp; tmp = tmp->next)
    if (waft_htable_cmp (table, tmp, start))
      {
        tmp->count++;
        return;
      }

  tmp               = waft_hnode_new ();
  tmp->seq          = start;
  tmp->count        = 1;
  tmp->next         = table->table[idx];
  table->table[idx] = tmp;
}

unsigned int
waft_htable_hash (WaftHTable *table,
                   WaftLetter *start)
{
  unsigned int ret = *start;
  unsigned int i;

  for (i = 1; i < table->word_size; i++)
    ret = (ret << table->shift) | *++start;

  return ret % WAFT_HTABLE_SIZE;
}

int
waft_htable_cmp (WaftHTable *table,
                 WaftHNode  *node,
                 WaftLetter *start)
{
  unsigned int i;
  WaftLetter  *tmp = node->seq - 1;

  --start;
  for (i = 0; i < table->word_size; i++)
    if (*++start != *++tmp)
      return 0;

  return 1;
}

unsigned int
waft_htable_d2 (WaftHTable *tab1,
                WaftHTable *tab2)
{
  WaftHNode   *tmp1;
  WaftHNode   *tmp2;
  unsigned int i;
  unsigned int d2 = 0;

  for (i = 0; i < tab1->size; i++)
    for (tmp1 = tab1->table[i]; tmp1; tmp1 = tmp1->next)
      {
        tmp2 = waft_htable_lookup (tab2, tmp1->seq);
        if (tmp2)
          d2 += tmp1->count * tmp2->count;
      }

  return d2;
}

WaftHNode*
waft_htable_lookup (WaftHTable *table,
                    WaftLetter *start)
{
  const unsigned int idx = waft_htable_hash (table, start);
  WaftHNode         *tmp;

  for (tmp = table->table[idx]; tmp; tmp = tmp->next)
    if (waft_htable_cmp (table, tmp, start))
      return tmp;

  return NULL;
}

/**
 * From JJ (Joerg Arndt)
 */
static unsigned int
waft_get_high_bit (unsigned int x)
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

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
