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
#include <stdint.h>
#include <string.h>

#include "wafthash.h"


#define WAFT_HTABLE_SIZE 1024

static unsigned int waft_get_high_bit (unsigned int x);


WaftHNode*
waft_hnode_new ()
{
  WaftHNode *node;

  node                = malloc (sizeof (*node));
  node->seq           = NULL;
  node->next          = NULL;
  node->count_query   = 0;
  node->count_subject = 0;

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
  WaftHTable *table;

  table            = malloc (sizeof (*table));
  table->size      = WAFT_HTABLE_SIZE;
  table->shift     = waft_get_high_bit (alphabet->size) + 1;
  table->table     = malloc (sizeof (*table->table) * table->size);
  table->word_size = word_size;
  table->hmask     = ~0;

  if (table->word_size * table->shift < sizeof (table->hmask) * 8)
    table->hmask >>= sizeof (table->hmask) * 8 - table->word_size * table->shift;

  memset (table->table, 0, sizeof(*table->table));

  return table;
}

void
waft_htable_clear (WaftHTable *table)
{
  unsigned int i;

  if (table)
    {
      for (i = 0; i < table->size; i++)
        waft_hnode_free (table->table[i]);
      memset (table->table, 0, sizeof(*table->table));
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
waft_htable_add_query (WaftHTable   *table,
                       WaftSequence *seq)
{
  WaftLetter  *start;
  unsigned int i;
  unsigned int hash = 0;

  if (table->word_size > seq->size)
    return;

  start = seq->seq - 1;
  for (i = 1; i < table->word_size; i++)
    hash = (hash << table->shift) | *++start;

  start = seq->seq;
  i     = table->word_size - 1;
  do
    {
      WaftHNode *node;

      hash <<= table->shift;
      hash  |= seq->seq[i];
      hash  &= table->hmask;

      for (node = table->table[hash]; node; node = node->next)
        if (waft_htable_cmp (table, node, start))
          {
            node->count_query++;
            goto found;
          }
      node               = waft_hnode_new ();
      node->seq          = start;
      node->count_query  = 1;
      node->next         = table->table[hash];
      table->table[hash] = node;
found:
      ++i;
      ++start;
    }
  while (i < seq->size);
}

void
waft_htable_add_subject (WaftHTable   *table,
                         WaftSequence *seq)
{
  WaftLetter  *start;
  unsigned int i;
  unsigned int hash       = 0;

  if (table->word_size > seq->size)
    return;

  start = seq->seq - 1;
  for (i = 1; i < table->word_size; i++)
    hash = (hash << table->shift) | *++start;

  start = seq->seq;
  i     = table->word_size - 1;
  do
    {
      WaftHNode *node;

      hash <<= table->shift;
      hash  |= seq->seq[i];
      hash  &= table->hmask;

      for (node = table->table[hash]; node; node = node->next)
        if (waft_htable_cmp (table, node, start))
          {
            node->count_subject++;
            break;
          }
      ++i;
      ++start;
    }
  while (i < seq->size);
}

void
waft_htable_clear_subject (WaftHTable *table)
{
  unsigned int i;

  for (i = 0; i < table->size; i++)
    {
      WaftHNode *node;
      for (node = table->table[i]; node; node = node->next)
        node->count_subject = 0;
    }
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
waft_htable_d2 (WaftHTable *tab)
{
  WaftHNode   *node;
  unsigned int i;
  unsigned int d2 = 0;

  for (i = 0; i < tab->size; i++)
    for (node = tab->table[i]; node; node = node->next)
      d2 += node->count_query * node->count_subject;

  return d2;
}

WaftHNode*
waft_htable_lookup (WaftHTable *table,
                    WaftLetter *start)
{
  const unsigned int idx = waft_htable_hash (table, start);
  WaftHNode         *node;

  for (node = table->table[idx]; node; node = node->next)
    if (waft_htable_cmp (table, node, start))
      return node;

  return NULL;
}

/**
 * From JJ (Joerg Arndt) www.jjj.de
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
