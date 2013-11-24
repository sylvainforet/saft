/* saftsearchenginednahash.c
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


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "saftfasta.h"
#include "safthash.h"
#include "saftsearchengines.h"
#include "saftstats.h"


/********************************/
/* Hash-based DNA Search Engine */
/********************************/

typedef struct _SearchEngineDNAHash SearchEngineDNAHash;

struct _SearchEngineDNAHash
{
  SaftSearchEngine  search_engine;

  SaftStatsContext *stats_context;
};


static void           search_engine_dna_hash_free                 (SaftSearchEngine     *engine);

static SaftHashTable* search_engine_dna_hash_hash_sequence        (SearchEngineDNAHash  *engine,
                                                                   SaftSequence         *sequence);

static unsigned long search_engine_dna_hash_d2                    (SearchEngineDNAHash  *engine,
                                                                   SaftHashTable        *counts1,
                                                                   SaftHashTable        *counts2);

static SaftSearch*    search_engine_dna_hash_search_two_sequences (SaftSearchEngine     *engine,
                                                                   SaftSequence         *query,
                                                                   SaftSequence         *subject);

static SaftSearch*    search_engine_dna_hash_search_all           (SaftSearchEngine     *engine,
                                                                   const char           *query_path,
                                                                   const char           *db_path);

SaftSearchEngine*
saft_search_engine_dna_hash_new (SaftOptions *options)
{
  SearchEngineDNAHash *engine;

  engine                                     = malloc (sizeof (*engine));
  engine->search_engine.options              = options;
  engine->search_engine.search_two_sequences = search_engine_dna_hash_search_two_sequences;
  engine->search_engine.search_all           = search_engine_dna_hash_search_all;
  engine->search_engine.free                 = search_engine_dna_hash_free;

  engine->stats_context                      = saft_stats_context_new (options->word_size,
                                                                       options->letter_frequencies, 4);

  return (SaftSearchEngine*)engine;
}

static void
search_engine_dna_hash_free (SaftSearchEngine *engine)
{
  SearchEngineDNAHash *se;

  se = (SearchEngineDNAHash*) engine;

  if (se->stats_context)
    saft_stats_context_free (se->stats_context);
}

static SaftHashTable*
search_engine_dna_hash_hash_sequence (SearchEngineDNAHash  *engine,
                                      SaftSequence         *sequence)
{
  SaftHashTable *table;
  const size_t   k = engine->search_engine.options->word_size;
  size_t         i;

  table = saft_hash_table_new (engine->search_engine.options->word_size);

  if (k <= KMER_VAL_NUCS)
    {
      const unsigned long mask = (~ 0L) >> (8 * sizeof (unsigned long) - (2 * k));
      unsigned long       w    = 0;

      for (i = 0; i < k; i++)
        {
          const unsigned char c = SaftAlphabetDNA.codes[(int)sequence->seq[i]];

          w <<= 2;
          w |= c;
        }
      saft_hash_table_increment (table, (unsigned char*) &w);
      for (i = k; i < sequence->seq_length; i++)
        {
          const unsigned char c = SaftAlphabetDNA.codes[(int)sequence->seq[i]];

          w <<= 2;
          w |= c;
          w &= mask;
          saft_hash_table_increment (table, (unsigned char*) &w);
        }
    }
  else
    {
      /* FIXME implement this */
    }

  return table;
}

static unsigned long
search_engine_dna_hash_d2 (SearchEngineDNAHash  *engine,
                           SaftHashTable        *counts1,
                           SaftHashTable        *counts2)
{
  SaftHashTableIter    iter;
  SaftHashNode        *node;
  unsigned long        d2 = 0;

  saft_hash_table_iter_init (&iter, counts1);

  /* FIXME iterate over the smallest table would probably be better */
  for (node = saft_hash_table_iter_next (&iter);
       node;
       node = saft_hash_table_iter_next (&iter))
    {
      SaftHashNode *res;

      res = saft_hash_table_lookup_with_key (counts2,
                                             node->kmer.kmer_ptr,
                                             node->key_hash);
      if (res)
        d2 += node->value.count * res->value.count;
    }

  return d2;
}

static SaftSearch*
search_engine_dna_hash_search_two_sequences (SaftSearchEngine *engine,
                                             SaftSequence     *query,
                                             SaftSequence     *subject)
{
  SearchEngineDNAHash *se;
  SaftSearch          *search;
  SaftResult          *result;
  SaftHashTable       *hash_query;
  SaftHashTable       *hash_subject;
  double               mean;
  double               var;

  se                   = (SearchEngineDNAHash*) engine;
  hash_query           = search_engine_dna_hash_hash_sequence (se, query);
  hash_subject         = search_engine_dna_hash_hash_sequence (se, query);
  result               = saft_result_new ();
  result->d2           = search_engine_dna_hash_d2 (se, hash_query, hash_subject);
  mean                 = saft_stats_mean (se->stats_context,
                                          query->seq_length,
                                          subject->seq_length);
  var                  = saft_stats_mean (se->stats_context,
                                          query->seq_length,
                                          subject->seq_length);
  result->name         = strdup(subject->name);
  result->p_value      = saft_stats_pgamma_m_v (result->d2, mean, var);
  result->p_value_adj  = result->p_value;
  search               = saft_search_new (1);
  search->name         = strdup(query->name);
  saft_search_add_result (search, result);

  return NULL;
}

static SaftSearch*
search_engine_dna_hash_search_all (SaftSearchEngine *engine,
                                   const char       *query_path,
                                   const char       *db_path)
{
  SearchEngineDNAHash *se;

  se = (SearchEngineDNAHash*) engine;
  se++;

  return NULL;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
