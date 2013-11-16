/* saftsearchengines.c
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
#include <string.h>

#include "saftsearchengines.h"
#include "saftstats.h"


typedef enum
{
  NUC_A  = 0,
  NUC_C  = 1,
  NUC_G  = 2,
  NUC_T  = 3,
  NUC_NB = 4
}
Nucleotide;

static unsigned char char_to_bin_table[128] =
{
  ['a'] = NUC_A,
  ['c'] = NUC_C,
  ['g'] = NUC_G,
  ['t'] = NUC_T,
  ['A'] = NUC_A,
  ['C'] = NUC_C,
  ['G'] = NUC_G,
  ['T'] = NUC_T,
};


/*************************/
/* Generic Search Engine */
/*************************/

typedef struct _SearchEngineGeneric SearchEngineGeneric;

struct _SearchEngineGeneric
{
  SaftSearchEngine search_engine;
};


static void        search_engine_generic_free                 (SaftSearchEngine *engine);

static SaftSearch* search_engine_generic_search_two_sequences (SaftSearchEngine *engine,
                                                               SaftSequence     *query,
                                                               SaftSequence     *subject);

static SaftSearch* search_engine_generic_search_all           (SaftSearchEngine *engine,
                                                               const char       *query_path,
                                                               const char       *db_path);

SaftSearchEngine*
saft_search_engine_generic (SaftOptions *options)
{
  SearchEngineGeneric *engine;

  engine                                     = malloc (sizeof (*engine));
  engine->search_engine.options              = options;
  engine->search_engine.search_two_sequences = search_engine_generic_search_two_sequences;
  engine->search_engine.search_all           = search_engine_generic_search_all;
  engine->search_engine.free                 = search_engine_generic_free;

  return (SaftSearchEngine*)engine;
}

static void
search_engine_generic_free (SaftSearchEngine *engine)
{
}

static SaftSearch*
search_engine_generic_search_two_sequences (SaftSearchEngine *engine,
                                            SaftSequence     *query,
                                            SaftSequence     *subject)
{
  SearchEngineGeneric *se;
  
  se = (SearchEngineGeneric*) engine;
  se++;

  return NULL;
}

static SaftSearch*
search_engine_generic_search_all (SaftSearchEngine *engine,
                                  const char       *query_path,
                                  const char       *db_path)
{
  SearchEngineGeneric *se;
  
  se = (SearchEngineGeneric*) engine;
  se++;

  return NULL;
}

/*********************************/
/* Array-based DNA Search Engine */
/*********************************/

typedef struct _SearchEngineDNAArray SearchEngineDNAArray;

struct _SearchEngineDNAArray
{
  SaftSearchEngine  search_engine;

  SaftStatsContext *stats_context;

  WordCount       **query_cache;
  WordCount       **db_cache;
};


static void        search_engine_dna_array_free                 (SaftSearchEngine     *engine);

static SaftSearch* search_engine_dna_array_search_two_sequences (SaftSearchEngine     *engine,
                                                                 SaftSequence         *query,
                                                                 SaftSequence         *subject);

static SaftSearch* search_engine_dna_array_search_all           (SaftSearchEngine     *engine,
                                                                 const char           *query_path,
                                                                 const char           *db_path);

static WordCount*  search_engine_dna_array_hash_sequence        (SearchEngineDNAArray *engine,
                                                                 SaftSequence         *sequence);

SaftSearchEngine*
saft_search_engine_dna_array_new (SaftOptions *options)
{
  SearchEngineDNAArray *engine;

  engine                                     = malloc (sizeof (*engine));
  engine->search_engine.options              = options;
  engine->search_engine.search_two_sequences = search_engine_dna_array_search_two_sequences;
  engine->search_engine.search_all           = search_engine_dna_array_search_all;
  engine->search_engine.free                 = search_engine_dna_array_free;

  engine->stats_context                      = saft_stats_context_new (options->word_size,
                                                                       options->letter_frequencies, 4);
  engine->query_cache                        = NULL;
  engine->db_cache                           = NULL;

  return (SaftSearchEngine*)engine;
}

static void
search_engine_dna_array_free (SaftSearchEngine *engine)
{
  SearchEngineDNAArray *se;

  se = (SearchEngineDNAArray*) engine;

  saft_stats_context_free (se->stats_context);
  if (se->query_cache)
    free (se->query_cache);
  if (se->db_cache)
    free (se->db_cache);
}

static SaftSearch*
search_engine_dna_array_search_two_sequences (SaftSearchEngine *engine,
                                              SaftSequence     *query,
                                              SaftSequence     *subject)
{
  SearchEngineDNAArray *se;
  SaftSearch           *search;
  SaftResult           *res;
  WordCount            *counts_query;
  WordCount            *counts_subject;
  double                mean;
  double                var;
  int                   i;

  se = (SearchEngineDNAArray*) engine;

  counts_query   = search_engine_dna_array_hash_sequence (se, query);
  counts_subject = search_engine_dna_array_hash_sequence (se, subject);

  res = saft_result_new ();
  for (i = (1 << (2 * engine->options->word_size)) - 1; i >= 0; i--)
    res->d2 += counts_query[i] * counts_subject[i];

  mean = saft_stats_mean (se->stats_context,
                          query->seq_length,
                          subject->seq_length);
  var  = saft_stats_mean (se->stats_context,
                          query->seq_length,
                          subject->seq_length);

  res->name         = subject->name;
  res->p_value      = saft_stats_pgamma_m_v (res->d2, mean, var);
  res->p_value_adj  = res->p_value;

  search                 = saft_search_new ();
  search->query          = query;
  search->results        = res;
  search->sorted_results = &search->results;
  search->n_results      = 1;

  free (counts_query);
  free (counts_subject);

  return search;
}

static SaftSearch*
search_engine_dna_array_search_all (SaftSearchEngine *engine,
                                    const char       *query_path,
                                    const char       *db_path)
{
  SearchEngineDNAArray *se;

  se = (SearchEngineDNAArray*) engine;
  se++;

  return NULL;
}

static WordCount*
search_engine_dna_array_hash_sequence (SearchEngineDNAArray *engine,
                                       SaftSequence         *sequence)
{
  const size_t   k    = engine->search_engine.options->word_size;
  const uint16_t mask = 0xffff >> (16 - (2 * k));
  WordCount     *counts;
  size_t         i;
  uint16_t       w = 0;

  counts = calloc ((1 << (2 * k)), sizeof (*counts));
  if (sequence->seq_length < k)
    return counts;

  /* TODO Compare speed of array lookup and switch conditional */
  /* FIXME Check that the letters are in [ATGCatgc] */
  /* FIXME Handle periodic boundary conditions */
  for (i = 0; i < k; i++)
    {
      const unsigned char c = char_to_bin_table[i];

      w <<= 2;
      w |= c;
    }
  ++counts[w];
  for (i = k; i < sequence->seq_length; i++)
    {
      const unsigned char c = char_to_bin_table[i];

      w <<= 2;
      w |= c;
      w &= mask;
      ++counts[w];
    }
  return counts;
}

/********************************/
/* Hash-based DNA Search Engine */
/********************************/

typedef struct _SearchEngineDNAHash SearchEngineDNAHash;

struct _SearchEngineDNAHash
{
  SaftSearchEngine search_engine;
};


static void        search_engine_dna_hash_free                 (SaftSearchEngine *engine);

static SaftSearch* search_engine_dna_hash_search_two_sequences (SaftSearchEngine *engine,
                                                                SaftSequence     *query,
                                                                SaftSequence     *subject);

static SaftSearch* search_engine_dna_hash_search_all           (SaftSearchEngine *engine,
                                                                const char       *query_path,
                                                                const char       *db_path);

SaftSearchEngine*
saft_search_engine_dna_hash_new (SaftOptions *options)
{
  SearchEngineDNAHash *engine;

  engine                                     = malloc (sizeof (*engine));
  engine->search_engine.options              = options;
  engine->search_engine.search_two_sequences = search_engine_dna_hash_search_two_sequences;
  engine->search_engine.search_all           = search_engine_dna_hash_search_all;
  engine->search_engine.free                 = search_engine_dna_hash_free;

  return (SaftSearchEngine*)engine;
}

static void
search_engine_dna_hash_free (SaftSearchEngine *engine)
{
}

static SaftSearch*
search_engine_dna_hash_search_two_sequences (SaftSearchEngine *engine,
                                             SaftSequence     *query,
                                             SaftSequence     *subject)
{
  SearchEngineDNAHash *se;

  se = (SearchEngineDNAHash*) engine;
  se++;

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
