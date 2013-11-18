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


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "saftfasta.h"
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

typedef struct _DNAArrayDBEntry DNAArrayDBEntry;

struct _DNAArrayDBEntry
{
  DNAArrayDBEntry *next;
  WordCount       *counts;
  char            *name;
  size_t           length;
};

static DNAArrayDBEntry* dna_array_db_entry_new      (void);

static void             dna_array_db_entry_free     (DNAArrayDBEntry *entry);

static void             dna_array_db_entry_free_all (DNAArrayDBEntry *entry);


typedef struct _SearchEngineDNAArray SearchEngineDNAArray;

struct _SearchEngineDNAArray
{
  SaftSearchEngine  search_engine;

  SaftStatsContext *stats_context;

  DNAArrayDBEntry  *query_cache;
  DNAArrayDBEntry  *db_cache;

  SaftSearch       *search;
  SaftSearch      **search_array;
  SaftSearch       *tmp_search;
  WordCount        *tmp_counts;

  size_t            n_queries;
  size_t            max_words;
  size_t            tmp_length;
};


static void        search_engine_dna_array_free                 (SaftSearchEngine     *engine);

static WordCount*  search_engine_dna_array_hash_sequence        (SearchEngineDNAArray *engine,
                                                                 SaftSequence         *sequence);

static SaftSearch* search_engine_dna_array_search_two_sequences (SaftSearchEngine     *engine,
                                                                 SaftSequence         *query,
                                                                 SaftSequence         *subject);

static SaftSearch* search_engine_dna_array_search_all           (SaftSearchEngine     *engine,
                                                                 const char           *query_path,
                                                                 const char           *db_path);

static SaftSearch* search_engine_dna_array_search_all_no_cache  (SearchEngineDNAArray *engine,
                                                                 const char           *query_path,
                                                                 const char           *db_path);

static SaftSearch* search_engine_dna_array_search_all_qcached   (SearchEngineDNAArray *engine,
                                                                 const char           *query_path,
                                                                 const char           *db_path);

static SaftSearch* search_engine_dna_array_search_all_dcached   (SearchEngineDNAArray *engine,
                                                                 const char           *query_path,
                                                                 const char           *db_path);

static int         search_engine_dna_array_cache_sequence       (SaftSequence         *sequence,
                                                                 void                 *data);

static int         search_engine_dna_array_search_query         (SaftSequence         *sequence,
                                                                 void                 *data);

static int         search_engine_dna_array_search_db            (SaftSequence         *sequence,
                                                                 void                 *data);

static int         search_engine_dna_array_queries_iter_func    (SaftSequence         *sequence,
                                                                 void                 *data);

static int         search_engine_dna_array_db_iter_func         (SaftSequence         *sequence,
                                                                 void                 *data);



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
  engine->search                             = NULL;
  engine->search_array                       = NULL;
  engine->tmp_counts                         = NULL;
  engine->n_queries                          = 0;
  engine->max_words                          = 1 << (2 * options->word_size);
  engine->tmp_length                         = 0;

  return (SaftSearchEngine*)engine;
}

static void
search_engine_dna_array_free (SaftSearchEngine *engine)
{
  SearchEngineDNAArray *se;

  se = (SearchEngineDNAArray*) engine;

  if (se->stats_context)
    saft_stats_context_free (se->stats_context);
  if (se->query_cache)
    dna_array_db_entry_free_all (se->query_cache);
  if (se->db_cache)
    dna_array_db_entry_free_all (se->db_cache);
  if (se->search)
    saft_search_free (se->search);
  if (se->search_array)
    free (se->search_array);

  free (se);
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

  counts = calloc (engine->max_words, sizeof (*counts));
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

static SaftSearch*
search_engine_dna_array_search_two_sequences (SaftSearchEngine *engine,
                                              SaftSequence     *query,
                                              SaftSequence     *subject)
{
  SearchEngineDNAArray *se;
  SaftSearch           *search;
  SaftResult           *result;
  WordCount            *counts_query;
  WordCount            *counts_subject;
  double                mean;
  double                var;
  int                   i;

  se = (SearchEngineDNAArray*) engine;

  counts_query   = search_engine_dna_array_hash_sequence (se, query);
  counts_subject = search_engine_dna_array_hash_sequence (se, subject);

  result = saft_result_new ();
  for (i = se->max_words - 1; i >= 0; i--)
    result->d2 += counts_query[i] * counts_subject[i];

  mean = saft_stats_mean (se->stats_context,
                          query->seq_length,
                          subject->seq_length);
  var  = saft_stats_mean (se->stats_context,
                          query->seq_length,
                          subject->seq_length);

  result->name         = strdup(subject->name);
  result->p_value      = saft_stats_pgamma_m_v (result->d2, mean, var);
  result->p_value_adj  = result->p_value;

  search                 = saft_search_new (1);
  search->name           = strdup(query->name);
  saft_search_add_result (search, result);

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
  SaftSearch           *ret = NULL;


  se = (SearchEngineDNAArray*) engine;

  if (engine->options->cache_db)
    ret = search_engine_dna_array_search_all_dcached (se, query_path, db_path);
  else if (engine->options->cache_queries)
    ret = search_engine_dna_array_search_all_qcached (se, query_path, db_path);
  else
    ret = search_engine_dna_array_search_all_no_cache (se, query_path, db_path);

  return ret;
}

static SaftSearch*
search_engine_dna_array_search_all_no_cache (SearchEngineDNAArray *engine,
                                             const char           *query_path,
                                             const char           *db_path)
{
  SaftSearch *ret;

  saft_fasta_iter (query_path,
                   search_engine_dna_array_queries_iter_func,
                   engine);

  ret = engine->search;
  engine->search = NULL;

  return ret;
}

static int
search_engine_dna_array_queries_iter_func (SaftSequence *sequence,
                                           void         *data)
{
  SearchEngineDNAArray *engine;

  engine = (SearchEngineDNAArray*)data;

  engine->tmp_counts       = search_engine_dna_array_hash_sequence (engine, sequence);
  engine->tmp_length       = sequence->seq_length;
  engine->tmp_search       = saft_search_new (engine->search_engine.options->max_results);
  engine->tmp_search->name = strdup (sequence->name);

  saft_fasta_iter (engine->search_engine.options->db_path,
                   search_engine_dna_array_db_iter_func,
                   engine);

  free (engine->tmp_counts);
  engine->tmp_counts = NULL;

  if (engine->tmp_search->n_results == 0)
    saft_search_free (engine->tmp_search);
  else
    {
      engine->tmp_search->next = engine->search;
      engine->search           = engine->tmp_search;
      engine->tmp_search       = NULL;
    }

  return 1;
}

static int
search_engine_dna_array_db_iter_func (SaftSequence *sequence,
                                      void         *data)
{
  SearchEngineDNAArray *engine;
  WordCount            *counts;
  unsigned long         d2 = 0;
  double                mean;
  double                var;
  int                   i;

  engine = (SearchEngineDNAArray*)data;
  counts = search_engine_dna_array_hash_sequence (engine, sequence);

  for (i = 0; i < engine->max_words; i++)
    d2 += engine->tmp_counts[i] * counts[i];

  free (counts);

  mean = saft_stats_mean (engine->stats_context,
                          sequence->seq_length,
                          engine->tmp_length);
  var  = saft_stats_mean (engine->stats_context,
                          sequence->seq_length,
                          engine->tmp_length);

  /* FIXME adjust this euristic depending on the user's required significance level */
  /* Fix this here and everywhere else in this file */
  if (d2 >= mean + 2 * sqrt (var))
    {
      SaftResult    *result;

      result          = saft_result_new ();
      result->d2      = d2;
      result->p_value = saft_stats_pgamma_m_v (result->d2, mean, var);
      result->name    = strdup(sequence->name);
      saft_search_add_result (engine->tmp_search, result);
    }

  return 1;
}

static int
search_engine_dna_array_cache_sequence (SaftSequence *sequence,
                                        void         *data)
{
  SearchEngineDNAArray *engine;
  DNAArrayDBEntry      *entry;

  engine        = (SearchEngineDNAArray*)data;
  entry         = dna_array_db_entry_new ();
  entry->counts = search_engine_dna_array_hash_sequence (engine, sequence);
  entry->name   = strdup (sequence->name);
  entry->length = sequence->seq_length;

  if (engine->search_engine.options->cache_db)
    {
      entry->next      = engine->db_cache;
      engine->db_cache = entry;
    }
  else
    {
      entry->next         = engine->query_cache;
      engine->query_cache = entry;
      engine->n_queries++;
    }

  return 1;
}

static SaftSearch*
search_engine_dna_array_search_all_dcached (SearchEngineDNAArray *engine,
                                            const char           *query_path,
                                            const char           *db_path)
{
  saft_fasta_iter (db_path,
                   search_engine_dna_array_cache_sequence,
                   engine);
  saft_fasta_iter (query_path,
                   search_engine_dna_array_search_query,
                   engine);

  return engine->search;
}

static int
search_engine_dna_array_search_query (SaftSequence *sequence,
                                      void         *data)
{
  SearchEngineDNAArray *engine;
  SaftSearch           *search;
  DNAArrayDBEntry      *entry;
  WordCount            *counts;

  engine       = (SearchEngineDNAArray*)data;
  counts       = search_engine_dna_array_hash_sequence (engine, sequence);
  search       = saft_search_new (engine->search_engine.options->max_results);
  search->name = strdup (sequence->name);

  for (entry = engine->db_cache; entry; entry = entry->next)
    {
      unsigned long d2 = 0;
      double        mean;
      double        var;
      int           i;

      for (i = 0; i < engine->max_words; i++)
        d2 += entry->counts[i] * counts[i];

      mean = saft_stats_mean (engine->stats_context,
                              sequence->seq_length,
                              entry->length);
      var  = saft_stats_mean (engine->stats_context,
                              sequence->seq_length,
                              entry->length);

      /* FIXME adjust this euristic depending on the user's required significance level */
      if (d2 < mean + 2 * sqrt (var))
          continue;
      else
        {
          SaftResult    *result;

          result          = saft_result_new ();
          result->d2      = d2;
          result->p_value = saft_stats_pgamma_m_v (result->d2, mean, var);
          result->name    = strdup (entry->name);
          saft_search_add_result (search, result);
        }
    }

  if (search->n_results > 0)
    {
      search->next = engine->search;
      engine->search = search;
    }
  else
    saft_search_free (search);

  return 1;
}

static SaftSearch*
search_engine_dna_array_search_all_qcached (SearchEngineDNAArray *engine,
                                            const char           *query_path,
                                            const char           *db_path)
{
  size_t i;

  saft_fasta_iter (query_path,
                   search_engine_dna_array_cache_sequence,
                   engine);
  engine->search_array = calloc (engine->n_queries, sizeof (*engine->search_array));
  saft_fasta_iter (db_path,
                   search_engine_dna_array_search_db,
                   engine);

  for (i = 0; i < engine->n_queries; i++)
    {
      SaftSearch *search;

      search = engine->search_array[i];
      if (!search || search->n_results == 0)
        continue;
      search->next   = engine->search;
      engine->search = search;
    }
  free (engine->search_array);
  engine->search_array = NULL;

  return engine->search;
}

static int
search_engine_dna_array_search_db (SaftSequence *sequence,
                                   void         *data)
{
  SearchEngineDNAArray *engine;
  DNAArrayDBEntry      *entry;
  WordCount            *counts;
  size_t                query_idx = 0;

  engine = (SearchEngineDNAArray*)data;
  counts = search_engine_dna_array_hash_sequence (engine, sequence);

  for (entry = engine->query_cache; entry; entry = entry->next)
    {
      unsigned long d2 = 0;
      double        mean;
      double        var;
      int           i;

      for (i = 0; i < engine->max_words; i++)
        d2 += entry->counts[i] * counts[i];

      mean = saft_stats_mean (engine->stats_context,
                              sequence->seq_length,
                              entry->length);
      var  = saft_stats_mean (engine->stats_context,
                              sequence->seq_length,
                              entry->length);

      if (d2 < mean + 2 * sqrt (var))
        {
          query_idx++;
          continue;
        }
      else
        {
          SaftSearch *search;
          SaftResult *result;

          if (!engine->search_array[query_idx])
            {
              search       = saft_search_new (engine->search_engine.options->max_results);
              /* TODO could use the same string as in the entry and deallocate carefully */
              search->name = strdup (entry->name);
              engine->search_array[query_idx] = search;
            }
          search = engine->search_array[query_idx];

          result          = saft_result_new ();
          result->d2      = d2;
          result->p_value = saft_stats_pgamma_m_v (result->d2, mean, var);
          result->name    = strdup (sequence->name);
          saft_search_add_result (search, result);
        }
      query_idx++;
    }

  return 1;
}

static DNAArrayDBEntry*
dna_array_db_entry_new ()
{
  DNAArrayDBEntry *entry;

  entry         = malloc (sizeof (*entry));
  entry->next   = NULL;
  entry->counts = NULL;
  entry->name   = NULL;
  entry->length = 0;

  return entry;
}

static void
dna_array_db_entry_free (DNAArrayDBEntry *entry)
{
  if (entry->counts)
    free (entry->counts);
  if (entry->name)
    free (entry->name);
  free (entry);
}

static void
dna_array_db_entry_free_all (DNAArrayDBEntry *entry)
{
  while (entry)
    {
      DNAArrayDBEntry *next;

      next  = entry->next;
      dna_array_db_entry_free (entry);
      entry = next;
    }
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
