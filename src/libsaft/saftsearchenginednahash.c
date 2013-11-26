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

#include "safterror.h"
#include "saftfasta.h"
#include "safthash.h"
#include "saftsearchengines.h"
#include "saftstats.h"


/********************************/
/* Hash-based DNA Search Engine */
/********************************/

typedef struct _DNAHashDBEntry DNAHashDBEntry;

struct _DNAHashDBEntry
{
  DNAHashDBEntry *next;
  SaftHashTable  *counts;
  char           *name;
  size_t          length;
};

static DNAHashDBEntry*  dna_hash_db_entry_new      (void);

static void             dna_hash_db_entry_free     (DNAHashDBEntry *entry);

static void             dna_hash_db_entry_free_all (DNAHashDBEntry *entry);


typedef struct _SearchEngineDNAHash SearchEngineDNAHash;

struct _SearchEngineDNAHash
{
  SaftSearchEngine  search_engine;

  SaftStatsContext *stats_context;

  DNAHashDBEntry   *query_cache;
  DNAHashDBEntry   *db_cache;

  SaftSearch       *search;
  SaftSearch      **search_array;
  SaftSearch       *tmp_search;
  SaftHashTable    *tmp_counts;

  size_t            n_queries;
  size_t            tmp_length;
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

static SaftSearch*   search_engine_dna_hash_search_all_no_cache   (SearchEngineDNAHash  *engine,
                                                                   const char           *query_path,
                                                                   const char           *db_path);

static SaftSearch*   search_engine_dna_hash_search_all_qcached    (SearchEngineDNAHash  *engine,
                                                                   const char           *query_path,
                                                                   const char           *db_path);

static SaftSearch*   search_engine_dna_hash_search_all_dcached    (SearchEngineDNAHash  *engine,
                                                                   const char           *query_path,
                                                                   const char           *db_path);

static int           search_engine_dna_hash_cache_sequence        (SaftSequence         *sequence,
                                                                   void                 *data);

static int           search_engine_dna_hash_search_query          (SaftSequence         *sequence,
                                                                   void                 *data);

static int           search_engine_dna_hash_search_db             (SaftSequence         *sequence,
                                                                   void                 *data);

static int           search_engine_dna_hash_queries_iter_func     (SaftSequence         *sequence,
                                                                   void                 *data);

static int           search_engine_dna_hash_db_iter_func          (SaftSequence         *sequence,
                                                                   void                 *data);

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

  engine->query_cache                        = NULL;
  engine->db_cache                           = NULL;
  engine->search                             = NULL;
  engine->search_array                       = NULL;
  engine->tmp_search                         = NULL;
  engine->tmp_counts                         = NULL;
  engine->n_queries                          = 0;
  engine->tmp_length                         = 0;

  return (SaftSearchEngine*)engine;
}

static void
search_engine_dna_hash_free (SaftSearchEngine *engine)
{
  SearchEngineDNAHash *se;

  se = (SearchEngineDNAHash*) engine;

  if (se->stats_context)
    saft_stats_context_free (se->stats_context);
  if (se->query_cache)
    dna_hash_db_entry_free_all (se->query_cache);
  if (se->db_cache)
    dna_hash_db_entry_free_all (se->db_cache);
  if (se->search)
    saft_search_free (se->search);
  if (se->search_array)
    free (se->search_array);

  free (se);
}

static SaftHashTable*
search_engine_dna_hash_hash_sequence (SearchEngineDNAHash *engine,
                                      SaftSequence        *sequence)
{
  SaftHashTable *table;
  const size_t   k = engine->search_engine.options->word_size;
  size_t         i;

  table = saft_hash_table_new (engine->search_engine.options->word_size);

  if (k <= KMER_VAL_NUCS)
    {
      /* FIXME this should be computed only once for an engine instance */
      const unsigned long mask = (~ 0ul) >> (8 * sizeof (unsigned long) - (2 * k));
      SaftHashKmer        kmer;

      for (i = 0; i < k; i++)
        {
          const unsigned char c = SaftAlphabetDNA.codes[(int)sequence->seq[i]];

          kmer.kmer_vall <<= 2;
          kmer.kmer_vall |= c;
        }
      saft_hash_table_increment (table, &kmer);
      for (i = k; i < sequence->seq_length; i++)
        {
          const unsigned char c = SaftAlphabetDNA.codes[(int)sequence->seq[i]];

          kmer.kmer_vall <<= 2;
          kmer.kmer_vall |= c;
          kmer.kmer_vall &= mask;
          saft_hash_table_increment (table, &kmer);
        }
    }
  else
    {
      saft_error ("[ERROR] saftn with words > %dbp not implemented", KMER_VAL_NUCS);
      exit (1);
    }

  return table;
}

static unsigned long
search_engine_dna_hash_d2 (SearchEngineDNAHash *engine,
                           SaftHashTable       *counts1,
                           SaftHashTable       *counts2)
{
  SaftHashTable *small_table;
  SaftHashTable *large_table;
  unsigned long  d2 = 0;
  size_t         i;

  small_table = counts1;
  large_table = counts2;
  if (small_table->nnodes > large_table->nnodes)
    {
      small_table = counts2;
      large_table = counts1;
    }

  for (i = 0; i < small_table->size; i++)
    {
      SaftHashNode  *node1;
      SaftHashNode  *node2;
      unsigned long  node_index;
      unsigned long  step = 0;

      node1 = small_table->nodes + i;
      if (node1->key_hash <= 1)
        continue;

      node_index = node1->key_hash % large_table->mod;
      node2      = large_table->nodes + node_index;

      while (node2->key_hash)
        {
          if (node2->key_hash == node1->key_hash)
            if (large_table->key_equal_func (&node1->kmer, &node2->kmer, large_table->kmer_bytes))
              {
                d2 += node1->value.count * node2->value.count;
                break;
              }
          step++;
          node_index += saft_hash_table_probe (step);
          node_index &= large_table->mask;
          node2       = large_table->nodes + node_index;
        }
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
  var                  = saft_stats_var (se->stats_context,
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
  SaftSearch          *search = NULL;

  se = (SearchEngineDNAHash*) engine;

  if (engine->options->cache_db)
    search = search_engine_dna_hash_search_all_dcached (se, query_path, db_path);
  else if (engine->options->cache_queries)
    search = search_engine_dna_hash_search_all_qcached (se, query_path, db_path);
  else
    search = search_engine_dna_hash_search_all_no_cache (se, query_path, db_path);

  se->search = NULL;
  return search;
}

static SaftSearch*
search_engine_dna_hash_search_all_no_cache (SearchEngineDNAHash *engine,
                                            const char          *query_path,
                                            const char          *db_path)
{
  saft_fasta_iter (query_path,
                   search_engine_dna_hash_queries_iter_func,
                   engine);

  engine->search = saft_search_reverse (engine->search);

  return engine->search;
}

static int
search_engine_dna_hash_queries_iter_func (SaftSequence *sequence,
                                          void         *data)
{
  SearchEngineDNAHash *engine;

  engine = (SearchEngineDNAHash*)data;

  engine->tmp_counts       = search_engine_dna_hash_hash_sequence (engine, sequence);
  engine->tmp_length       = sequence->seq_length;
  engine->tmp_search       = saft_search_new (engine->search_engine.options->max_results);
  engine->tmp_search->name = strdup (sequence->name);

  saft_fasta_iter (engine->search_engine.options->db_path,
                   search_engine_dna_hash_db_iter_func,
                   engine);

  saft_hash_table_destroy (engine->tmp_counts);
  engine->tmp_counts = NULL;

  if (engine->tmp_search->n_results == 0)
    saft_search_free (engine->tmp_search);
  else
    {
      engine->tmp_search->next = engine->search;
      engine->search           = engine->tmp_search;
    }
  engine->tmp_search = NULL;

  return 1;
}

static int
search_engine_dna_hash_db_iter_func (SaftSequence *sequence,
                                     void         *data)
{
  SearchEngineDNAHash *engine;
  SaftHashTable       *counts;
  unsigned long        d2;
  double               mean;
  double               var;

  engine = (SearchEngineDNAHash*)data;
  counts = search_engine_dna_hash_hash_sequence (engine, sequence);
  d2     = search_engine_dna_hash_d2 (engine,
                                       engine->tmp_counts,
                                       counts);
  mean   = saft_stats_mean (engine->stats_context,
                            sequence->seq_length,
                            engine->tmp_length);
  var    = saft_stats_var (engine->stats_context,
                           sequence->seq_length,
                           engine->tmp_length);
  saft_hash_table_destroy (counts);

  /* FIXME adjust this euristic depending on the user's required significance level */
  /* Fix this here and everywhere else in this file */
  if (d2 > mean + 2 * sqrt (var))
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
search_engine_dna_hash_cache_sequence (SaftSequence *sequence,
                                       void         *data)
{
  SearchEngineDNAHash *engine;
  DNAHashDBEntry     *entry;

  engine        = (SearchEngineDNAHash*)data;
  entry         = dna_hash_db_entry_new ();
  entry->counts = search_engine_dna_hash_hash_sequence (engine, sequence);
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
search_engine_dna_hash_search_all_qcached (SearchEngineDNAHash *engine,
                                           const char          *query_path,
                                           const char          *db_path)
{
  unsigned long i;

  saft_fasta_iter (query_path,
                   search_engine_dna_hash_cache_sequence,
                   engine);
  engine->search_array = calloc (engine->n_queries, sizeof (*engine->search_array));
  saft_fasta_iter (db_path,
                   search_engine_dna_hash_search_db,
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
search_engine_dna_hash_search_db (SaftSequence *sequence,
                                  void         *data)
{
  SearchEngineDNAHash *engine;
  DNAHashDBEntry      *entry;
  SaftHashTable       *counts;
  size_t               query_idx = 0;

  engine = (SearchEngineDNAHash*)data;
  counts = search_engine_dna_hash_hash_sequence (engine, sequence);

  for (entry = engine->query_cache; entry; entry = entry->next)
    {
      unsigned long d2;
      double        mean;
      double        var;

      d2   = search_engine_dna_hash_d2 (engine,
                                        entry->counts,
                                        counts);
      mean = saft_stats_mean (engine->stats_context,
                              sequence->seq_length,
                              entry->length);
      var  = saft_stats_var (engine->stats_context,
                             sequence->seq_length,
                             entry->length);

      if (d2 > mean + 2 * sqrt (var))
        {
          SaftSearch *search;
          SaftResult *result;

          if (!engine->search_array[query_idx])
            {
              search       = saft_search_new (engine->search_engine.options->max_results);
              /* TODO could use the same string as in the entry and deallocate
               * carefully (probably not worth the trouble) */
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

  saft_hash_table_destroy (counts);

  return 1;
}

static SaftSearch*
search_engine_dna_hash_search_all_dcached (SearchEngineDNAHash *engine,
                                           const char          *query_path,
                                           const char          *db_path)
{
  saft_fasta_iter (db_path,
                   search_engine_dna_hash_cache_sequence,
                   engine);
  saft_fasta_iter (query_path,
                   search_engine_dna_hash_search_query,
                   engine);

  engine->search = saft_search_reverse (engine->search);

  return engine->search;
}

static int
search_engine_dna_hash_search_query (SaftSequence         *sequence,
                                     void                 *data)
{
  SearchEngineDNAHash *engine;
  SaftSearch          *search;
  DNAHashDBEntry      *entry;
  SaftHashTable       *counts;

  engine       = (SearchEngineDNAHash*) data;
  counts       = search_engine_dna_hash_hash_sequence (engine, sequence);
  search       = saft_search_new (engine->search_engine.options->max_results);
  search->name = strdup (sequence->name);

  for (entry = engine->db_cache; entry; entry = entry->next)
    {
      unsigned long d2;
      double        mean;
      double        var;

      d2   = search_engine_dna_hash_d2 (engine,
                                        entry->counts,
                                        counts);
      mean = saft_stats_mean (engine->stats_context,
                              sequence->seq_length,
                              entry->length);
      var  = saft_stats_var (engine->stats_context,
                             sequence->seq_length,
                             entry->length);

      /* FIXME adjust this euristic depending on the user's required significance level */
      if (d2 > mean + 2 * sqrt (var))
        {
          SaftResult    *result;

          result          = saft_result_new ();
          result->d2      = d2;
          result->p_value = saft_stats_pgamma_m_v (result->d2, mean, var);
          result->name    = strdup (entry->name);
          saft_search_add_result (search, result);
        }
    }

  saft_hash_table_destroy (counts);

  if (search->n_results == 0)
    saft_search_free (search);
  else
    {
      search->next = engine->search;
      engine->search = search;
    }

  return 1;
}

static DNAHashDBEntry*
dna_hash_db_entry_new ()
{
  DNAHashDBEntry *entry;

  entry         = malloc (sizeof (*entry));
  entry->next   = NULL;
  entry->counts = NULL;
  entry->name   = NULL;
  entry->length = 0;

  return entry;
}

static void
dna_hash_db_entry_free (DNAHashDBEntry *entry)
{
  if (entry->counts)
    saft_hash_table_destroy (entry->counts);
  if (entry->name)
    free (entry->name);
  free (entry);
}

static void
dna_hash_db_entry_free_all (DNAHashDBEntry *entry)
{
  while (entry)
    {
      DNAHashDBEntry *next;

      next  = entry->next;
      dna_hash_db_entry_free (entry);
      entry = next;
    }
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
