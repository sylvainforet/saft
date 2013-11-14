/* saftsearch.c
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

#include "safterror.h"
#include "saftsearch.h"
#include "saftstats.h"


static void saft_search_adjust_pvalues (SaftSearch *search);

static void saft_search_sort_results   (SaftSearch *search);

static int  saft_result_pvalue_cmp_inc (const void *p1,
                                        const void *p2);


/******************/
/* Search Options */
/******************/


static SaftOptions*
saft_options_new ()
{
  SaftOptions *options;

  options                      = malloc (sizeof (*options));
  options->input_path          = NULL;
  options->db_path             = NULL;
  options->output_path         = NULL;
  options->alphabet            = NULL;
  options->letter_frequencies  = NULL;
  options->p_max               = 5e-2;
  options->word_size           = 0;
  options->verbosity           = 0;
  options->show_max            = 50;
  options->program             = SAFT_UNKNOWN_PROGRAM;
  options->freq_type           = SAFT_FREQ_UNIFORM;
  options->cache_db            = 0;
  options->cache_queries       = 0;

  return options;
}

static void
saft_options_free (SaftOptions *options)
{
  if (options)
    free (options);
}


/**********/
/* Result */
/**********/

SaftResult*
saft_result_new ()
{
  SaftResult *result;

  result               = malloc (sizeof (*result));
  result->next         = NULL;
  result->name         = NULL;
  result->d2           = 0;
  result->subject_size = 0;
  result->p_value      = 1;
  result->p_value_adj  = 1;

  return result;
}

void
saft_result_free (SaftResult *result)
{
  if (result)
    {
      saft_result_free (result->next);
      if (result->name)
        free (result->name);
      free (result);
    }
}

/**********/
/* Search */
/**********/

SaftSearch*
saft_search_new (SaftSequence *query,
                 unsigned int  word_size,
                 SaftFreqType  freq_type,
                 const double *letter_frequencies)
{
  SaftSearch *search;

  search                      = malloc (sizeof (*search));
  search->query               = query;
  search->word_size           = word_size;
  search->freq_type           = freq_type;
  search->htable              = saft_htable_new (query->alphabet, search->word_size);
  search->letter_frequencies  = malloc (query->alphabet->size * sizeof (*search->letter_frequencies));
  search->letters_counts      = malloc (query->alphabet->size * sizeof (*search->letters_counts));
  search->results             = NULL;
  search->sorted_results      = NULL;
  search->n_results           = 0;
  saft_htable_add_query (search->htable, search->query);

  if (search->freq_type == SAFT_FREQ_QUERY ||
      search->freq_type == SAFT_FREQ_QUERY_SUBJECTS)
    {
      SaftSegment *segment;

      for (segment = query->segments; segment; segment = segment->next)
        {
          unsigned int i;

          for (i = 0; i < segment->size; i++)
            /* (segment->seq[i] - 1) because there should not be any `0' letter
             * left at this stage */
            search->letters_counts[segment->seq[i] - 1]++;
        }
    }
  else if (search->freq_type == SAFT_FREQ_USER)
    {
      unsigned int i;

      for (i = 0; i < query->alphabet->size; i++)
        search->letter_frequencies[i] = letter_frequencies[i];
    }
  else /* Assume uniform (SAFT_FREQ_UNIFORM) by default */
    {
      unsigned int i;

      for (i = 0; i < query->alphabet->size; i++)
        search->letters_counts[i] = 1;
    }

  return search;
}

void
saft_search_free (SaftSearch *search)
{
  if (search)
    {
      if (search->query)
        saft_sequence_free (search->query);
      if (search->htable)
        saft_htable_free (search->htable);
      if (search->letter_frequencies)
        free (search->letter_frequencies);
      if (search->letters_counts)
        free (search->letters_counts);
      if (search->results)
        saft_result_free (search->results);
      if (search->sorted_results)
        free (search->sorted_results);
      free (search);
    }
}

void
saft_search_add_subject (SaftSearch   *search,
                         SaftSequence *subject)
{
  SaftResult *result;

  saft_htable_add_subject (search->htable, subject);
  result               = saft_result_new ();
  result->d2           = saft_htable_d2 (search->htable);
  result->name         = strdup (subject->name);
  result->subject_size = subject->size;
  result->next         = search->results;
  search->results      = result;
  search->n_results++;

  if (search->freq_type == SAFT_FREQ_SUBJECTS ||
      search->freq_type == SAFT_FREQ_QUERY_SUBJECTS)
    {
      SaftSegment *segment;

      for (segment = subject->segments; segment; segment = segment->next)
        {
          unsigned int i;

          for (i = 0; i < segment->size; i++)
            /* (segment->seq[i] - 1) because there should not be any `0' letter
             * left at this stage */
            search->letters_counts[segment->seq[i] - 1]++;
        }
    }
  saft_htable_clear_subject (search->htable);
}

void
saft_search_compute_pvalues (SaftSearch *search)
{
  SaftStatsContext *context;
  SaftResult       *result;
  unsigned int      i;

  if (search->freq_type != SAFT_FREQ_USER)
    {
      unsigned int total = 0;

      for (i = 0; i < search->query->alphabet->size; i++)
        total += search->letters_counts[i];
      for (i = 0; i < search->query->alphabet->size; i++)
        search->letter_frequencies[i] = ((double)search->letters_counts[i]) / total;
    }

  context = saft_stats_context_new (search->word_size,
                                    search->letter_frequencies,
                                    search->query->alphabet->size);

  for (result = search->results; result; result = result->next)
    {
      const double mean = saft_stats_mean (context,
                                           search->query->size,
                                           result->subject_size);
      const double var  = saft_stats_var  (context,
                                           search->query->size,
                                           result->subject_size);
      result->p_value   = saft_stats_pgamma_m_v (result->d2, mean, var);
    }

  saft_search_adjust_pvalues (search);
}

static void
saft_search_adjust_pvalues (SaftSearch *search)
{
  int i;

  saft_search_sort_results (search);

  search->sorted_results[search->n_results - 1]->p_value_adj = search->sorted_results[search->n_results - 1]->p_value;
  for (i = search->n_results - 2; i >= 0; i--)
    search->sorted_results[i]->p_value_adj = saft_stats_BH_element (search->sorted_results[i]->p_value,
                                                                    search->sorted_results[i + 1]->p_value_adj,
                                                                    i,
                                                                    search->n_results);
}

static void
saft_search_sort_results (SaftSearch *search)
{
  SaftResult  *result;
  SaftResult **tmp;

  search->sorted_results = malloc (search->n_results * sizeof (*search->sorted_results));
  tmp                    = search->sorted_results - 1;

  for (result = search->results; result; result = result->next)
    *++tmp = result;

  qsort (search->sorted_results,
         search->n_results,
         sizeof (*search->sorted_results),
         saft_result_pvalue_cmp_inc);
}

static int
saft_result_pvalue_cmp_inc (const void *p1,
                            const void *p2)
{
  const SaftResult **res1 = (const SaftResult**)p1;
  const SaftResult **res2 = (const SaftResult**)p2;

  if ((*res1)->p_value > (*res2)->p_value)
    return 1;
  return -1;
}


/********************/
/* SaftSearchEngine */
/********************/

SaftSearchEngine*
saft_search_engine_new (SaftOptions *options)
{
  if (options->alphabet == SaftAlphabetDNA)
    {
      if (options->word_size <= 8)
        {
          /* Array based DNA engine */
        }
      else
        {
          /* Hash-Table based DNA engine */
        }
    }
  else
    {
      /* Generic engine */
    }
}

void
saft_search_engine_free (SaftSearchEngine *engine)
{
  engine->free (engine);
}

SaftSearch*
saft_search_two_sequences (SaftSearchEngine *engine,
                           SaftSequence     *query,
                           SaftSequence     *subject)
{
  return engine->search_two_sequences (engine,
                                       query,
                                       subject);
}

SaftSearch**
saft_search_all (SaftSearchEngine *engine,
                 const char       *query_path,
                 const char       *db_path)
{
  return engine->search_all (engine,
                             query_path,
                             db_path);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
