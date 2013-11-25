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
#include "saftsearchengines.h"
#include "saftstats.h"


const char *saft_program_names[NB_SAFT_PROGRAMS] =
{
  [SAFT]   = "saft",
  [SAFTN]  = "saftn",
  [SAFTP]  = "saftp",
  [SAFTX]  = "saftx",
  [TSAFTN] = "tsaftn",
  [TSAFTX] = "tsaftx"
};


/* Priority queue functions and macros */

#define results_heap_parent(i) ((i) / 2)

#define results_heap_left(i)   (2 * (i))

#define results_heap_right(i)  (2 * (i) + 1)

static void results_heap_insert      (SaftSearch *search,
                                      SaftResult *result);

static void results_heap_heapify     (SaftSearch *search);

static void results_heap_sort        (SaftSearch *search);


/******************/
/* Search Options */
/******************/


SaftOptions*
saft_options_new ()
{
  SaftOptions *options;

  options                               = malloc (sizeof (*options));
  options->input_path                   = NULL;
  options->db_path                      = NULL;
  options->output_path                  = NULL;
  options->alphabet                     = NULL;
  options->letter_frequencies           = NULL;
  options->p_max                        = 5e-2;
  options->word_size                    = 0;
  options->verbosity                    = 0;
  options->max_results                  = 50;
  options->program                      = SAFT_UNKNOWN_PROGRAM;
  options->freq_type                    = SAFT_FREQ_UNIFORM;
  options->cache_db                     = 0;
  options->cache_queries                = 0;
  options->periodic_boundary_conditions = 0;

  return options;
}

void
saft_options_free (SaftOptions *options)
{
  if (options)
    {
      if (options->input_path)
        free (options->input_path);
      if (options->db_path)
        free (options->db_path);
      if (options->output_path)
        free (options->output_path);
      if (options->letter_frequencies)
        free (options->letter_frequencies);

      free (options);
    }
}


/**********/
/* Result */
/**********/

SaftResult*
saft_result_new ()
{
  SaftResult *result;

  result               = malloc (sizeof (*result));
  result->name         = NULL;
  result->d2           = 0;
  result->p_value      = 1;
  result->p_value_adj  = 1;

  return result;
}

void
saft_result_free (SaftResult *result)
{
  if (result)
    {
      if (result->name)
        free (result->name);

      free (result);
    }
}

/**********/
/* Search */
/**********/

SaftSearch*
saft_search_new (unsigned int max_results)
{
  SaftSearch *search;

  search                 = malloc (sizeof (*search));
  search->results        = calloc (max_results, sizeof (*search->results));
  search->name           = NULL;
  search->n_results      = 0;
  search->max_results    = max_results;

  return search;
}

void
saft_search_free (SaftSearch *search)
{
  if (search)
    {
      if (search->name)
        free (search->name);

      if (search->results)
        {
          int i;

          for (i = 0; i < search->max_results; i++)
            if (search->results[i])
              saft_result_free (search->results[i]);

          free (search->results);
        }

      free (search);
    }
}

void
saft_search_free_all (SaftSearch *search)
{
  while (search)
    {
      SaftSearch *next;

      next = search->next;
      saft_search_free (search);
      search = next;
    }
}

SaftSearch*
saft_search_reverse (SaftSearch *search)
{
  SaftSearch *prev = NULL;

  while (search)
    {
      SaftSearch *next;

      next         = search->next;
      search->next = prev;
      prev         = search;
      search       = next;
    }

  return prev;
}

void
saft_search_add_result (SaftSearch *search,
                        SaftResult *result)
{
  /* TODO Deal with ties */

  if (search->n_results == search->max_results)
    {
      if (result->p_value > search->results[0]->p_value)
        {
          saft_result_free (result);
          return;
        }
      saft_result_free (search->results[0]);
      search->results[0] = search->results[search->n_results - 1];
      search->results[search->n_results - 1] = NULL;
      search->n_results--;
      results_heap_heapify (search);
    }
  /* Insert new value */
  results_heap_insert (search, result);
}

void
saft_search_adjust_pvalues (SaftSearch *search)
{
  int i;

  if (search->n_results == 0)
    return;

  results_heap_sort (search);

  /* FIXME Make sure the direction of the sorting is correct */
  /* FIXME Have a conservative adjustment of the first p-value */
  search->results[search->n_results - 1]->p_value_adj = search->results[search->n_results - 1]->p_value;
  for (i = search->n_results - 2; i >= 0; i--)
    search->results[i]->p_value_adj = saft_stats_BH_element (search->results[i]->p_value,
                                                                    search->results[i + 1]->p_value_adj,
                                                                    i,
                                                                    search->n_results);
}

/* TODO Try alternative data structure for the results, maybe a Fibonacci heap */

/* FIXME Fix up the indices to have a consistent zero-based indexing system.
 * Importantly, make sure that this does not break the parent/left/right macros
 * */

static void
results_heap_insert (SaftSearch *search,
                     SaftResult *result)
{
  int i;
  int p;

  search->results[search->n_results] = result;
  search->n_results++;

  /* Heap-Increase-Key */
  i = search->n_results;
  p = results_heap_parent (i);
  while (i > 1 && search->results[p - 1]->p_value < search->results[i - 1]->p_value)
    {
      SaftResult *tmp;

      tmp                    = search->results[i - 1];
      search->results[i - 1] = search->results[p - 1];
      search->results[p - 1] = tmp;
      i                      = p;
      p                      = results_heap_parent (i);
    }
}

static void
results_heap_heapify (SaftSearch *search)
{
  int i;

  i = 1;
  while (i <= search->n_results)
    {
      const int l   = results_heap_left (i);
      const int r   = results_heap_right (i);
      int       max = i;

      if (l <= search->n_results && search->results[l - 1]->p_value > search->results[i - 1]->p_value)
        max = l;
      if (r <= search->n_results && search->results[r - 1]->p_value > search->results[max - 1]->p_value)
        max = r;

      if (max == i)
        break;
      else
        {
          SaftResult *tmp;

          tmp                      = search->results[i - 1];
          search->results[i - 1]   = search->results[max - 1];
          search->results[max - 1] = tmp;
          i                        = max;
        }
    }
}

static void
results_heap_sort (SaftSearch *search)
{
  const int n = search->n_results;
  int       i;

  for (i = n; i >= 2; i--)
    {
      SaftResult *tmp;

      tmp                    = search->results[i - 1];
      search->results[i - 1] = search->results[0];
      search->results[0]     = tmp;
      search->n_results--;
      results_heap_heapify (search);
    }
  search->n_results = n;
}

/********************/
/* SaftSearchEngine */
/********************/

SaftSearchEngine*
saft_search_engine_new (SaftOptions *options)
{
  if (options->program == SAFTN)
    {
      if (options->word_size <= 8)
        {
          return saft_search_engine_dna_array_new (options);
        }
      else
        {
          /* Hash-Table based DNA engine */
          saft_error ("Word sizes > 8 not implemented yet");
        }
    }
  else
    {
      /* Generic engine */
      saft_error ("Only DNA alphabet is implemented");
    }
  return NULL;
}

void
saft_search_engine_free (SaftSearchEngine *engine)
{
  if (engine->free)
    engine->free (engine);
}

SaftSearch*
saft_search_two_sequences (SaftSearchEngine *engine,
                           SaftSequence     *query,
                           SaftSequence     *subject)
{
  if (engine->search_two_sequences)
    return engine->search_two_sequences (engine,
                                         query,
                                         subject);
  return NULL;
}

SaftSearch*
saft_search_all (SaftSearchEngine *engine,
                 const char       *query_path,
                 const char       *db_path)
{
  if (engine->search_all)
    return engine->search_all (engine,
                               query_path,
                               db_path);
  return NULL;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
