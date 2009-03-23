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

#include "saftsearch.h"


/**********/
/* Result */
/**********/

SaftResult*
saft_result_new ()
{
  SaftResult *result;

  result               = malloc (sizeof (*result));
  result->next         = NULL;
  result->d2           = 0;
  result->subject_size = 0;

  return result;
}

void
saft_result_free (SaftResult *result)
{
  if (result)
    {
      saft_result_free (result->next);
      free (result);
    }
}

/**********/
/* Search */
/**********/

SaftSearch*
saft_search_new (SaftSequence *query,
                 unsigned int  word_size,
                 SaftFreqType  freq_type)
{
  SaftSearch *search;

  search                      = malloc (sizeof (*search));
  search->query               = query;
  search->word_size           = word_size;
  search->freq_type           = freq_type;
  search->htable              = saft_htable_new (query->alphabet, search->word_size);
  search->letters_frequencies = malloc (query->alphabet->size * sizeof (*search->letters_frequencies));
  search->letters_counts      = malloc (query->alphabet->size * sizeof (*search->letters_counts));
  saft_htable_add_query (search->htable, search->query);

  if (search->freq_type == SAFT_FREQ_QUERY ||
      search->freq_type == SAFT_FREQ_QUERY_SUBJECTS)
    {
      SaftSegment *segment;

      for (segment = query->segments; segment; segment = segment->next)
        {
          unsigned int i;

          for (i = 0; i < segment->size; i++)
            search->letters_counts[segment->seq[i]]++;
        }
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
      if (search->letters_frequencies)
        free (search->letters_frequencies);
      if (search->letters_counts)
        free (search->letters_counts);
      if (search->results)
        saft_result_free (search->results);
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
  result->subject_size = subject->size;
  result->next         = search->results;
  search->results      = result;

  if (search->freq_type == SAFT_FREQ_SUBJECTS ||
      search->freq_type == SAFT_FREQ_QUERY_SUBJECTS)
    {
      SaftSegment *segment;

      for (segment = subject->segments; segment; segment = segment->next)
        {
          unsigned int i;

          for (i = 0; i < segment->size; i++)
            search->letters_counts[segment->seq[i]]++;
        }
    }
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
