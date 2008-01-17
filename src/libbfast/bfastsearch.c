/* bfastsearch.c
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

#include "bfastsearch.h"


BfastSearchOptions*
bfast_search_options_new ()
{
  BfastSearchOptions *options;

  options            = malloc (sizeof (*options));
  options->word_size = 0;

  return options;
}

void
bfast_search_options_free (BfastSearchOptions *options)
{
  if (options)
    free (options);
}

BfastSearch*
bfast_search_new ()
{
  BfastSearch *search;

  search            = malloc (sizeof (*search));
  search->options   = NULL;
  search->query     = NULL;
  search->subject   = NULL;
  search->query_h   = NULL;
  search->subject_h = NULL;
  search->d2        = 0;

  return search;
}

void
bfast_search_free (BfastSearch *search)
{
  if (search)
    free (search);
}

void
bfast_search_process (BfastSearch *search)
{
  if (!search->query_h)
    {
      search->query_h = bfast_htable_new (search->query->alphabet,
                                          search->options->word_size);
      bfast_htable_add_seq (search->query_h, search->query);
    }
  if (!search->subject_h)
    {
      search->subject_h = bfast_htable_new (search->query->alphabet,
                                            search->options->word_size);
      bfast_htable_add_seq (search->subject_h, search->subject);
    }
  search->d2 = bfast_htable_d2 (search->subject_h, search->query_h);
}
