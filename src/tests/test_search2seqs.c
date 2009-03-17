/* test_search2seqs.c
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

#include <stdio.h>
#include <stdlib.h>

#include "wafterror.h"
#include "waftfasta.h"
#include "waftsearch.h"

int
main (int argc, char **argv)
{
  WaftFasta   **seqs;
  WaftFasta   **tmp;
  unsigned int   word_size;
  unsigned int   n;

  if (argc < 2)
    {
      waft_error ("A fasta file and a word size need to be given as argument");
      return 1;
    }

  word_size = atoi (argv[2]);
  seqs = waft_fasta_read (argv[1], &n);

  if (n < 2)
    {
      waft_error ("The fasta file must contain at least two sequences");
      return 1;
    }
  else
    {
      WaftSequence      *s1  = waft_fasta_to_seq (seqs[0], &WaftAlphabetDNA);
      WaftSequence      *s2  = waft_fasta_to_seq (seqs[1], &WaftAlphabetDNA);
      WaftSearchOptions *opt = waft_search_options_new ();
      WaftSearch        *s   = waft_search_new ();

      opt->word_size = word_size;
      s->options     = opt;
      s->query       = s1;
      s->subject     = s2;

      waft_search_process (s);
      printf ("D2 == %ld\n", s->d2);

      waft_sequence_free (s1);
      waft_sequence_free (s2);
      waft_search_options_free (opt);
      waft_htable_free (s->query_h);
      waft_htable_free (s->subject_h);
      waft_search_free (s);
    }

  tmp = seqs - 1;
  while (*++tmp)
    waft_fasta_free (*tmp);
  free (seqs);

  return 0;
}
