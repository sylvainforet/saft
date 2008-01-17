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

#include "bfasterror.h"
#include "bfastfasta.h"
#include "bfastsearch.h"

int
main (int argc, char **argv)
{
  BfastFasta   **seqs;
  BfastFasta   **tmp;
  unsigned int   word_size;
  unsigned int   n;

  if (argc < 2)
    {
      bfast_error ("A fasta file and a word size need to be given as argument");
      return 1;
    }

  word_size = atoi (argv[2]);
  seqs = bfast_fasta_read (argv[1], &n);

  if (n < 2)
    {
      bfast_error ("The fasta file must contain at least two sequences");
      return 1;
    }
  else
    {
      BfastSequence      *s1  = bfast_fasta_to_seq (seqs[0], &BfastAlphabetDNA);
      BfastSequence      *s2  = bfast_fasta_to_seq (seqs[1], &BfastAlphabetDNA);
      BfastSearchOptions *opt = bfast_search_options_new ();
      BfastSearch        *s   = bfast_search_new ();

      opt->word_size = word_size;
      s->options     = opt;
      s->query       = s1;
      s->subject     = s2;

      bfast_search_process (s);
      printf ("D2 == %ld\n", s->d2);

      bfast_sequence_free (s1);
      bfast_sequence_free (s2);
      bfast_search_options_free (opt);
      bfast_htable_free (s->query_h);
      bfast_htable_free (s->subject_h);
      bfast_search_free (s);
    }

  tmp = seqs - 1;
  while (*++tmp)
    bfast_fasta_free (*tmp);
  free (seqs);

  return 0;
}
