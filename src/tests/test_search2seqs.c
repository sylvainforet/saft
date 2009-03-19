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

#include "safterror.h"
#include "saftfasta.h"
#include "saftsearch.h"

int
main (int argc, char **argv)
{
  SaftFasta   **seqs;
  SaftFasta   **tmp;
  unsigned int  word_size;
  unsigned int  n;

  if (argc < 3)
    {
      saft_error ("Usage: %s FASTA_FILE WORD_SIZE", argv[0]);
      return 1;
    }

  seqs = saft_fasta_read (argv[1], &n);
  word_size = atoi (argv[2]);

  if (n < 2)
    {
      saft_error ("The fasta file must contain at least two sequences");
      return 1;
    }
  else
    {
      SaftSequence      *s1  = saft_fasta_to_seq (seqs[0], &SaftAlphabetDNA);
      SaftSequence      *s2  = saft_fasta_to_seq (seqs[1], &SaftAlphabetDNA);
      SaftSearchOptions *opt = saft_search_options_new ();
      SaftSearch        *s   = saft_search_new ();

      opt->word_size = word_size;
      s->options     = opt;
      s->query       = s1;
      s->subject     = s2;

      saft_search_process (s);
      printf ("D2 == %ld\n", s->d2);

      saft_sequence_free (s1);
      saft_sequence_free (s2);
      saft_search_options_free (opt);
      saft_htable_free (s->htable);
      saft_search_free (s);
    }

  tmp = seqs - 1;
  while (*++tmp)
    saft_fasta_free (*tmp);
  free (seqs);

  return 0;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
