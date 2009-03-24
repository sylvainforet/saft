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
      SaftSequence *s1     = saft_fasta_to_seq (seqs[0], &SaftAlphabetDNA);
      SaftSequence *s2     = saft_fasta_to_seq (seqs[1], &SaftAlphabetDNA);
      SaftSearch   *search = saft_search_new (s1, word_size, SAFT_FREQ_UNIFORM, NULL);

      saft_search_add_subject (search, s2);
      saft_search_compute_pvalues (search);
      printf ("D2 = %d ; p-value = %.5e\n", search->results->d2, search->results->p_value);

      saft_sequence_free (s2);
      saft_search_free (search);
    }

  tmp = seqs - 1;
  while (*++tmp)
    saft_fasta_free (*tmp);
  free (seqs);

  return 0;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
