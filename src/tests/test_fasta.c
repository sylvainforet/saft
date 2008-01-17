/* test_fasta.c
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

int
main (int argc, char **argv)
{
  BfastFasta   **seqs;
  BfastFasta   **tmp;
  unsigned int   n;

  if (argc < 2)
    {
      bfast_error ("A fasta file need to be given as argument");
      return 1;
    }

  seqs = bfast_fasta_read (argv[1], &n);
  tmp  = seqs - 1;
  while (*++tmp)
    {
      BfastSequence *seq;
      char          *str;

      seq = bfast_fasta_to_seq (*tmp, &BfastAlphabetDNA);
      printf (">%s\n", seq->name);
      str = bfast_sequence_to_string (seq);
      printf ("%s\n", str);
      free (str);
      bfast_sequence_free (seq);
    }

  tmp = seqs - 1;
  while (*++tmp)
    bfast_fasta_free (*tmp);
  free (seqs);

  return 0;
}
