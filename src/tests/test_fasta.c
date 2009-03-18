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

#include "wafterror.h"
#include "waftfasta.h"

int
main (int argc, char **argv)
{
  WaftFasta   **seqs;
  WaftFasta   **tmp;
  unsigned int  n;

  if (argc < 2)
    {
      waft_error ("A fasta file need to be given as argument");
      return 1;
    }

  seqs = waft_fasta_read (argv[1], &n);
  tmp  = seqs - 1;
  while (*++tmp)
    {
      WaftSequence *seq;
      char          *str;

      seq = waft_fasta_to_seq (*tmp, &WaftAlphabetDNA);
      printf (">%s\n", seq->name);
      str = waft_sequence_to_string (seq);
      printf ("%s\n", str);
      free (str);
      waft_sequence_free (seq);
    }

  tmp = seqs - 1;
  while (*++tmp)
    waft_fasta_free (*tmp);
  free (seqs);

  return 0;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
