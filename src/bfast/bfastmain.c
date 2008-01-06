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
      bfast_error ("No argument given");
      return 1;
    }
  seqs = bfast_fasta_read (argv[1], &n);

  tmp = seqs - 1;
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
