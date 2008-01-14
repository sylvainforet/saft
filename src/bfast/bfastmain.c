/**
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "bfasterror.h"
#include "bfastfasta.h"
#include "bfastsearch.h"


static void test_load_print (const char  *path);

static void test_d2         (const char  *path,
                             unsigned int word_size);


int
main (int argc, char **argv)
{
  if (argc < 2)
    {
      bfast_error ("No argument given");
      return 1;
    }
  if (argc == 2)
    test_load_print (argv[1]);
  if (argc == 3)
    {
      int word_size = atoi (argv[2]);
      test_d2 (argv[1], word_size);
    }

  return 0;
}

static void
test_load_print (const char* path)
{
  BfastFasta   **seqs;
  BfastFasta   **tmp;
  unsigned int   n;

  seqs = bfast_fasta_read (path, &n);

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
}

static void
test_d2 (const char  *path,
         unsigned int word_size)
{
  BfastFasta   **seqs;
  BfastFasta   **tmp;
  unsigned int   n;

  seqs = bfast_fasta_read (path, &n);

  if (n < 2)
    {
      bfast_error ("Two sequences are needed");
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
}
