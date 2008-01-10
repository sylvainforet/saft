/**
 *
 */

#include <stdlib.h>

#include "bfastsequence.h"


BfastAlphabet BfastAlphabetDNA =
{
  .name    = "DNA",
  .letters = "NATGC",
  .size    = 4,
  .codes   = 
    {
      ['A'] = 1,
      ['T'] = 2,
      ['G'] = 3,
      ['C'] = 4,
      ['a'] = 1,
      ['t'] = 2,
      ['g'] = 3,
      ['c'] = 4
    }
};

BfastAlphabet BfastAlphabetProtein =
{
  .name    = "Protein",
  .letters = "XACDEFGHIKLMNPQRSTVWY",
  .size    = 20,
  .codes   = 
    {
      ['A'] = 1,
      ['C'] = 2,
      ['D'] = 3,
      ['E'] = 4,
      ['F'] = 5,
      ['G'] = 6,
      ['H'] = 7,
      ['I'] = 8,
      ['K'] = 9,
      ['L'] = 10,
      ['M'] = 11,
      ['N'] = 12,
      ['P'] = 13,
      ['Q'] = 14,
      ['R'] = 15,
      ['S'] = 16,
      ['T'] = 17,
      ['V'] = 18,
      ['W'] = 19,
      ['Y'] = 20,
      ['a'] = 1,
      ['c'] = 2,
      ['d'] = 3,
      ['e'] = 4,
      ['f'] = 5,
      ['g'] = 6,
      ['h'] = 7,
      ['i'] = 8,
      ['k'] = 9,
      ['l'] = 10,
      ['m'] = 11,
      ['n'] = 12,
      ['p'] = 13,
      ['q'] = 14,
      ['r'] = 15,
      ['s'] = 16,
      ['t'] = 17,
      ['v'] = 18,
      ['w'] = 19,
      ['y'] = 20
    }
};

BfastAlphabet*
bfast_alphabet_new ()
{
  BfastAlphabet *alphabet;

  alphabet          = malloc (sizeof (*alphabet));
  alphabet->name    = NULL;
  alphabet->size    = 0;
  alphabet->letters = NULL;

  return alphabet;
}

void
bfast_alphabet_free (BfastAlphabet *alphabet)
{
  if (alphabet)
    {
      if (alphabet->name)
        free (alphabet->name);
      if (alphabet->letters)
        free (alphabet->letters);
      free (alphabet);
    }
}

BfastSequence*
bfast_sequence_new (void)
{
  BfastSequence *seq;

  seq           = malloc (sizeof (*seq));
  seq->name     = NULL;
  seq->seq      = NULL;
  seq->alphabet = NULL;
  seq->size     = 0;

  return seq;
}

void
bfast_sequence_free (BfastSequence *seq)
{
  if (seq)
    {
      if (seq->name)
        free (seq->name);
      if (seq->seq)
        free (seq->seq);
      free (seq);
    }
}

char*
bfast_sequence_to_string (BfastSequence *seq)
{
  char*        ret;
  unsigned int i;

  ret = malloc ((seq->size + 1) * sizeof (*ret));
  for (i = 0; i < seq->size; i++)
    ret[i] = seq->alphabet->letters[seq->seq[i]];
  ret[i] = '\0';

  return ret;
}

void
bfast_sequences_iter (BfastSequence  **sequences,
                      BfastSeqIterFunc func,
                      void            *data)
{
  while (*sequences)
    {
      if (!func (*sequences, data))
        return;
      sequences++;
    }
}
