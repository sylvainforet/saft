/* waftsequence.c
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

#include "waftsequence.h"


WaftAlphabet WaftAlphabetDNA =
{
  .name    = "DNA",
  .letters = "NATGC",
  .size    = 4,
  .codes   = 
    {
      ['N'] = 0,
      ['A'] = 1,
      ['T'] = 2,
      ['G'] = 3,
      ['C'] = 4,
      ['n'] = 0,
      ['a'] = 1,
      ['t'] = 2,
      ['g'] = 3,
      ['c'] = 4
    }
};

WaftAlphabet WaftAlphabetProtein =
{
  .name    = "Protein",
  .letters = "XACDEFGHIKLMNPQRSTVWY",
  .size    = 20,
  .codes   = 
    {
      ['X'] = 0,
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
      ['n'] = 0,
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

WaftAlphabet*
waft_alphabet_new ()
{
  WaftAlphabet *alphabet;

  alphabet          = malloc (sizeof (*alphabet));
  alphabet->name    = NULL;
  alphabet->size    = 0;
  alphabet->letters = NULL;

  return alphabet;
}

void
waft_alphabet_free (WaftAlphabet *alphabet)
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

WaftSequence*
waft_sequence_new (void)
{
  WaftSequence *seq;

  seq           = malloc (sizeof (*seq));
  seq->name     = NULL;
  seq->seq      = NULL;
  seq->alphabet = NULL;
  seq->size     = 0;

  return seq;
}

void
waft_sequence_free (WaftSequence *seq)
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
waft_sequence_to_string (WaftSequence *seq)
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
waft_sequences_iter (WaftSequence  **sequences,
                      WaftSeqIterFunc func,
                      void            *data)
{
  while (*sequences)
    {
      if (!func (*sequences, data))
        return;
      sequences++;
    }
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
