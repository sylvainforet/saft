/* bfastsequence.h
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

#ifndef __BFAST_SEQUENCE_H__
#define __BFAST_SEQUENCE_H__

#ifdef __cplusplus
extern "C"
{
#endif

/**********/
/* Letter */
/**********/

typedef unsigned char BfastLetter;

/************/
/* Alphabet */
/************/

typedef struct _BfastAlphabet BfastAlphabet;

struct _BfastAlphabet
{
  char         *name;
  /* The letters for printing, indexed as in codes
   * The first letter is this the "unknown letter",
   * and is not really part of the alphabet */
  char         *letters;
  /* The codes indexing the letters
   * Letters of the aphabet start at 1
   * The 0 is for all the unknown codes */
  BfastLetter   codes[256];
  unsigned int  size;
};

BfastAlphabet* bfast_alphabet_new (void);

void           bfast_alphabet_free (BfastAlphabet *alphabet);

/* Statically predefined alphabets */

extern BfastAlphabet BfastAlphabetDNA;
extern BfastAlphabet BfastAlphabetProtein;

/************/
/* Sequence */
/************/

typedef struct _BfastSequence BfastSequence;

struct _BfastSequence
{
  char          *name;
  BfastLetter   *seq;
  BfastAlphabet *alphabet;
  unsigned int   size;
};

BfastSequence* bfast_sequence_new       (void);

void           bfast_sequence_free      (BfastSequence   *seq);

char*          bfast_sequence_to_string (BfastSequence   *seq);

/*************/
/* Sequences */
/*************/

typedef int    (*BfastSeqIterFunc)      (BfastSequence   *seq,
                                         void            *data);

void           bfast_sequences_iter     (BfastSequence  **sequences,
                                         BfastSeqIterFunc func,
                                         void            *data);

#ifdef __cplusplus
}
#endif

#endif /* __BFAST_SEQUENCE_H__ */
