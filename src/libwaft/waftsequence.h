/* waftsequence.h
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

#ifndef __WAFT_SEQUENCE_H__
#define __WAFT_SEQUENCE_H__

#ifdef __cplusplus
extern "C"
{
#endif

/**********/
/* Letter */
/**********/

typedef unsigned char WaftLetter;

/************/
/* Alphabet */
/************/

typedef struct _WaftAlphabet WaftAlphabet;

struct _WaftAlphabet
{
  char         *name;
  /* The letters for printing, indexed as in codes
   * The first letter is this the "unknown letter",
   * and is not really part of the alphabet */
  char         *letters;
  /* The codes indexing the letters
   * Letters of the aphabet start at 1
   * The 0 is for all the unknown codes */
  WaftLetter    codes[256];
  unsigned int  size;
};

WaftAlphabet* waft_alphabet_new  (void);

void          waft_alphabet_free (WaftAlphabet *alphabet);

/* Statically predefined alphabets */

extern WaftAlphabet WaftAlphabetDNA;
extern WaftAlphabet WaftAlphabetProtein;

/************/
/* Sequence */
/************/

typedef struct _WaftSequence WaftSequence;

struct _WaftSequence
{
  char          *name;
  WaftLetter    *seq;
  WaftAlphabet  *alphabet;
  unsigned int   size;
};

WaftSequence* waft_sequence_new       (void);

void          waft_sequence_free      (WaftSequence *seq);

char*         waft_sequence_to_string (WaftSequence *seq);

/*************/
/* Sequences */
/*************/

typedef int (*WaftSeqIterFunc)  (WaftSequence   *seq,
                                 void           *data);

void        waft_sequences_iter (WaftSequence  **sequences,
                                 WaftSeqIterFunc func,
                                 void           *data);

#ifdef __cplusplus
}
#endif

#endif /* __WAFT_SEQUENCE_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
