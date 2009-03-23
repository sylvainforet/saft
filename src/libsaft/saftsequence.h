/* saftsequence.h
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

#ifndef __SAFT_SEQUENCE_H__
#define __SAFT_SEQUENCE_H__

#ifdef __cplusplus
extern "C"
{
#endif

/**********/
/* Letter */
/**********/

typedef unsigned char SaftLetter;

/************/
/* Alphabet */
/************/

typedef struct _SaftAlphabet SaftAlphabet;

struct _SaftAlphabet
{
  char         *name;
  /* The letters for printing, indexed as in codes
   * The first letter is this the "unknown letter",
   * and is not really part of the alphabet */
  char         *letters;
  /* The codes indexing the letters
   * Letters of the aphabet start at 1
   * The 0 is for all the unknown codes */
  SaftLetter    codes[256];
  unsigned int  size;
};

SaftAlphabet* saft_alphabet_new  (void);

void          saft_alphabet_free (SaftAlphabet *alphabet);

/* Statically predefined alphabets */

extern SaftAlphabet SaftAlphabetDNA;
extern SaftAlphabet SaftAlphabetProtein;

/***********/
/* Segment */
/***********/

typedef struct _SaftSegment SaftSegment;

struct _SaftSegment
{
  SaftLetter  *seq;
  unsigned int size;
  SaftSegment *next;
};

SaftSegment* saft_segment_new  (void);

void         saft_segment_free (SaftSegment *segment);

/************/
/* Sequence */
/************/

typedef struct _SaftSequence SaftSequence;

struct _SaftSequence
{
  char          *name;
  SaftAlphabet  *alphabet;
  SaftLetter    *seq;
  SaftSegment   *segments;
  unsigned int   size;
};

SaftSequence* saft_sequence_new       (void);

void          saft_sequence_free      (SaftSequence *seq);

char*         saft_sequence_to_string (SaftSequence *seq);

/*************/
/* Sequences */
/*************/

typedef int (*SaftSeqIterFunc)  (SaftSequence   *seq,
                                 void           *data);

void        saft_sequences_iter (SaftSequence  **sequences,
                                 SaftSeqIterFunc func,
                                 void           *data);

#ifdef __cplusplus
}
#endif

#endif /* __SAFT_SEQUENCE_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
