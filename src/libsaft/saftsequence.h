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

/***************/
/* Nucleotides */
/***************/

#define BITS_PER_BYTE   8
#define BITS_PER_NUC    2
#define NUCS_PER_BYTE   (BITS_PER_BYTE / BITS_PER_NUC)

typedef enum
{
  NUC_A  = 0,
  NUC_C  = 1,
  NUC_G  = 2,
  NUC_T  = 3,
  NUC_NB = 4
}
Nucleotide;

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
  SaftLetter    codes[128];
  unsigned int  size;
};

SaftAlphabet* saft_alphabet_new  (void);

void          saft_alphabet_free (SaftAlphabet *alphabet);

/* Statically predefined alphabets */

extern SaftAlphabet SaftAlphabetDNA;
extern SaftAlphabet SaftAlphabetProtein;

/************/
/* Sequence */
/************/

typedef struct _SaftSequence SaftSequence;

struct _SaftSequence
{
  char   *name;
  char   *seq;

  size_t  name_length;
  size_t  seq_length;

  size_t  name_alloc;
  size_t  seq_alloc;
};

SaftSequence* saft_sequence_new  (void);

void          saft_sequence_free (SaftSequence *seq);

SaftSequence* saft_sequence_copy (SaftSequence *seq);

/* TODO add functions for six frame translations */

#ifdef __cplusplus
}
#endif

#endif /* __SAFT_SEQUENCE_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
