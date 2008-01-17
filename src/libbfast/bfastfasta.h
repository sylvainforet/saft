/* bfastfasta.h
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

#ifndef __BFAST_FASTA_H__
#define __BFAST_FASTA_H__

#include "bfastsequence.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* A sequence read from a fasta formated file */

typedef struct _BfastFasta BfastFasta;

struct _BfastFasta
{
  char *name;
  char *seq;
};

typedef int    (*BfastFastaIterFunc) (BfastFasta        *fasta,
                                      void              *data);

BfastFasta*    bfast_fasta_new       (void);

void           bfast_fasta_free      (BfastFasta        *fasta);

BfastFasta**   bfast_fasta_read      (const char        *filename,
                                      unsigned int      *n);

void           bfast_fasta_iter      (const char        *filename,
                                      BfastFastaIterFunc func,
                                      void              *data);

BfastSequence* bfast_fasta_to_seq    (BfastFasta        *fasta,
                                      BfastAlphabet     *alphabet);

#ifdef __cplusplus
}
#endif

#endif /* __BFAST_FASTA_H__ */
