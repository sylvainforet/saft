/* waftfasta.h
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

#ifndef __WAFT_FASTA_H__
#define __WAFT_FASTA_H__

#include "waftsequence.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* A sequence read from a fasta formated file */

typedef struct _WaftFasta WaftFasta;

struct _WaftFasta
{
  char *name;
  char *seq;
};

typedef int    (*WaftFastaIterFunc) (WaftFasta        *fasta,
                                      void              *data);

WaftFasta*    waft_fasta_new       (void);

void           waft_fasta_free      (WaftFasta        *fasta);

WaftFasta**   waft_fasta_read      (const char        *filename,
                                      unsigned int      *n);

void           waft_fasta_iter      (const char        *filename,
                                      WaftFastaIterFunc func,
                                      void              *data);

WaftSequence* waft_fasta_to_seq    (WaftFasta        *fasta,
                                      WaftAlphabet     *alphabet);

#ifdef __cplusplus
}
#endif

#endif /* __WAFT_FASTA_H__ */
