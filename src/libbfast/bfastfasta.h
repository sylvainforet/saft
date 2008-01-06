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
