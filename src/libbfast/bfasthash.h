/**
 *
 */

#ifndef __BFAST_HASH__H__
#define __BFAST_HASH__H__


#include "bfastsequence.h"


#ifdef __cplusplus
extern "C"
{
#endif

/*************/
/* Hash Node */
/*************/

typedef struct _BfastHNode BfastHNode;

struct _BfastHNode
{
  BfastHNode    *next;
  unsigned char *seq;
  unsigned int   count;
};


BfastHNode *bfast_hnode_new  (void);

void        bfast_hnode_free (BfastHNode *node);


/**************/
/* Hash Table */
/**************/

typedef struct _BfastHTable BfastHTable;

struct _BfastHTable
{
  BfastHNode   **table;
  unsigned int   shift;
  unsigned int   size;
  unsigned int   word_size;
};


BfastHTable  *bfast_htable_new     (BfastAlphabet   *alphabet,
                                    unsigned int     word_size);

void          bfast_htable_free    (BfastHTable     *table);

void          bfast_htable_clear   (BfastHTable     *table);

void          bfast_htable_add_seq (BfastHTable     *table,
                                    BfastSequence   *seq);

void          bfast_htable_add     (BfastHTable     *table,
                                    BfastLetter     *start);

unsigned int  bfast_htable_hash    (BfastHTable     *table,
                                    BfastLetter     *start);

int           bfast_htable_cmp     (BfastHTable     *table,
                                    BfastHNode      *node,
                                    BfastLetter     *start);

unsigned int  bfast_htable_d2      (BfastHTable     *tab1,
                                    BfastHTable     *tab2);

BfastHNode   *bfast_htable_lookup  (BfastHTable     *table,
                                    BfastLetter     *start);

#ifdef __cplusplus
}
#endif

#endif /* __BFAST_HASH__H__ */
