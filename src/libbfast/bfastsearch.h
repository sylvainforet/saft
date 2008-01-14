/**
 *
 */

#ifndef __BFAST_SEARCH_H__
#define __BFAST_SEARCH_H__


#include "bfasthash.h"


#ifdef __cplusplus
extern "C"
{
#endif

/******************/
/* Search Options */
/******************/

typedef struct _BfastSearchOptions BfastSearchOptions;

struct _BfastSearchOptions
{
  unsigned int word_size;
};

BfastSearchOptions *bfast_search_options_new  (void);

void                bfast_search_options_free (BfastSearchOptions *options);

/**********/
/* Search */
/**********/

typedef struct _BfastSearch BfastSearch;

struct _BfastSearch
{
  BfastSearchOptions *options;
  BfastSequence      *query;
  BfastSequence      *subject;
  BfastHTable        *query_h;
  BfastHTable        *subject_h;
  unsigned long int   d2;
};

BfastSearch *bfast_search_new     (void);

void         bfast_search_free    (BfastSearch *search);

void         bfast_search_process (BfastSearch *search);

#ifdef __cplusplus
}
#endif

#endif /* __BFAST_SEARCH_H__ */
