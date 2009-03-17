/* wafthash.h
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

#ifndef __WAFT_HASH__H__
#define __WAFT_HASH__H__


#include "waftsequence.h"


#ifdef __cplusplus
extern "C"
{
#endif

/*************/
/* Hash Node */
/*************/

typedef struct _WaftHNode WaftHNode;

struct _WaftHNode
{
  WaftHNode    *next;
  unsigned char *seq;
  unsigned int   count;
};


WaftHNode *waft_hnode_new  (void);

void        waft_hnode_free (WaftHNode *node);


/**************/
/* Hash Table */
/**************/

typedef struct _WaftHTable WaftHTable;

struct _WaftHTable
{
  WaftHNode   **table;
  unsigned int   shift;
  unsigned int   size;
  unsigned int   word_size;
};


WaftHTable  *waft_htable_new     (WaftAlphabet   *alphabet,
                                    unsigned int     word_size);

void          waft_htable_free    (WaftHTable     *table);

void          waft_htable_clear   (WaftHTable     *table);

void          waft_htable_add_seq (WaftHTable     *table,
                                    WaftSequence   *seq);

void          waft_htable_add     (WaftHTable     *table,
                                    WaftLetter     *start);

unsigned int  waft_htable_hash    (WaftHTable     *table,
                                    WaftLetter     *start);

int           waft_htable_cmp     (WaftHTable     *table,
                                    WaftHNode      *node,
                                    WaftLetter     *start);

unsigned int  waft_htable_d2      (WaftHTable     *tab1,
                                    WaftHTable     *tab2);

WaftHNode   *waft_htable_lookup  (WaftHTable     *table,
                                    WaftLetter     *start);

#ifdef __cplusplus
}
#endif

#endif /* __WAFT_HASH__H__ */
