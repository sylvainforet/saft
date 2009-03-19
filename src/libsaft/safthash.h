/* safthash.h
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

#ifndef __SAFT_HASH__H__
#define __SAFT_HASH__H__


#include <saftsequence.h>


#ifdef __cplusplus
extern "C"
{
#endif

/*************/
/* Hash Node */
/*************/

typedef struct _SaftHNode SaftHNode;

struct _SaftHNode
{
  SaftHNode     *next;
  SaftLetter    *seq;
  unsigned int   count_query;
  unsigned int   count_subject;
};


SaftHNode *saft_hnode_new  (void);

void       saft_hnode_free (SaftHNode *node);


/**************/
/* Hash Table */
/**************/

typedef struct _SaftHTable SaftHTable;

struct _SaftHTable
{
  SaftHNode    **table;
  unsigned int   shift;
  unsigned int   size;
  unsigned int   word_size;
  unsigned int   hmask;
};


SaftHTable*   saft_htable_new            (SaftAlphabet   *alphabet,
                                          unsigned int    word_size);

void          saft_htable_free           (SaftHTable     *table);

void          saft_htable_clear          (SaftHTable     *table);

void          saft_htable_add_query      (SaftHTable     *table,
                                          SaftSequence   *seq);

void          saft_htable_add_subject    (SaftHTable     *table,
                                          SaftSequence   *seq);

void          saft_htable_clear_subject  (SaftHTable     *table);

unsigned int  saft_htable_hash           (SaftHTable     *table,
                                          SaftLetter     *start);

int           saft_htable_cmp            (SaftHTable     *table,
                                          SaftHNode      *node,
                                          SaftLetter     *start);

unsigned int  saft_htable_d2             (SaftHTable     *tab);

SaftHNode*    saft_htable_lookup         (SaftHTable     *table,
                                          SaftLetter     *start);

#ifdef __cplusplus
}
#endif

#endif /* __SAFT_HASH__H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
