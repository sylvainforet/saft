/* waftsearch.h
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

#ifndef __WAFT_SEARCH_H__
#define __WAFT_SEARCH_H__


#include "wafthash.h"


#ifdef __cplusplus
extern "C"
{
#endif

/******************/
/* Search Options */
/******************/

typedef struct _WaftSearchOptions WaftSearchOptions;

struct _WaftSearchOptions
{
  unsigned int word_size;
};

WaftSearchOptions *waft_search_options_new  (void);

void                waft_search_options_free (WaftSearchOptions *options);

/**********/
/* Search */
/**********/

typedef struct _WaftSearch WaftSearch;

struct _WaftSearch
{
  WaftSearchOptions *options;
  WaftSequence      *query;
  WaftSequence      *subject;
  WaftHTable        *query_h;
  WaftHTable        *subject_h;
  unsigned long int   d2;
};

WaftSearch *waft_search_new     (void);

void         waft_search_free    (WaftSearch *search);

void         waft_search_process (WaftSearch *search);

#ifdef __cplusplus
}
#endif

#endif /* __WAFT_SEARCH_H__ */
