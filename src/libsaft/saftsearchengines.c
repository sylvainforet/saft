/* saftsearchengines.c
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


#include <stdlib.h>
#include <string.h>

#include "saftsearchengines.h"


/*************************/
/* Generic Search Engine */
/*************************/

typedef struct _SearchEngineGeneric SearchEngineGeneric;

struct _SearchEngineGeneric
{
  SaftSearchEngine search_engine;
};


static void        search_engine_generic_free                 (SaftSearchEngine *engine);

static SaftSearch* search_engine_generic_search_two_sequences (SaftSearchEngine *engine,
                                                               SaftSequence     *query,
                                                               SaftSequence     *subject);

static SaftSearch* search_engine_generic_search_all           (SaftSearchEngine *engine,
                                                               const char       *query_path,
                                                               const char       *db_path);

SaftSearchEngine*
saft_search_engine_generic (SaftOptions *options)
{
  SearchEngineGeneric *engine;

  engine                                     = malloc (sizeof (*engine));
  engine->search_engine.options              = options;
  engine->search_engine.search_two_sequences = search_engine_generic_search_two_sequences;
  engine->search_engine.search_all           = search_engine_generic_search_all;
  engine->search_engine.free                 = search_engine_generic_free;

  return (SaftSearchEngine*)engine;
}

static void
search_engine_generic_free (SaftSearchEngine *engine)
{
}

static SaftSearch*
search_engine_generic_search_two_sequences (SaftSearchEngine *engine,
                                            SaftSequence     *query,
                                            SaftSequence     *subject)
{
  SearchEngineGeneric *se;
  
  se = (SearchEngineGeneric*) engine;
  se++;

  return NULL;
}

static SaftSearch*
search_engine_generic_search_all (SaftSearchEngine *engine,
                                  const char       *query_path,
                                  const char       *db_path)
{
  SearchEngineGeneric *se;
  
  se = (SearchEngineGeneric*) engine;
  se++;

  return NULL;
}

/*********************************/
/* Array-based DNA Search Engine */
/*********************************/

typedef struct _SearchEngineDNAArray SearchEngineDNAArray;

struct _SearchEngineDNAArray
{
  SaftSearchEngine search_engine;
};


static void        search_engine_dna_array_free                 (SaftSearchEngine *engine);

static SaftSearch* search_engine_dna_array_search_two_sequences (SaftSearchEngine *engine,
                                                                 SaftSequence     *query,
                                                                 SaftSequence     *subject);

static SaftSearch* search_engine_dna_array_search_all           (SaftSearchEngine *engine,
                                                                 const char       *query_path,
                                                                 const char       *db_path);

SaftSearchEngine*
saft_search_engine_dna_array_new (SaftOptions *options)
{
  SearchEngineDNAArray *engine;

  engine                                     = malloc (sizeof (*engine));
  engine->search_engine.options              = options;
  engine->search_engine.search_two_sequences = search_engine_dna_array_search_two_sequences;
  engine->search_engine.search_all           = search_engine_dna_array_search_all;
  engine->search_engine.free                 = search_engine_dna_array_free;

  return (SaftSearchEngine*)engine;
}

static void
search_engine_dna_array_free (SaftSearchEngine *engine)
{
}

static SaftSearch*
search_engine_dna_array_search_two_sequences (SaftSearchEngine *engine,
                                              SaftSequence     *query,
                                              SaftSequence     *subject)
{
  SearchEngineDNAArray *se;

  se = (SearchEngineDNAArray*) engine;
  se++;

  return NULL;
}

static SaftSearch*
search_engine_dna_array_search_all (SaftSearchEngine *engine,
                                    const char       *query_path,
                                    const char       *db_path)
{
  SearchEngineDNAArray *se;

  se = (SearchEngineDNAArray*) engine;
  se++;

  return NULL;
}

/********************************/
/* Hash-based DNA Search Engine */
/********************************/

typedef struct _SearchEngineDNAHash SearchEngineDNAHash;

struct _SearchEngineDNAHash
{
  SaftSearchEngine search_engine;
};


static void        search_engine_dna_hash_free                 (SaftSearchEngine *engine);

static SaftSearch* search_engine_dna_hash_search_two_sequences (SaftSearchEngine *engine,
                                                                SaftSequence     *query,
                                                                SaftSequence     *subject);

static SaftSearch* search_engine_dna_hash_search_all           (SaftSearchEngine *engine,
                                                                const char       *query_path,
                                                                const char       *db_path);

SaftSearchEngine*
saft_search_engine_dna_hash_new (SaftOptions *options)
{
  SearchEngineDNAHash *engine;

  engine                                     = malloc (sizeof (*engine));
  engine->search_engine.options              = options;
  engine->search_engine.search_two_sequences = search_engine_dna_hash_search_two_sequences;
  engine->search_engine.search_all           = search_engine_dna_hash_search_all;
  engine->search_engine.free                 = search_engine_dna_hash_free;

  return (SaftSearchEngine*)engine;
}

static void
search_engine_dna_hash_free (SaftSearchEngine *engine)
{
}

static SaftSearch*
search_engine_dna_hash_search_two_sequences (SaftSearchEngine *engine,
                                             SaftSequence     *query,
                                             SaftSequence     *subject)
{
  SearchEngineDNAHash *se;

  se = (SearchEngineDNAHash*) engine;
  se++;

  return NULL;
}

static SaftSearch*
search_engine_dna_hash_search_all (SaftSearchEngine *engine,
                                   const char       *query_path,
                                   const char       *db_path)
{
  SearchEngineDNAHash *se;

  se = (SearchEngineDNAHash*) engine;
  se++;

  return NULL;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
