/* saftsearch.h
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

#ifndef __SAFT_SEARCH_H__
#define __SAFT_SEARCH_H__


#include "saftsequence.h"


#ifdef __cplusplus
extern "C"
{
#endif

typedef enum
{
  /* Frequencies are calculated based on the compostion of the query */
  SAFT_FREQ_QUERY = 0,
  /* Frequencies are calculated based on the compostion of the database */
  SAFT_FREQ_SUBJECTS,
  /* Frequencies are calculated based on the compostion of the database and the
   * query */
  SAFT_FREQ_QUERY_SUBJECTS,
  /* Frequencies are provided by the user */
  SAFT_FREQ_USER,
  /* Uniform frequencies across the alphabet */
  SAFT_FREQ_UNIFORM,
  NB_SAFT_FREQUENCIES,
  SAFT_UNKNOWN_FREQUENCY
}
SaftFreqType;

typedef enum
{
  /* Generic saft program type to deal with any possible alphabet. In this case
   * alphabet must be specified */
  SAFT = 0,
  /* Nucleotide query against nucleotide database */
  SAFTN,
  /* Protein query against protein database */
  SAFTP,
  /* Six frame translations of a nucleotide query against protein database */
  SAFTX,
  /* Protein query against six frame translation of a nucleotide database */
  TSAFTN,
  /* Six frame translations of a nucleotide query against six frame translation
   * of a nucleotide database */
  TSAFTX,
  NB_SAFT_PROGRAMS,
  SAFT_UNKNOWN_PROGRAM
}
SaftProgramType;

extern const char *saft_program_names[NB_SAFT_PROGRAMS];

/* TODO add an option for the type of p-value approximation.
 * Gamma is generally better, but a normal approximation might be desirable if
 * the sequences are large enough and speed is important */

/******************/
/* Search Options */
/******************/

typedef struct _SaftOptions SaftOptions;

struct _SaftOptions
{
  char           *input_path;
  char           *db_path;
  char           *output_path;

  SaftAlphabet   *alphabet;

  double         *letter_frequencies;
  double          p_max;

  unsigned int    word_size;

  unsigned int    verbosity;
  int             show_max;

  SaftProgramType program;
  SaftFreqType    freq_type;

  unsigned int    cache_db: 1;
  unsigned int    cache_queries: 1;
};

SaftOptions* saft_options_new  (void);

void         saft_options_free (SaftOptions *options);

/**********/
/* Result */
/**********/

/* The result of a pairwise comparison */

typedef struct _SaftResult SaftResult;

struct _SaftResult
{
  SaftResult  *next;
  char        *name;
  double       p_value;
  double       p_value_adj;
  unsigned int d2;
  unsigned int subject_size;
  char         frame;
};

SaftResult* saft_result_new  (void);

void        saft_result_free (SaftResult *result);

/**********/
/* Search */
/**********/

/* The results of scanning a sequence against a database */

typedef struct _SaftSearch SaftSearch;

struct _SaftSearch
{
  SaftSearch    *next;
  SaftSequence  *query;
  SaftResult    *results;
  SaftResult   **sorted_results;
  unsigned int   n_results;
};

SaftSearch* saft_search_new             (void);

void        saft_search_free            (SaftSearch   *search);

void        saft_search_compute_pvalues (SaftSearch   *search,
                                         SaftOptions  *options);

/********************/
/* SaftSearchEngine */
/********************/

typedef struct _SaftSearchEngine SaftSearchEngine;

struct _SaftSearchEngine
{
  /* Member variables */
  SaftOptions  *options;

  /* Virtual methods table */
  SaftSearch* (*search_two_sequences) (SaftSearchEngine *engine,
                                       SaftSequence     *query,
                                       SaftSequence     *subject);
  SaftSearch* (*search_all)           (SaftSearchEngine *engine,
                                       const char       *queries_path,
                                       const char       *db_path);
  void        (*free)                 (SaftSearchEngine *engine);

};

SaftSearchEngine* saft_search_engine_new    (SaftOptions      *options);

void              saft_search_engine_free   (SaftSearchEngine *engine);

SaftSearch*       saft_search_two_sequences (SaftSearchEngine *engine,
                                             SaftSequence     *query,
                                             SaftSequence     *subject);

SaftSearch*       saft_search_all           (SaftSearchEngine *engine,
                                             const char       *query_path,
                                             const char       *db_path);

#ifdef __cplusplus
}
#endif

#endif /* __SAFT_SEARCH_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
