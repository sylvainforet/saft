/* saftmain.c
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


#define _GNU_SOURCE
#include <errno.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "safterror.h"
#include "saftfasta.h"
#include "saftsearch.h"


typedef struct _SaftOptDesc SaftOptDesc;

struct _SaftOptDesc
{
  char *name;
  int   has_arg;
  int   val;
  char *description;
};

static SaftOptDesc opt_desc[] =
{
    /* Program information */
    {"help",        no_argument,       'h', "Prints a short help"},
    {"version",     no_argument,       'V', "Prints the program's version"},
    /* General options */
    {"verbose",     no_argument,       'v', "Increases the program's verbosity"},
    /* TODO multi threading options ? */
    /* Input / Output */
    {"input",       required_argument, 'i', "Path to the input file"},
    {"database",    required_argument, 'd', "Path to the database to search"},
    {"output",      required_argument, 'o', "Path to the output file"},
    /* Search setup */
    {"program",     required_argument, 'p', "Program to use"},
    {"wordsize",    required_argument, 'w', "Word size"},
    {"showmax",     required_argument, 'b', "Maximum number of results to show"},
    {"pmax",        required_argument, 'e', "Show results with a p-value smaller than this"},
    {"letter_freq", required_argument, 'f', "Comma separated list of letter frequencies"},
    /* FIXME add option to choose the strand(s) of the query */

    /* TODO We could have an extra mechanism to add engine specific options */
    /* The two options below could be implemented with that mechanism */
    /* This would be analogous to the -o option of 'mount' */
    {"cacheq",      no_argument,       'q', "Cache the queries in memory (only used by some engines)"},
    {"cached",      no_argument,       'a', "Cache the database in memory (only used by some engines)"},

    {NULL, 0, 0, NULL}
};

static struct option*   saft_main_get_options   (void);

static char*            saft_main_get_optstring (void);

static void             saft_main_usage         (char        *argv0);

static void             saft_main_help          (char        *argv0);

static void             saft_main_version       (void);

static SaftProgramType  saft_main_program_type  (char        *program);

static int              saft_main_search        (SaftOptions *options);

static void             saft_main_write_search  (SaftOptions *options,
                                                 SaftSearch  *search,
                                                 FILE        *stream);

int
main (int    argc,
      char **argv)
{
  struct option  *long_options = saft_main_get_options ();
  SaftOptions    *options      = saft_options_new ();
  int             ret          = 0;
  char           *tmp_freqs    = NULL;
  char           *optstring;
  char           *endptr;

  optstring = saft_main_get_optstring ();
  while (1)
    {
      int option_index = 0;
      int c;

      c = getopt_long (argc,
                       argv,
                       optstring,
                       long_options,
                       &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'h':
              saft_main_help (argv[0]);
              goto cleanup;
          case 'V':
              saft_main_version ();
              goto cleanup;
          case 'v':
              options->verbosity++;
              break;
          case 'i':
              options->input_path = strdup (optarg);
              break;
          case 'd':
              options->db_path = strdup (optarg);
              break;
          case 'o':
              options->output_path = strdup (optarg);
              break;
          case 'p':
              options->program = saft_main_program_type (optarg);
              if (options->program == SAFT_UNKNOWN_PROGRAM)
                {
                  ret = 1;
                  saft_error ("Wrong `--program (-p)' argument: unknown program `%s'", optarg);
                  goto cleanup;
                }
              break;
          case 'w':
              options->word_size = strtol (optarg, &endptr, 10);
              if (errno == ERANGE || errno == EINVAL || *endptr != '\0')
                {
                  saft_error ("Wrong `--wordsize (-w)' argument: could not convert `%s' to an integer", optarg);
                  ret = 1;
                  goto cleanup;
                }
              break;
          case 'b':
              options->max_results = strtol (optarg, &endptr, 10);
              if (errno == ERANGE || errno == EINVAL || *endptr != '\0')
                {
                  saft_error ("Wrong `--showmax (-b)' argument: could not convert `%s' to an integer", optarg);
                  ret = 1;
                  goto cleanup;
                }
              break;
          case 'e':
              options->p_max = strtod (optarg, &endptr);
              if (errno == ERANGE || *endptr != '\0')
                {
                  saft_error ("Wrong `--pmax (-e)' argument: could not convert `%s' to a double", optarg);
                  ret = 1;
                  goto cleanup;
                }
              break;
          case 'f':
              tmp_freqs = optarg;
              break;
          case 'q':
              options->cache_queries = 1;
              break;
          case 'a':
              options->cache_db = 1;
              break;
          default:
              saft_main_usage (argv[0]);
              ret = 1;
              goto cleanup;
        }
    }
  if (options->program == SAFT_UNKNOWN_PROGRAM)
    {
      saft_error ("Program type was not provided, use the `--program' or `-p' option");
      ret = 1;
      goto cleanup;
    }
  else
    {
      switch (options->program)
        {
          case SAFTN:
              options->alphabet = &SaftAlphabetDNA;
              break;
          case SAFTP:
          case SAFTX:
          case TSAFTN:
          case TSAFTX:
              options->alphabet = &SaftAlphabetProtein;
              break;
          case SAFT:
              saft_error ("Generic SAFT program not implemented");
              break;
          default:
              saft_error ("Unknown SAFT program");
              break;
        }
    }
  if (!options->input_path)
    {
      saft_error ("Query file was not provided, use the `--input' or `-i' option'");
      ret = 1;
      goto cleanup;
    }
  if (!options->db_path)
    {
      saft_error ("Database file was not provided, use the `--database' or `-d' option'");
      ret = 1;
      goto cleanup;
    }
  if (options->word_size == 0)
    {
      if (options->program == SAFTN)
        options->word_size = 7;
      else
        options->word_size = 3;
    }
  if (options->cache_queries && options->cache_db)
    {
      saft_error ("Can't cache both the queries (-q) and the database (-a)");
    }

  options->letter_frequencies = malloc (options->alphabet->size * sizeof (*options->letter_frequencies));
  if (tmp_freqs)
    {
      const double epsilon = 1e-6;
      double       tot = 0.0;
      char        *saveptr;
      int          i;

      for (i = 0; i < options->alphabet->size; i++)
        {
          char *token;

          token = strtok_r (tmp_freqs, ",", &saveptr);
          if (!token)
            {
              saft_error ("Failed to parse frequencies: `%s'", tmp_freqs);
              ret = 1;
              goto cleanup;
            }

          options->letter_frequencies[i] = strtod (token, &endptr);
          if (errno == ERANGE || *endptr != '\0')
            {
              saft_error ("Failed to convert frequency `%s' to double", token);
              ret = 1;
              goto cleanup;
            }
          tot += options->letter_frequencies[i];
        }

      if (tot + epsilon < 1.0 || tot - epsilon > 1.0)
        {
          saft_error ("Frequencies did not sum to 1.0: `%s'", tmp_freqs);
          ret = 1;
          goto cleanup;
        }
    }
  else
    {
      const double f = 1. / options->alphabet->size;
      int          i;

      for (i = 0; i < options->alphabet->size; i++)
        options->letter_frequencies[i] = f;
    }

  saft_main_search (options);

cleanup:
  free (optstring);
  free (long_options);
  saft_options_free (options);

  return ret;
}

static struct option*
saft_main_get_options ()
{
  struct option *long_options;
  int            nb_options;
  int            i;

  for (nb_options = 0; opt_desc[nb_options].name; nb_options++);

  long_options = malloc ((nb_options + 1) * sizeof (*long_options));

  for (i = 0; i < nb_options; i++)
    {
      long_options[i].name    = opt_desc[i].name;
      long_options[i].has_arg = opt_desc[i].has_arg;
      long_options[i].flag    = NULL;
      long_options[i].val     = opt_desc[i].val;
    }
  long_options[i].name    = NULL;
  long_options[i].has_arg = 0;
  long_options[i].flag    = NULL;
  long_options[i].val     = 0;

  return long_options;
}

static char*
saft_main_get_optstring ()
{
  char *optstring;
  int   nb_options;
  int   i;
  int   j;

  for (nb_options = 0; opt_desc[nb_options].name; nb_options++);

  optstring = malloc (2 * nb_options + 1);

  for (i = 0, j = 0; i < nb_options; i++)
    {
      optstring[j] = (char)opt_desc[i].val;
      j++;
      if (opt_desc[i].has_arg == required_argument)
        {
          optstring[j] = ':';
          j++;
        }
    }
  optstring[j] = '\0';

  return optstring;
}

static void
saft_main_usage (char *argv0)
{
  char *prog = basename (argv0);
  printf ("Usage: %s OPTIONS\n", prog);
  printf ("Try %s -h for help\n", prog);
}

static void
saft_main_help (char *argv0)
{
  char        *prog    = basename (argv0);
  size_t       longest = 0;
  unsigned int i;

  for (i = 0; opt_desc[i].name; i++)
    {
      const size_t l = strlen (opt_desc[i].name);
      if (l > longest)
        longest = l;
    }

  printf ("SAFT (Sequence Alignment Free Tool)\n");
  printf ("Usage: %s OPTIONS\n", prog);
  printf ("Where OPTIONS are:\n");

  for (i = 0; opt_desc[i].name; i++)
    printf ("  --%-*s (-%c) : %s\n",
            (int)longest,
            opt_desc[i].name,
            opt_desc[i].val,
            opt_desc[i].description);
}

static void
saft_main_version (void)
{
  printf ("SAFT version " SAFT_VERSION "\n");
}

static SaftProgramType
saft_main_program_type (char *program)
{
  unsigned int i;
  for (i = 0; i < NB_SAFT_PROGRAMS; i++)
    if (!strcmp (program, saft_program_names[i]))
      return i;
  return SAFT_UNKNOWN_PROGRAM;
}

static int
saft_main_search (SaftOptions *options)
{
  SaftSearchEngine *engine;
  SaftSearch       *searches;
  SaftSearch       *search;
  FILE             *out_stream;

  if (options->output_path == NULL)
    out_stream = stdout;
  else if ((out_stream = fopen (options->output_path, "w")) == NULL) 
    {
      saft_error ("Could not open output file");
      return 1;
    }

  engine = saft_search_engine_new (options);
  if (!engine)
    {
      saft_error ("Could not create search engine");
      return 1;
    }

  searches = saft_search_all (engine, options->input_path, options->db_path);
  saft_search_engine_free (engine);

  for (search = searches; search; search = search->next)
    {
      saft_main_write_search (options,
                              search,
                              out_stream);
    }
  saft_search_free_all (search);

  if (options->output_path == NULL)
    fclose (out_stream);
  return 0;
}

static void
saft_main_write_search (SaftOptions *options,
                        SaftSearch  *search,
                        FILE        *stream)
{
  int i;

  fprintf (stream, "Query: %s program: %s word size: %ld\n",
           search->name,
           saft_program_names[options->program],
           (long)options->word_size);

  saft_search_adjust_pvalues (search);
  for (i = 0; i < search->n_results; i++)
    {
      fprintf (stream, "  Hit: %s D2: %ld adj.p.val: %.5e p.val: %.5e\n",
               search->results[i]->name,
               search->results[i]->d2,
               search->results[i]->p_value_adj,
               search->results[i]->p_value);
    }
  if (i == 0)
    fprintf (stream, "No hit found\n");
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
