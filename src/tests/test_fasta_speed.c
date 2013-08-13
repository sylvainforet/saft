/* test_fasta.c
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

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "safterror.h"
#include "saftfasta.h"


#define ITERATIONS 1000
#define INPUT      "nuc_100.fasta"


static int iter_func (SaftFasta *fasta,
                      void      *data);

int
main (int    argc,
      char **argv)
{
  struct rusage  ru_start;
  struct rusage  ru_end;
  char          *path       = INPUT;
  long int       iterations = ITERATIONS;
  long int       i;

  if (argc > 1)
    path = argv[1];
  if (argc > 2)
    iterations = atol(argv[2]);
  if (argc > 3)
    saft_error ("Usage: test_fasta_speed [FASTA [ITERATIONS]]");

  getrusage (RUSAGE_SELF, &ru_start);
  for (i = 0; i < iterations; i++)
    saft_fasta_iter (path, iter_func, NULL);
  getrusage (RUSAGE_SELF, &ru_end);

  printf ("Fasta speed results (%ld iterations on file `%s')\n",
          iterations, path);
  printf ("User  : %ld.%06lds\n",
          ru_end.ru_utime.tv_sec - ru_start.ru_utime.tv_sec,
          ru_end.ru_utime.tv_usec - ru_start.ru_utime.tv_usec);
  printf ("System: %ld.%06lds\n",
          ru_end.ru_stime.tv_sec - ru_start.ru_stime.tv_sec,
          ru_end.ru_stime.tv_usec - ru_start.ru_stime.tv_usec);

  return 0;
}

static int
iter_func (SaftFasta *fasta,
           void      *data)
{
  return 1;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
