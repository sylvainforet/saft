/* test_pgamma.c
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

#include "saftstats.h"


int
main (int    argc,
      char **argv)
{
  unsigned int seq_sizes[]     = {10, 20, 40, 60, 100, 200, 400, 800};
  unsigned int nb_seq_sizes    = sizeof seq_sizes / sizeof (*seq_sizes);
  unsigned int word_sizes[]    = {1, 2, 4, 6, 8, 10};
  unsigned int nb_word_sizes   = sizeof word_sizes / sizeof (*word_sizes);
  /*double       letters_freqs[] = {1./4, 1./4, 1./4, 1./4};*/
  double       letters_freqs[] = {1./3, 1./3, 1./6, 1./6};
  unsigned int alphabet_size   = sizeof letters_freqs / sizeof (*letters_freqs);
  unsigned int i;
  unsigned int j;

  for (i = 0; i < nb_seq_sizes; i++)
    {
      unsigned int seq_size = seq_sizes[i];

      for (j = 0; j < nb_word_sizes; j++)
        {
          unsigned int      word_size = word_sizes[j];
          SaftStatsContext *context   = saft_stats_context_new (word_size,
                                                                letters_freqs,
                                                                alphabet_size);
          double            mean      = saft_stats_mean (context,
                                                         seq_size,
                                                         seq_size);
          double            var       = saft_stats_var  (context,
                                                         seq_size,
                                                         seq_size);
          double            pgamma    = saft_stats_pgamma_m_v (mean, mean, var);
          printf ("n = m = %-4d ; k = %-3d ; mean = %.5e ; var = %.5e ; pgamma (mean) = %.5e\n",
                  seq_size, word_size, mean, var, pgamma);
        }
    }

  return 0;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
