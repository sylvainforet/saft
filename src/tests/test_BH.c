/* test_BH.c
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
#include <stdio.h>

#include "saftstats.h"

int
main (int    argc,
      char **argv)
{
  double       p_values[] = {0.0001, 0.001, 0.01, 0.1, 0.5};
  unsigned int n_p_values = sizeof p_values / sizeof (*p_values);
  double      *adj_p_values;
  double       prev;
  int          i;

  adj_p_values = malloc (n_p_values * sizeof (*adj_p_values));

  for (i = 0; i < n_p_values; i++)
    adj_p_values[i] = p_values[i];

  printf ("Unadjusted:\n");
  for (i = 0; i < n_p_values; i++)
    printf ("%.4f ", p_values[i]);
  printf ("\n");

  adj_p_values = saft_stats_BH_array (adj_p_values,
                                      n_p_values);

  printf ("Adjusted:\n");
  for (i = 0; i < n_p_values; i++)
    printf ("%.4f ", adj_p_values[i]);
  printf ("\n");

  for (i = 0; i < n_p_values; i++)
    adj_p_values[i] = p_values[i];

  printf ("\n");
  printf ("Unadjusted:\n");
  for (i = 0; i < n_p_values; i++)
    printf ("%.4f ", adj_p_values[i]);
  printf ("\n");

  prev = 1;
  for (i = n_p_values - 1; i >= 0; i--)
    {
      adj_p_values[i] = saft_stats_BH_element (adj_p_values[i],
                                               prev,
                                               i,
                                               n_p_values);
      prev = adj_p_values[i];
    }

  printf ("Adjusted:\n");
  for (i = 0; i < n_p_values; i++)
    printf ("%.4f ", adj_p_values[i]);
  printf ("\n");

  free (adj_p_values);

  return 0;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
