/* saftstats.c
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

#include <math.h>

#include "saftstats.h"

static double saft_stats_sum_freq_pow (double      *f,
                                       unsigned int l,
                                       unsigned int freq_pow,
                                       unsigned int sum_pow);


double
saft_stats_mean (unsigned int m,
                 unsigned int n,
                 unsigned int k,
                 double      *f,
                 unsigned int l)
{
  double mean;

  mean = m * n * saft_stats_sum_freq_pow (f, l, 2, k);

  return mean;
}

double
saft_stats_var (unsigned int m,
                unsigned int n,
                unsigned int k,
                double      *f,
                unsigned int l)
{
#define p(freq_pow, sum_pow) saft_stats_sum_freq_pow(f, l, freq_pow, sum_pow)
  double       sum_var_Yu;
  double       cov_crab;
  double       cov_diag;
  double       cov_ac1;
  double       cov_ac2;
  double       f0 = f[0];
  unsigned int unif;
  unsigned int i;
  unsigned int j;

  unif = 1;
  for (i = 1; i < l; i++)
    if (f[i] != f0)
      {
        unif = 0;
        break;
      }

  sum_var_Yu = m * n * (p (2, k) - p (2, 2 * k));

  cov_crab = 0;
  if (!unif)
    {
      int sm = m;
      int sn = n;
      int sk = k;
      cov_crab = sm * sn * (2 + sm + sn - 4 * sk) *
          (p(3, k) +
           2 * p(2, 2) * p(3, 1) *
           ((p(3, k - 1) - p(2, 2 * (k - 1))) / (p(3, 1) - p(2, 2))) -
           (2 * sk - 1) * p(2, 2 * k));
    }

  if (k == 1)
    return sum_var_Yu + cov_crab;

  cov_diag = 2 * m * n * (p(2, k + 1) * ((1 - p(2, k - 1)) / (1 - p(2, 1))) - (k - 1) * p(2, 2 * k));

  cov_ac1 = 0;
  for (i = 1; i < k; i++)
    for (j = 0; j < i; j++)
      {
        unsigned int nu = (k - j) / (i - j);
        unsigned int ro = (k - j) % (i - j);

        cov_ac1 += (p(2, 2 * j) * p(2 * nu + 3, ro) * p(2 * nu + 1, i - j - ro) - p(2, 2 * k));
      }
  cov_ac1 *= 4 * m * n;

  cov_ac2 = 0;
  for (i = 1; i < k; i++)
    {
      for (j = 1; j < k; j++)
        {
          unsigned int x;
          unsigned int nu    = k / (i + j);
          unsigned int ro    = k % (i + j);
          double       prod1 = 1;
          double       prod2 = 1;

          for (x = 1; x <= j; x++)
            {
              unsigned int t = 1 + 2 * nu;

              if (x <= ro)
                t++;
              if (x + i <= ro)
                t++;
              prod1 *= p(t, 1);
            }
          for (x = 1; x <= i; x++)
            {
              unsigned int t = 1 + 2 * nu;

              if (x <= ro)
                t++;
              if (x + j <= ro)
                t++;
              prod2 *= p(t, 1);
            }
          cov_ac2 += prod1 * prod2;
        }
    }
  cov_ac2 -= (k - 1) * (k - 1) * p(2, 2 * k);
  cov_ac2 *= 2. * m * n;

  return sum_var_Yu + cov_crab + cov_diag + cov_ac1 + cov_ac2;
#undef p
}

/* FIXME These should be computed once and for all for a given search */
static double
saft_stats_sum_freq_pow (double      *f,
                         unsigned int l,
                         unsigned int freq_pow,
                         unsigned int sum_pow)
{
  double       res = 0;
  unsigned int i;

  for (i = 0; i < l; i++)
    res += pow (f[i], (double)freq_pow);

  res = pow (res, (double)sum_pow);

  return res;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
