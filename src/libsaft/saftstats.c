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
#include <stdlib.h>

#include <gsl/gsl_cdf.h>

#include "saftstats.h"

static double saft_stats_sum_freq_pow (double      *f,
                                       unsigned int l,
                                       unsigned int freq_pow,
                                       unsigned int sum_pow);


/* FIXME There's quite a bit of optimisation and checking for numeric stability
 * that could be done here */
SaftStatsContext*
saft_stats_context_new (unsigned int word_size,
                        double      *letters_frequencies,
                        unsigned int n_letters)
{
#define p(freq_pow, sum_pow) saft_stats_sum_freq_pow(letters_frequencies, n_letters, freq_pow, sum_pow)

  SaftStatsContext *context;
  int               k = word_size;
  unsigned int      i;
  unsigned int      j;

  context           = malloc (sizeof (*context));
  context->word_size = word_size;
  context->p_2_k     = p(2, k);
  context->cov_crab  = 0;
  context->cov_diag  = 0;
  context->cov_ac1   = 0;
  context->cov_ac2   = 0;
  context->unif      = 1;

  for (i = 1; i < n_letters; i++)
    if (letters_frequencies[i] != letters_frequencies[0])
      {
        context->unif = 0;
        break;
      }

  context->sum_var_Yu = p (2, k) - p (2, 2 * k);

  if (!context->unif)
    context->cov_crab = (p(3, k) +
                         2 * p(2, 2) * p(3, 1) *
                         ((p(3, k - 1) - p(2, 2 * (k - 1))) / (p(3, 1) - p(2, 2))) -
                         (2 * k - 1) * p(2, 2 * k));

  if (context->word_size == 1)
    return context;

  context->cov_diag = (p(2, k + 1) *
                      ((1 - p(2, k - 1)) / (1 - p(2, 1))) -
                      (k - 1) * p(2, 2 * k)) * 2;

  for (i = 1; i < k; i++)
    for (j = 0; j < i; j++)
      {
        unsigned int nu = (k - j) / (i - j);
        unsigned int ro = (k - j) % (i - j);

        context->cov_ac1 += (p(2, 2 * j) *
                            p(2 * nu + 3, ro) *
                            p(2 * nu + 1, i - j - ro) -
                            p(2, 2 * k));
      }
  context->cov_ac1 *= 4;

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
          context->cov_ac2 += prod1 * prod2;
        }
    }
  context->cov_ac2 -= (k - 1) * (k - 1) * p(2, 2 * k);
  context->cov_ac2 *= 2;

  return context;

#undef p
}

void
saft_stats_context_free (SaftStatsContext *context)
{
  if (context)
    free (context);
}

double
saft_stats_mean (SaftStatsContext *context,
                 unsigned int      query_size,
                 unsigned int      subject_size)
{
  return query_size * subject_size * context->p_2_k;
}

double
saft_stats_var (SaftStatsContext *context,
                unsigned int      query_size,
                unsigned int      subject_size)
{
  double sum_var_Yu;
  double cov_crab;
  double cov_diag;
  double cov_ac1;
  double cov_ac2;
  int    m = query_size;
  int    n = subject_size;
  int    k = context->word_size;

  sum_var_Yu = m * n * context->sum_var_Yu;
  cov_crab   = m * n * (n + m - 4 * k + 2) * context->cov_crab;

  if (context->word_size == 1)
    return sum_var_Yu + cov_crab;

  cov_diag = m * n * context->cov_diag;
  cov_ac1  = m * n * context->cov_ac1;
  cov_ac2  = m * n * context->cov_ac2;

  return sum_var_Yu + cov_crab + cov_diag + cov_ac1 + cov_ac2;
}

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

double
saft_stats_pgamma_m_v (double d2,
                       double mean,
                       double var)
{
  const double scale = var / mean;
  const double shape = mean / scale;

  return saft_stats_pgamma (d2, shape, scale);
}

double
saft_stats_pgamma (double d2,
                   double shape,
                   double scale)
{
  return gsl_cdf_gamma_Q (d2, shape, scale);
}

/**
 * Benjamini and Hochberg method
 * p_values are expected to be already sorted in increasing order
 * p_values are sorted in place
 */
double*
saft_stats_BH_array (double       *p_values,
                     unsigned int  n_p_values)
{
  int i;

  /* We start at the second p-value (from the end) because the first one is
   * always unchanged
   */
  for (i = n_p_values - 2; i >= 0; i--)
    {
      p_values[i] = (p_values[i] * n_p_values) / (i + 1);
      if (p_values[i] > p_values[i + 1])
        p_values[i] = p_values[i + 1];
    }
  return p_values;
}

/**
 * This version is designed to be called on structures where p-values to avoid
 * having to allocate an array of p_values
 * Index is assumed to be zero based
 */
double
saft_stats_BH_element (double       p_value,
                       double       p_previous,
                       unsigned int index,
                       unsigned int n_p_values)
{
  double adjusted = (p_value * n_p_values) / (index + 1);

  if (adjusted > p_previous)
    return p_previous;

  return adjusted;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
