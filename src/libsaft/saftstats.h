/* saftstats.h
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

#ifndef __SAFT_STATS_H__
#define __SAFT_STATS_H__

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct _SaftStatsContext SaftStatsContext;

struct _SaftStatsContext
{
  double p_2_k;

  double sum_var_Yu;
  double cov_crab;
  double cov_diag;
  double cov_ac1;
  double cov_ac2;

  unsigned int word_size;

  unsigned int unif: 1;
};

SaftStatsContext* saft_stats_context_new  (unsigned int      word_size,
                                           double           *letters_frequencies,
                                           unsigned int      n_letters);

void              saft_stats_context_free (SaftStatsContext *context);

double            saft_stats_mean         (SaftStatsContext *context,
                                           unsigned int      query_size,
                                           unsigned int      subject_size);

double            saft_stats_var          (SaftStatsContext *context,
                                           unsigned int      query_size,
                                           unsigned int      subject_size);

double            saft_stats_pgamma_m_v   (double            d2,
                                           double            mean,
                                           double            var);

double            saft_stats_pgamma       (double            d2,
                                           double            shape,
                                           double            scale);

#ifdef __cplusplus
}
#endif

#endif /* __SAFT_STATS_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
