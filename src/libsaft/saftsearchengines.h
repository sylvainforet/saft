/* saftsearchengine.h
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

#ifndef __SAFT_SEARCH_ENGINE_H__
#define __SAFT_SEARCH_ENGINE_H__

#include "saftsearch.h"

#ifdef __cplusplus
extern "C"
{
#endif

SaftSearchEngine* saft_search_engine_generic       (SaftOptions *options);

SaftSearchEngine* saft_search_engine_dna_array_new (SaftOptions *options);

SaftSearchEngine* saft_search_engine_dna_hash_new  (SaftOptions *options);

#ifdef __cplusplus
}
#endif

#endif /* __SAFT_SEARCH_ENGINE_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */

