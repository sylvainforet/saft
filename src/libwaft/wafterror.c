/* wafterror.c
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

#include "wafterror.h"


static void waft_error_handler_default (const char *fmt,
                                         va_list     ap);


WaftErrorHandler waft_error_handler = waft_error_handler_default;


WaftErrorHandler
waft_set_error_handler (WaftErrorHandler handler)
{
  WaftErrorHandler ret = waft_error_handler;
  waft_error_handler = handler;
  return ret;
}

void
waft_error (const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);

  if (waft_error_handler)
    (waft_error_handler)(fmt, ap);

  va_end(ap);
}


static void
waft_error_handler_default (const char *fmt,
                             va_list     ap)
{
  fprintf (stderr, "[ERROR] ");
  vfprintf (stderr, fmt, ap);
  fprintf (stderr, ".\n");
}
