/**
 *
 */

#include <stdio.h>

#include "bfasterror.h"


static void bfast_error_handler_default (const char *fmt,
                                         va_list     ap);


BfastErrorHandler bfast_error_handler = bfast_error_handler_default;


BfastErrorHandler
bfast_set_error_handler (BfastErrorHandler handler)
{
  BfastErrorHandler ret = bfast_error_handler;
  bfast_error_handler = handler;
  return ret;
}

void
bfast_error (const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);

  if (bfast_error_handler)
    (bfast_error_handler)(fmt, ap);

  va_end(ap);
}


static void
bfast_error_handler_default (const char *fmt,
                             va_list     ap)
{
  fprintf (stderr, "[ERROR] ");
  vfprintf (stderr, fmt, ap);
  fprintf (stderr, ".\n");
}
