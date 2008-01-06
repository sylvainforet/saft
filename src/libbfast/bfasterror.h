/**
 *
 */

#ifndef __BFAST_ERROR_H__
#define __BFAST_ERROR_H__

#include "stdarg.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef void      (*BfastErrorHandler)    (const char *fmt, va_list ap);

BfastErrorHandler bfast_set_error_handler (BfastErrorHandler handler);

void              bfast_error             (const char *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif /* __BFAST_ERROR_H__ */
