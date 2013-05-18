/*!
\file gk_externs.h
\brief This file contains definitions of external variables created by GKlib

\date   Started 3/27/2007
\author George
\version\verbatim $Id: gk_externs.h 10711 2011-08-31 22:23:04Z karypis $ \endverbatim
*/

#include <config.h>

#ifndef _GK_EXTERNS_H_
#define _GK_EXTERNS_H_


/*************************************************************************
* Extern variable definition. Hopefully, the TLS_SPECIFIER makes them thread-safe.
**************************************************************************/
#ifndef _GK_ERROR_C_
/* declared in error.c */
extern TLS_SPECIFIER int gk_cur_jbufs;
extern TLS_SPECIFIER jmp_buf gk_jbufs[];
extern TLS_SPECIFIER jmp_buf gk_jbuf;

#endif

#endif
