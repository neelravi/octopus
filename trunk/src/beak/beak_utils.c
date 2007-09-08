/*
 Copyright (C) 2006 X. Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: operate.c 2146 2006-05-23 17:36:00Z xavier $
*/

#include <stdio.h>

#include <config.h>

#include "beak.h"

#define _XOPEN_SOURCE 600
#include <stdlib.h>

#if defined(USE_FAKE_MALLOC)
#include <errno.h>

/* 
Override calls to malloc with calls to posix_memalign. 

Not the most elegant thing in the world, but it appears that most x86
Fortran compilers can't be instructed to align allocated memory to a
16 bytes boundary as required by SSE.
*/

void * malloc (size_t size){
  int err;
  void * ptr;
  err = posix_memalign((void *) &ptr, 16, size);
  if( err == 0 ) return ptr;
  else {
    errno = ENOMEM;
    return NULL;
  }
}

#endif /* USE_FAKE_MALLOC */

int FC_FUNC_(op_is_available, OP_IS_AVAILABLE)
  (int * opid, int * type){
  int result = 1;
  
  if( *opid == OP_OLU_C && *type == M_CMPLX) result = 0;

#if !defined(USE_VECTORS)
  if( *opid == OP_SSE  ) result = 0;
#elif defined(SINGLE_PRECISION)
  if( *opid == OP_SSE   && *type == M_CMPLX ) result = 0;
#endif

  return result;
}

