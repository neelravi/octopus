!! Copyright (C) 2003 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.

#include "config_F90.h"

#define __STRING(x)     #x

#if defined(DEBUG)
#  define ASSERT(expr) \
  if(.not.(expr)) call assert_die (__STRING(expr), __FILE__, __LINE__)
#else
#  define ASSERT(expr)
!!#  define IN inout
#endif

#if defined(SINGLE_PRECISION)
#  define PRECISION 4
#  define FLOAT     real(4)
#  define CMPLX     complex(4) 
#  define PREC(x)   s ## x
#  define CNST(x)   x ## _4
#else
#  define PRECISION 8
#  define FLOAT     real(8)
#  define CMPLX     complex(8)
#  define PREC(x)   d ## x
#  define CNST(x)   x ## _8
#endif

! what do you wish for dinner, dear?
#if defined(COMPLEX_WFNS)
#  include "complex.F90"
#else
#  include "real.F90"
#endif
