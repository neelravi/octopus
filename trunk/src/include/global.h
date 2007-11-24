!! Copyright (C) 2003-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!!
!! $Id$

#include "config_F90.h"

#define NP      gr%m%np
#define NP_PART gr%m%np_part

#define MAX_DIM 3
#define MAX_SPIN 4
#define NDIM    gr%sb%dim
#define LAP     f_der%der_discr%lapl

#if defined(F90_ACCEPTS_LINE_NUMBERS)
# define CARDINAL \newline\cardinal __LINE__ __FILE__
#else
# define CARDINAL \newline
#endif

#define __STRING(x)     #x

#if !defined(NDEBUG) 
#  define ASSERT(expr) \
     if(.not.(expr)) then                    \newline \
       call assert_die (__STRING(expr), &    \newline \
                        __FILE__, __LINE__)  \newline \
     end if                                  \
     CARDINAL

#else
#  define ASSERT(expr)
#endif


# define ALLOCATE(x, size) \
    allocate(x, stat=global_alloc_err)                 \newline \
     if(global_alloc_err.ne.0) then                    \newline \
       call alloc_error((size), __FILE__, __LINE__)    \newline \
     end if                                            \
     CARDINAL

#define REAL_DOUBLE real(8)
#define REAL_SINGLE real(4)

#if defined(SINGLE_PRECISION)
#  define REAL_PRECISION 4
#  define FLOAT     	 real(4)
#  define MPI_FLOAT 	 MPI_REAL
#  define CMPLX     	 complex(4)
#  define MPI_CMPLX 	 MPI_COMPLEX
#  define PREC(x)   	 s ## x
#  define CNST(x)   	 x ## _4
#else
#  define REAL_PRECISION 8
#  define FLOAT          real(8)
#  define MPI_FLOAT 	 MPI_DOUBLE_PRECISION
#  define CMPLX     	 complex(8)
#  define MPI_CMPLX 	 MPI_DOUBLE_COMPLEX
#  define PREC(x)   	 d ## x
#  define CNST(x)   	 x ## _8
#endif

#define   TOFLOAT(x) real(x, REAL_PRECISION)
