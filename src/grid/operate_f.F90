!! Copyright (C) 2012 D. Strubbe
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#include "global.h"

! -----------------------------------------------------------------------
!> This module contains interfaces for routines in operate.c
! -----------------------------------------------------------------------
module operate_f_m

  implicit none

  public ! only interfaces in this module

  interface
    subroutine doperate_ri_vec(opn, w, opnri, opri, rimap_inv, rimap_inv_max, fi, ldfp, fo)
      implicit none
      integer, intent(in)    :: opn
      FLOAT,   intent(in)    :: w
      integer, intent(in)    :: opnri
      integer, intent(in)    :: opri
      integer, intent(in)    :: rimap_inv
      integer, intent(in)    :: rimap_inv_max
      FLOAT,   intent(in)    :: fi
      integer, intent(in)    :: ldfp
      FLOAT,   intent(inout) :: fo
    end subroutine doperate_ri_vec
  end interface

  interface
    subroutine zoperate_ri_vec(opn, w, opnri, opri, rimap_inv, rimap_inv_max, fi, ldfp, fo)
      implicit none
      integer, intent(in)    :: opn
      FLOAT,   intent(in)    :: w
      integer, intent(in)    :: opnri
      integer, intent(in)    :: opri
      integer, intent(in)    :: rimap_inv
      integer, intent(in)    :: rimap_inv_max
      CMPLX,   intent(in)    :: fi
      integer, intent(in)    :: ldfp
      CMPLX,   intent(inout) :: fo
    end subroutine zoperate_ri_vec
  end interface

  interface
    subroutine dgauss_seidel(opn, w, opnri, opri, rimap_inv, rimap_inv_max, factor, pot, rho)
      implicit none
      integer, intent(in)    :: opn
      FLOAT,   intent(in)    :: w
      integer, intent(in)    :: opnri
      integer, intent(in)    :: opri
      integer, intent(in)    :: rimap_inv
      integer, intent(in)    :: rimap_inv_max
      FLOAT,   intent(in)    :: factor
      FLOAT,   intent(inout) :: pot
      FLOAT,   intent(in)    :: rho
    end subroutine dgauss_seidel
  end interface

end module operate_f_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
