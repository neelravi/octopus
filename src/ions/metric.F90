!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
  
module metric_m
  use datasets_m
  use geometry_m
  use global_m
  use loct_m
  use math_m
  use messages_m
  use parser_m
  use profiling_m
  use unit_m
  use unit_system_m
  use utils_m
  
  implicit none
  
  private
  
  public ::                       &
    metric_t,                     & 
    metric_init,                  &
    metric_end,                   &
    metric_copy,                  &
    metric_write_info
                   
  type metric_t
    FLOAT :: tensor(1:3, 1:3)
    FLOAT :: br_vecs(1:3)
    FLOAT :: br_angles(1:3)
  end type metric_t

contains

  subroutine metric_init(this)
    type(metric_t),    intent(out) :: this

    PUSH_SUB(metric_init)

    this%tensor(:,:) = M_ZERO 
    this%tensor(1,1) = M_ONE
    this%tensor(2,2) = M_ONE
    this%tensor(3,3) = M_ONE
    
    this%br_vecs(1:3)   = M_ONE
    this%br_angles(1:3) = M_PI*M_HALF


    POP_SUB(metric_init)
  end subroutine metric_init

  ! ---------------------------------------------------------
  subroutine metric_end(this)
    type(kpoints_grid_t), intent(inout) :: this

    PUSH_SUB(metric_end)


    POP_SUB(metric_end)
  end subroutine metric_end


  ! ---------------------------------------------------------
  subroutine metric_copy(bb, aa)
    type(kpoints_grid_t), intent(in)  :: bb
    type(kpoints_grid_t), intent(out) :: aa

    PUSH_SUB(metric_copy)
    

    POP_SUB(metric_copy)
  end subroutine metric_copy


  


  ! ---------------------------------------------------------
  subroutine metric_write_info(this, iunit)
    type(kpoints_t),    intent(in) :: this
    integer,            intent(in) :: iunit
    
    integer :: ik, idir
    character(len=100) :: str_tmp
    character :: index
    
    PUSH_SUB(metric_write_info)
    
    call messages_print_stress(iunit, 'Brillouin zone sampling')

    if(this%method == KPOINTS_MONKH_PACK) then

      call messages_write('Dimensions of the k-point grid      =')
      do idir = 1, this%full%dim
        call messages_write(this%nik_axis(idir), fmt = '(i3,1x)')
      end do
      call messages_new_line()
      
      call messages_write('Total number of k-points            =')
      call messages_write(this%full%npoints)
      call messages_new_line()

      call messages_write('Number of symmetry-reduced k-points =')
      call messages_write(this%reduced%npoints)

      call messages_info(iunit = iunit)
      
    else

      call messages_write('Total number of k-points            =')
      call messages_write(this%full%npoints)
      call messages_new_line()
      call messages_info(iunit = iunit)

    endif

    call messages_new_line()
    call messages_write('List of k-points:')
    call messages_info(iunit = iunit)

    write(message(1), '(6x,a)') 'ik'
    do idir = 1, this%full%dim
      index = index2axis(idir)
      write(str_tmp, '(9x,2a)') 'k_', index
      message(1) = trim(message(1)) // trim(str_tmp)
    enddo
    write(str_tmp, '(6x,a)') 'Weight'
    message(1) = trim(message(1)) // trim(str_tmp)
    message(2) = '---------------------------------------------------------'
    call messages_info(2, iunit)
    
    do ik = 1, kpoints_number(this)
      write(message(1),'(i8,1x)') ik
      do idir = 1, this%full%dim
        write(str_tmp,'(f12.4)') this%reduced%red_point(idir, ik)
        message(1) = trim(message(1)) // trim(str_tmp)
      enddo
      write(str_tmp,'(f12.4)') kpoints_get_weight(this, ik)
      message(1) = trim(message(1)) // trim(str_tmp)
      call messages_info(1, iunit)
    end do

    call messages_info(iunit = iunit)

    call messages_print_stress(iunit)

    POP_SUB(metric_write_info)
  end subroutine metric_write_info
  


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
