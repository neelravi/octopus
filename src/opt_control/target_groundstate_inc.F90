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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: target_groundstate_inc.F90 $


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_groundstate(gr, tg)
    type(grid_t),     intent(in)    :: gr
    type(target_t),   intent(inout) :: tg

    integer :: ierr
    PUSH_SUB(target_init_groundstate)

    message(1) =  'Info: Using Ground State for TargetOperator'
    call messages_info(1)
    call restart_read(trim(restart_dir)//GS_DIR, tg%st, gr, ierr, exact = .true.)

    POP_SUB(target_init_groundstate)
  end subroutine target_init_groundstate


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_groundstate(tg, gr, psi) result(j1)
    type(target_t),   intent(inout) :: tg
    type(grid_t),     intent(inout) :: gr
    type(states_t),   intent(inout) :: psi

    integer :: ist, ik
    PUSH_SUB(target_j1_groundstate)

    do ik = 1, psi%d%nik
      do ist = psi%st_start, psi%st_end
        j1 = j1 + psi%occ(ist, ik) * &
          abs(zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, ist, ik), &
              tg%st%zpsi(:, :, ist, ik)))**2
      end do
    end do

    POP_SUB(target_j1_groundstate)
  end function target_j1_groundstate



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End: