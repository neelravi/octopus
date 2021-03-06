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

module unocc_m
  use datasets_m
  use density_m
  use eigensolver_m
  use global_m
  use output_m
  use hamiltonian_m
  use io_m
  use lcao_m
  use mesh_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use restart_m
  use simul_box_m
  use states_m
  use states_io_m
  use states_dim_m
  use states_restart_m
  use system_m
  use v_ks_m
  use xc_m

  implicit none

  private
  public :: &
    unocc_run


contains

  ! ---------------------------------------------------------
  subroutine unocc_run(sys, hm, fromscratch)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    logical,             intent(inout) :: fromscratch

    type(eigensolver_t) :: eigens
    integer :: iunit, ierr, iter, ierr_rho, ik
    logical :: read_gs, converged, forced_finish, showoccstates, is_orbital_dependent, occ_missing
    integer :: max_iter, nst_calculated, showstart
    integer :: n_filled, n_partially_filled, n_half_filled
    integer, allocatable :: lowest_missing(:, :), occ_states(:)
    character(len=50) :: str
    character(len=10) :: dirname
    type(restart_t) :: restart_load_unocc, restart_load_gs, restart_dump
    logical :: write_density

    PUSH_SUB(unocc_run)

    !%Variable MaximumIter
    !%Type integer
    !%Default 50
    !%Section Calculation Modes::Unoccupied States
    !%Description
    !% Maximum number of eigensolver iterations. The code will stop even if convergence
    !% has not been achieved. -1 means unlimited. 0 means just do LCAO or read from
    !% restart, and stop.
    !%End
    call parse_integer(datasets_check('MaximumIter'), 50, max_iter)
    call messages_obsolete_variable('UnoccMaximumIter', 'MaximumIter')
    if(max_iter < 0) max_iter = huge(max_iter)

    !%Variable UnoccShowOccStates
    !%Type logical
    !%Default false
    !%Section Calculation Modes::Unoccupied States
    !%Description
    !% If true, the convergence for the occupied states will be shown too in the output.
    !% This is useful for testing, or if the occupied states fail to converge.
    !% It will be enabled automatically if only occupied states are being calculated.
    !%End
    call parse_logical(datasets_check('UnoccShowOccStates'), .false., showoccstates)

    SAFE_ALLOCATE(occ_states(1:sys%st%d%nik))
    do ik = 1, sys%st%d%nik
      call occupied_states(sys%st, ik, n_filled, n_partially_filled, n_half_filled)
      occ_states(ik) = n_filled + n_partially_filled + n_half_filled
    enddo

    call init_(sys%gr%mesh, sys%st)
    converged = .false.

    SAFE_ALLOCATE(lowest_missing(1:sys%st%d%dim, 1:sys%st%d%nik))

    read_gs = .true.
    if (.not. fromScratch) then
      call restart_init(restart_load_unocc, RESTART_UNOCC, RESTART_TYPE_LOAD, sys%st%dom_st_kpt_mpi_grp, &
                        ierr, mesh=sys%gr%mesh)

      if(ierr == 0) then
        call states_load(restart_load_unocc, sys%st, sys%gr, ierr, lowest_missing = lowest_missing)
        call restart_end(restart_load_unocc)
      endif
      
      ! If successfully read states from unocc, do not read from gs.
      ! If RESTART_GS and RESTART_UNOCC have the same directory (the default), and we tried RESTART_UNOCC
      ! already and failed, it is a waste of time to try to read again.
      if(ierr == 0 .or. restart_are_basedirs_equal(RESTART_GS, RESTART_UNOCC)) &
        read_gs = .false.
    end if

    call restart_init(restart_load_gs, RESTART_GS, RESTART_TYPE_LOAD, sys%st%dom_st_kpt_mpi_grp, &
                      ierr_rho, mesh=sys%gr%mesh)

    if(ierr_rho == 0) then
      if (read_gs) then
        call states_load(restart_load_gs, sys%st, sys%gr, ierr, lowest_missing = lowest_missing)
      end if
      call states_load_rho(restart_load_gs, sys%st, sys%gr, ierr_rho)
      write_density = restart_has_map(restart_load_gs)
      call restart_end(restart_load_gs)
    else
      write_density = .true.
    endif

    is_orbital_dependent = (sys%ks%theory_level == HARTREE .or. sys%ks%theory_level == HARTREE_FOCK .or. &
      (sys%ks%theory_level == KOHN_SHAM_DFT .and. xc_is_orbital_dependent(sys%ks%xc)))

    if(is_orbital_dependent) then
      message(1) = "Be sure your gs run is well converged since you have an orbital-dependent functional."
      message(2) = "Otherwise, the occupied states may change in CalculationMode = unocc, and your"
      message(3) = "unoccupied states will not be consistent with the gs run."
      call messages_warning(3)
    endif

    if(ierr_rho /= 0 .or. is_orbital_dependent) then
      occ_missing = .false.
      do ik = 1, sys%st%d%nik
        if(any(lowest_missing(1:sys%st%d%dim, ik) <= occ_states(ik))) then
          occ_missing = .true.
        end if
      enddo

      if(occ_missing) then
        if(is_orbital_dependent) then
          message(1) = "For an orbital-dependent functional, all occupied orbitals must be provided."
        else if(ierr_rho /= 0) then
          message(1) = "Since density could not be read, all occupied orbitals must be provided."
        endif

        message(2) = "Not all the occupied orbitals could be read."
        message(3) = "Please run a ground-state calculation first!"
        call messages_fatal(3)
      endif

      message(1) = "Unable to read density: Building density from wavefunctions."
      call messages_info(1)

      if(.not. hm%cmplxscl%space) then
        call density_calc(sys%st, sys%gr, sys%st%rho)
      else
        call density_calc(sys%st, sys%gr, sys%st%zrho%Re, sys%st%zrho%Im)
      endif
    end if

    if (states_are_real(sys%st)) then
      message(1) = 'Info: Using real wavefunctions.'
    else
      message(1) = 'Info: Using complex wavefunctions.'
    end if
    call messages_info(1)

    if(fromScratch .or. ierr /= 0) then
      if(fromScratch) then
        ! do not use previously calculated occupied states
        nst_calculated = min(maxval(occ_states), minval(lowest_missing) - 1)
      else
        ! or, use as many states as have been calculated
        nst_calculated = minval(lowest_missing) - 1
      end if
      showstart = max(nst_calculated + 1, 1)
      call lcao_run(sys, hm, st_start = showstart)
    else
      ! we successfully read all the states and are planning to use them, no need for LCAO
      call v_ks_calc(sys%ks, hm, sys%st, sys%geo, calc_eigenval = .false.)
      showstart = minval(occ_states(:)) + 1
    end if

    ! it is strange and useless to see no eigenvalues written if you are only calculating
    ! occupied states, on a different k-point.
    if(showstart > sys%st%nst) showstart = 1
    
    SAFE_DEALLOCATE_A(lowest_missing)

    if(showoccstates) showstart = 1

    ! Restart dump should be initialized after restart_load, as the mesh might have changed
    call restart_init(restart_dump, RESTART_UNOCC, RESTART_TYPE_DUMP, sys%st%dom_st_kpt_mpi_grp, &
                      ierr, mesh=sys%gr%mesh)

    ! make sure the density is defined on the same mesh as the wavefunctions that will be written
    if(write_density) &
      call states_dump_rho(restart_dump, sys%st, sys%gr, ierr_rho)

    message(1) = "Info: Starting calculation of unoccupied states."
    call messages_info(1)

    ! reset this variable, so that the eigensolver passes through all states
    eigens%converged(:) = 0

    do iter = 1, max_iter
      call eigensolver_run(eigens, sys%gr, sys%st, hm, 1, converged)

      write(str, '(a,i5)') 'Unoccupied states iteration #', iter
      call messages_print_stress(stdout, trim(str))

      ! If not all gs wavefunctions were read when starting, in particular for nscf with different k-points,
      ! the occupations must be recalculated each time, though they do not affect the result of course.
      ! FIXME: This is wrong for metals where we must use the Fermi level from the original calculation!
      call states_fermi(sys%st, sys%gr%mesh)
      call states_write_eigenvalues(stdout, sys%st%nst, sys%st, sys%gr%sb, eigens%diff, st_start = showstart)
      call messages_print_stress(stdout)

      ! write output file
      if(mpi_grp_is_root(mpi_world)) then
        call io_mkdir(STATIC_DIR)
        iunit = io_open(STATIC_DIR//'/eigenvalues', action='write')
        
        if(converged) then
          write(iunit,'(a)') 'All states converged.'
        else
          write(iunit,'(a)') 'Some of the states are not fully converged!'
        end if
        write(iunit,'(a, e17.6)') 'Criterion = ', eigens%tolerance
        write(iunit,'(1x)')
        call states_write_eigenvalues(iunit, sys%st%nst, sys%st, sys%gr%sb, eigens%diff)
        call io_close(iunit)
      end if

      forced_finish = clean_stop(sys%mc%master_comm)
      
      ! write restart information.
      if(converged .or. (modulo(iter, sys%outp%restart_write_interval) == 0) .or. iter == max_iter .or. forced_finish) then
        call states_dump(restart_dump, sys%st, sys%gr, ierr, iter=iter)
        if(ierr /= 0) then
          message(1) = "Unable to write states wavefunctions."
          call messages_warning(1)
        end if
      end if

      if(sys%outp%output_interval /= 0 .and. mod(iter, sys%outp%output_interval) == 0) then
        write(dirname,'(a,i4.4)') "unocc.",iter
        call output_all(sys%outp, sys%gr, sys%geo, sys%st, hm, sys%ks, dirname)
      end if

      if(converged .or. forced_finish) exit
    end do

    call restart_end(restart_dump)

    if(any(eigens%converged(:) < occ_states(:))) then
      write(message(1),'(a)') 'Some of the occupied states are not fully converged!'
      call messages_warning(1)
    endif

    SAFE_DEALLOCATE_A(occ_states)

    if(.not. converged) then
      write(message(1),'(a)') 'Some of the unoccupied states are not fully converged!'
      call messages_warning(1)
    endif

    if(simul_box_is_periodic(sys%gr%sb).and. sys%st%d%nik > sys%st%d%nspin) then
      call states_write_bands(STATIC_DIR, sys%st%nst, sys%st, sys%gr%sb)
      call states_write_fermi_energy(STATIC_DIR, sys%st, sys%gr%mesh, sys%gr%sb)
    end if

    call output_all(sys%outp, sys%gr, sys%geo, sys%st, hm, sys%ks, STATIC_DIR)

    call end_()
    POP_SUB(unocc_run)

  contains

    ! ---------------------------------------------------------
    subroutine init_(mesh, st)
      type(mesh_t),   intent(in)    :: mesh
      type(states_t), intent(inout) :: st

      PUSH_SUB(unocc_run.init_)

      call messages_obsolete_variable("NumberUnoccStates", "ExtraStates")

      call states_allocate_wfns(st, mesh)

      ! now the eigensolver stuff
      call eigensolver_init(eigens, sys%gr, st)

      POP_SUB(unocc_run.init_)
    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      PUSH_SUB(unocc_run.end_)

      call eigensolver_end(eigens)

      POP_SUB(unocc_run.end_)
    end subroutine end_

  end subroutine unocc_run


end module unocc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
