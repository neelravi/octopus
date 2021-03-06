!! Copyright (C) 2008 X. Andrade
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

module gauge_field_m
  use datasets_m
  use derivatives_m
  use geometry_m
  use global_m
  use grid_m
  use io_m
  use lalg_basic_m
  use logrid_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use projector_m
  use ps_m
  use restart_m
  use simul_box_m
  use species_m
  use splines_m
  use states_m
  use states_dim_m
  use submesh_m
  use symmetries_m
  use symmetrizer_m
  use symm_op_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private

  public ::                               &
    gauge_force_t,                        &   
    gauge_field_t,                        &
    gauge_field_nullify,                  &
    gauge_field_init,                     &
    gauge_field_init_vec_pot,             &
    gauge_field_is_applied,               &
    gauge_field_set_vec_pot,              &
    gauge_field_set_vec_pot_vel,          &
    gauge_field_get_vec_pot,              &
    gauge_field_get_vec_pot_vel,          &
    gauge_field_get_vec_pot_acc,          &
    gauge_field_propagate,                &
    gauge_field_propagate_vel,            &
    gauge_field_get_energy,               &
    gauge_field_dump,                     &
    gauge_field_load,                     &
    gauge_field_end,                      &
    gauge_field_get_force

  type gauge_force_t
    FLOAT   :: vecpot(1:MAX_DIM)   
  end type gauge_force_t

  type gauge_field_t
    FLOAT   :: vecpot(1:MAX_DIM)   
    FLOAT   :: vecpot_vel(1:MAX_DIM)
    FLOAT   :: vecpot_acc(1:MAX_DIM)    
    FLOAT   :: wp2
    logical :: with_gauge_field
  end type gauge_field_t

contains

  subroutine gauge_field_nullify(this)
    type(gauge_field_t), intent(out) :: this

    PUSH_SUB(gauge_field_nullify)
    this%with_gauge_field = .false.

    POP_SUB(gauge_field_nullify)
  end subroutine gauge_field_nullify

  ! ---------------------------------------------------------
  subroutine gauge_field_init(this, sb)
    type(gauge_field_t),     intent(out)   :: this
    type(simul_box_t),       intent(in)    :: sb

    integer :: ii, iop, iop2, ik
    type(block_t) :: blk

    PUSH_SUB(gauge_field_init)

    this%with_gauge_field = .false.
    this%vecpot = M_ZERO
    this%vecpot_vel = M_ZERO
    this%vecpot_acc = M_ZERO

    !%Variable GaugeVectorField
    !%Type block
    !%Section Hamiltonian
    !%Description
    !% The gauge vector field is used to include a uniform (but time-dependent)
    !% external electric field in a time-dependent run for
    !% a periodic system. An optional second row specifies the initial
    !% value for the time derivative of the gauge field (which is set
    !% to zero by default). By default this field is not included.
    !% If <tt>KPointsUseSymmetries = yes</tt>, then <tt>SymmetryBreakDir</tt>
    !% must be set in the same direction.
    !% This is used with utility <tt>oct-dielectric_function</tt>
    !% according to GF Bertsch, J-I Iwata, A Rubio, and K Yabana,
    !% <i>Phys. Rev. B</i> <b>62</b>, 7998-8002 (2000).
    !%End
    
    ! Read the initial gauge vector field
    
    if(parse_block(datasets_check('GaugeVectorField'), blk) == 0) then
      
      this%with_gauge_field = .true.
      
      do ii = 1, sb%dim
        call parse_block_float(blk, 0, ii - 1, this%vecpot(ii))
      end do
      
      call parse_block_end(blk)
      
      if(.not. simul_box_is_periodic(sb)) then
        message(1) = "GaugeVectorField is intended for periodic systems."
        call messages_warning(1)
      endif

      if(sb%kpoints%use_symmetries) then
        do ik = 1, sb%kpoints%reduced%npoints
          do iop = 1, sb%kpoints%num_symmetry_ops(ik)
            iop2 = sb%kpoints%symmetry_ops(ik, iop)
            if(.not. symm_op_invariant(sb%symm%ops(iop2), this%vecpot, CNST(1e-5))) then
              message(1) = "The GaugeVectorField breaks (at least) one of the symmetries used to reduce the k-points."
              message(2) = "Set SymmetryBreakDir equal to GaugeVectorField."
              call messages_fatal(2)
            endif
          enddo
        enddo
      endif
    end if

    POP_SUB(gauge_field_init)
  end subroutine gauge_field_init


  ! ---------------------------------------------------------
  subroutine gauge_field_end(this)
    type(gauge_field_t),     intent(inout) :: this

    PUSH_SUB(gauge_field_end)
    this%with_gauge_field = .false.

    POP_SUB(gauge_field_end)
  end subroutine gauge_field_end


  ! ---------------------------------------------------------
  logical pure function gauge_field_is_applied(this) result(is_applied)
    type(gauge_field_t),  intent(in) :: this

    is_applied = this%with_gauge_field
  end function gauge_field_is_applied


  ! ---------------------------------------------------------
  subroutine gauge_field_set_vec_pot(this, vec_pot)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: vec_pot(1:MAX_DIM)

    PUSH_SUB(gauge_field_set_vec_pot)
    this%vecpot = vec_pot

    POP_SUB(gauge_field_set_vec_pot)
  end subroutine gauge_field_set_vec_pot


  ! ---------------------------------------------------------
  subroutine gauge_field_set_vec_pot_vel(this, vec_pot_vel)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: vec_pot_vel(1:MAX_DIM)

    PUSH_SUB(gauge_field_set_vec_pot_vel)
    this%vecpot_vel = vec_pot_vel

    POP_SUB(gauge_field_set_vec_pot_vel)
  end subroutine gauge_field_set_vec_pot_vel


  ! ---------------------------------------------------------
  function gauge_field_get_vec_pot(this) result(vec_pot)
    type(gauge_field_t),  intent(in) :: this
    FLOAT :: vec_pot(1:MAX_DIM)

    PUSH_SUB(gauge_field_get_vec_pot)
    vec_pot = this%vecpot

    POP_SUB(gauge_field_get_vec_pot)
  end function gauge_field_get_vec_pot


  ! ---------------------------------------------------------
  function gauge_field_get_vec_pot_vel(this) result(vec_pot_vel)
    type(gauge_field_t),  intent(in) :: this
    FLOAT :: vec_pot_vel(1:MAX_DIM)

    PUSH_SUB(gauge_field_get_vec_pot_vel)
    vec_pot_vel = this%vecpot_vel

    POP_SUB(gauge_field_get_vec_pot_vel)
  end function gauge_field_get_vec_pot_vel


  ! ---------------------------------------------------------
  function gauge_field_get_vec_pot_acc(this) result(vec_pot_acc)
    type(gauge_field_t),  intent(in) :: this
    FLOAT :: vec_pot_acc(1:MAX_DIM)

    PUSH_SUB(gauge_field_get_vec_pot_acc)
    vec_pot_acc = this%vecpot_acc

    POP_SUB(gauge_field_get_vec_pot_acc)
  end function gauge_field_get_vec_pot_acc


  ! ---------------------------------------------------------
  subroutine gauge_field_propagate(this, force, dt)
    type(gauge_field_t),  intent(inout) :: this
    type(gauge_force_t),  intent(in)    :: force
    FLOAT,                intent(in)    :: dt

    PUSH_SUB(gauge_field_propagate)

    this%vecpot_acc(1:MAX_DIM) = force%vecpot(1:MAX_DIM)

    this%vecpot = this%vecpot + dt * this%vecpot_vel + M_HALF * dt**2 * force%vecpot

    POP_SUB(gauge_field_propagate)
  end subroutine gauge_field_propagate


  ! ---------------------------------------------------------
  subroutine gauge_field_propagate_vel(this, force, dt)
    type(gauge_field_t),  intent(inout) :: this
    type(gauge_force_t),  intent(in)    :: force
    FLOAT,                intent(in)    :: dt

    PUSH_SUB(gauge_field_propagate_vel)
    this%vecpot_vel = this%vecpot_vel + M_HALF * dt * (this%vecpot_acc + force%vecpot)

    POP_SUB(gauge_field_propagate_vel)
  end subroutine gauge_field_propagate_vel


  ! ---------------------------------------------------------
  subroutine gauge_field_init_vec_pot(this, sb, st)
    type(gauge_field_t),  intent(inout) :: this
    type(simul_box_t),    intent(in)    :: sb
    type(states_t),       intent(in)    :: st
    
    PUSH_SUB(gauge_field_init_vec_pot)

    this%wp2 = M_FOUR*M_PI*st%qtot/sb%rcell_volume

    write (message(1), '(a,f12.6,a)') "Info: Electron-gas plasmon frequency", &
         units_from_atomic(units_out%energy, sqrt(this%wp2)), " ["//trim(units_abbrev(units_out%energy))//"]"
    call messages_info(1)

    POP_SUB(gauge_field_init_vec_pot)
  end subroutine gauge_field_init_vec_pot

  ! ---------------------------------------------------------
  FLOAT function gauge_field_get_energy(this, sb) result(energy)
    type(gauge_field_t),  intent(in)    :: this
    type(simul_box_t),    intent(in)    :: sb

    PUSH_SUB(gauge_field_get_energy)
    energy = sb%rcell_volume / (M_EIGHT * M_PI * P_c**2) * sum(this%vecpot_vel(1:MAX_DIM)**2)

    POP_SUB(gauge_field_get_energy)
  end function gauge_field_get_energy


  ! ---------------------------------------------------------
  subroutine gauge_field_dump(restart, gfield, ierr)
    type(restart_t),      intent(in)  :: restart
    type(gauge_field_t),  intent(in)  :: gfield
    integer,              intent(out) :: ierr

    integer :: err
    FLOAT :: vecpot(MAX_DIM, 2)
    
    PUSH_SUB(gauge_field_dump)

    ierr = 0
    
    if (restart_skip(restart)) then
      POP_SUB(gauge_field_dump)
      return
    end if

    if (in_debug_mode) then
      message(1) = "Debug: Writing gauge field restart."
      call messages_info(1)
    end if

    vecpot = M_ZERO
    vecpot(:,1) = gauge_field_get_vec_pot(gfield)
    vecpot(:,2) = gauge_field_get_vec_pot_vel(gfield)

    call drestart_write_binary(restart, "gauge_field", 2*MAX_DIM, vecpot, err)
    if (err /= 0) ierr = ierr + 1

    if (in_debug_mode) then
      message(1) = "Debug: Writing gauge field restart done."
      call messages_info(1)
    end if

    POP_SUB(gauge_field_dump)
  end subroutine gauge_field_dump


  ! ---------------------------------------------------------
  subroutine gauge_field_load(restart, gfield, ierr)
    type(restart_t),      intent(in)    :: restart
    type(gauge_field_t),  intent(inout) :: gfield
    integer,              intent(out)   :: ierr

    integer :: err
    FLOAT :: vecpot(MAX_DIM, 2)
    
    PUSH_SUB(gauge_field_load)

    ierr = 0
    
    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(gauge_field_load)
      return
    end if

    if (in_debug_mode) then
      message(1) = "Debug: Reading gauge field restart."
      call messages_info(1)
    end if

    call drestart_read_binary(restart, "gauge_field", 2*MAX_DIM, vecpot, err)
    if (err /= 0) ierr = ierr + 1

    call gauge_field_set_vec_pot(gfield, vecpot(:,1))
    call gauge_field_set_vec_pot_vel(gfield, vecpot(:,2))

    if (in_debug_mode) then
      message(1) = "Debug: Reading gauge field restart done."
      call messages_info(1)
    end if

    POP_SUB(gauge_field_load)
  end subroutine gauge_field_load

  ! ---------------------------------------------------------

  subroutine gauge_field_get_force(gr, st, force)
    type(grid_t),         intent(in)    :: gr
    type(states_t),       intent(in)    :: st
    type(gauge_force_t),  intent(out)   :: force

    integer :: idir

    PUSH_SUB(gauge_field_get_force)

    ASSERT(st%d%nspin == 1)

    do idir = 1, gr%sb%dim
      force%vecpot(idir) = CNST(4.0)*M_PI*P_c/gr%sb%rcell_volume*dmf_integrate(gr%mesh, st%current(:, idir, 1))
    end do

    POP_SUB(gauge_field_get_force)
  end subroutine gauge_field_get_force


end module gauge_field_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
