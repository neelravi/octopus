!! Copyright (C) 2014 M. Marques, A. Castro, A. Rubio, G. Bertsch, J.Jornet-Somoza
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
!! $Id: local_write.F90 11872 2014-03-12 19:43:35Z dstrubbe $

#include "global.h"

module local_write_m
  use box_union_m
  use c_pointer_m
  use datasets_m
  use geometry_m
  use global_m
  use grid_m
  use io_m
  use io_function_m
  use kick_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use mpi_debug_m
  use mpi_lib_m
  use parser_m
  use profiling_m
  use species_m
  use states_m
  use unit_m
  use unit_system_m
  use varinfo_m
  use write_iter_m

  implicit none

  private
  public ::         &
    local_write_t,     &
    local_write_init,  &
    local_write_end,   &
    local_write_iter

  integer, parameter ::   &
    LOCAL_OUT_MULTIPOLES  =  1, &
    LOCAL_OUT_DENSITY     =  2, &
    LOCAL_OUT_POTENTIAL   =  3, &
    LOCAL_OUT_MAX         =  3
  
  type local_write_prop_t
    private
    type(c_ptr) :: handle
    logical :: write = .false.
  end type local_write_prop_t

  type local_write_t
    private
    type(local_write_prop_t),allocatable :: out(:,:)
    integer                  :: how              !< how to output
    integer                  :: lmax             !< maximum multipole moment to output
  end type local_write_t

contains

  ! ---------------------------------------------------------
  subroutine local_write_init(writ, nd, lab, iter, dt)
    type(local_write_t), intent(out)   :: writ
    integer,             intent(in)    :: nd 
    character(len=15),   intent(in)    :: lab(:)
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt

    integer :: first, id, flags, iout, default

    PUSH_SUB(local_write_init)

    ! FIXME: if and when these routines are called from a normal run, the Section can be Output.
    ! but then it will need to be in a different folder, since src/util is not linked by the other folders.

    !%Variable LDOutput
    !%Type flag
    !%Default multipoles 
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Defines what should be output during the local domains 
    !% simulation. Many of the options can increase the computational
    !% cost of the simulation, so only use the ones that you need. In
    !% most cases the default value is enough, as it is adapted to the
    !% details. 
    !%Option multipoles 1
    !% Outputs the (electric) multipole moments of the density to the file <tt>td.general/multipoles</tt>.
    !% This is required to, <i>e.g.</i>, calculate optical absorption spectra of finite systems. The
    !% maximum value of <math>l</math> can be set with the variable <tt>LDMultipoleLmax</tt>.
    !%Option density 2
    !% If set (and if the atoms are allowed to move), outputs the coordinates, velocities,
    !% and forces of the atoms to the the file <tt>td.general/coordinates</tt>. On by default if <tt>MoveIons = yes</tt>.
    !%Option local_v 128
    !% If set, <tt>octopus</tt> outputs the different components of the Coulomb potential
    !% to the file <tt>td.general/energy</tt>. Will be zero except for every <tt>TDEnergyUpdateIter</tt> iterations.
    !%End

    default = 2**(LOCAL_OUT_MULTIPOLES - 1) 

    call parse_integer(datasets_check('LDOutput'), default, flags)

    if(.not.varinfo_valid_option('LDOutput', flags, is_flag = .true.)) call input_error('LDOutput')

    SAFE_ALLOCATE(writ%out(LOCAL_OUT_MAX, nd))
    do iout = 1, LOCAL_OUT_MAX
      writ%out(iout,:)%write = (iand(flags, 2**(iout - 1)) /= 0)
    end do

    !%Variable LDOutputHow
    !%Type flag
    !%Default 0
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Describes the format of the output files (see <tt>Output</tt>).
    !% It can take the same values as OutputHow flag.
    !%End
    call parse_integer(datasets_check('LDOutputHow'), 0, writ%how)
    if(.not.varinfo_valid_option('OutputHow', writ%how, is_flag=.true.)) then
      call input_error('LDOutputHow')
    end if

    !%Variable LDMultipoleLmax
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Maximum electric multipole of the density output to the file <tt>local.multipoles/<>domain%<>.multipoles</tt>
    !% during a time-dependent simulation. Must be non-negative.
    !%End

    call parse_integer(datasets_check('LDMultipoleLmax'), 1, writ%lmax)
    if (writ%lmax < 0) then
      write(message(1), '(a,i6,a)') "Input: '", writ%lmax, "' is not a valid LDMultipoleLmax."
      message(2) = '(Must be LDMultipoleLmax >= 0 )'
      call messages_fatal(2)
    end if

    call io_mkdir('local.general')

    if(mpi_grp_is_root(mpi_world)) then
      do id = 1, nd
        if(writ%out(LOCAL_OUT_MULTIPOLES, id)%write) then 
          call io_mkdir('local.general/multipoles')
          call write_iter_init(writ%out(LOCAL_OUT_MULTIPOLES,id)%handle, &
            iter, units_from_atomic(units_out%time, dt), &
          trim(io_workpath("local.general/multipoles/"//trim(lab(id))//".multipoles")))
        end if

        if(writ%out(LOCAL_OUT_POTENTIAL, id)%write) then
          message(1) = 'Potential output is still not implemented'
          call messages_warning(1)
          !call io_mkdir('local.general/potential')
          !call write_iter_init(writ%out(LOCAL_OUT_POTENTIAL)%handle(id), first, &
          !units_from_atomic(units_out%time, dt), trim(io_workpath("local.general/potential/"/trim(lab(id))/".potential")))
        end if

        if(writ%out(LOCAL_OUT_DENSITY,id)%write) then
          message(1) = 'Density output is still not implemented'
          call messages_warning(1)
          !call io_mkdir('local.general/densities')
          !call write_iter_init(writ%out(LOCAL_OUT_DENSITY)%handle(id), first, &
          !units_from_atomic(units_out%time, dt), trim(io_workpath("local.general/densities/"/trim(lab(id))/".densities")))
        end if
      end do
    end if

    POP_SUB(local_write_init)
  end subroutine local_write_init

  ! ---------------------------------------------------------
  subroutine local_write_end(writ)
    type(local_write_t), intent(inout) :: writ
    
    PUSH_SUB(local_write_end)

    SAFE_DEALLOCATE_A(writ%out)

    POP_SUB(local_write_end)
  end subroutine local_write_end

  ! ---------------------------------------------------------
  subroutine local_write_iter(writ, nd, domain, lab, inside, center, gr, st, & 
                              geo, kick, iter, dt, l_start, ldoverwrite)
    type(local_write_t),    intent(inout) :: writ
    integer,                intent(in)    :: nd 
    type(box_union_t),      intent(in)    :: domain(:)
    character(len=15),      intent(in)    :: lab(:)
    logical,                intent(in)    :: inside(:,:)
    FLOAT  ,                intent(in)    :: center(:,:)
    type(grid_t),           intent(inout) :: gr
    type(states_t),         intent(inout) :: st
    type(geometry_t),       intent(inout) :: geo
    type(kick_t),           intent(inout) :: kick
    integer,                intent(in)    :: iter
    FLOAT,                  intent(in)    :: dt
    integer,                intent(in)    :: l_start
    logical,                intent(in)    :: ldoverwrite

    type(profile_t), save :: prof
    integer :: id

    PUSH_SUB(local_write_iter)
    call profiling_in(prof, "LOCAL_WRITE_ITER")

    if(any(writ%out(LOCAL_OUT_MULTIPOLES,:)%write)) then
      if(.not.ldoverwrite)then
      end if
      call local_write_multipole(writ%out(LOCAL_OUT_MULTIPOLES,:), nd, domain, lab, inside, center, & 
                        gr, geo, st, writ%lmax, kick, iter, l_start, ldoverwrite, writ%how)
      do id = 1, nd
        call write_iter_flush(writ%out(LOCAL_OUT_MULTIPOLES, id)%handle)
      end do
    end if

    !if(writ%out(LOCAL_OUT_DENSITY)%write) &
    !  call local_write_density(writ%out(LOCAL_OUT_DENSITY)%handle, hm, iter, geo%kinetic_energy)
    
    !if(writ%out(LOCAL_OUT_POTENTIAL)%write) &
    !  call local_write_potential(writ%out(LOCAL_OUT_POTENTIAL)%handle, hm, iter, geo%kinetic_energy)

    call profiling_out(prof)
    POP_SUB(local_write_iter)
  end subroutine local_write_iter

  ! ---------------------------------------------------------
  subroutine local_write_multipole(out_multip, nd, domain, lab, inside, center, & 
                                gr, geo, st, lmax, kick, iter, l_start, start, how)
    type(local_write_prop_t),      intent(inout) :: out_multip(:)
    integer,                  intent(in)    :: nd 
    type(box_union_t),        intent(in)    :: domain(:)
    character(len=15),        intent(in)    :: lab(:)
    logical,                  intent(in)    :: inside(:,:)
    FLOAT,                    intent(in)    :: center(:,:)
    type(grid_t),         intent(in) :: gr
    type(geometry_t),     intent(in) :: geo
    type(states_t),       intent(in) :: st
    integer,              intent(in) :: lmax
    type(kick_t),         intent(in) :: kick
    integer,              intent(in) :: iter
    integer,              intent(in) :: l_start
    logical,              intent(in) :: start
    integer,              intent(in) :: how

    integer :: id, is, ll, mm, add_lm
    character(len=120) :: aux
    FLOAT, allocatable :: ionic_dipole(:,:), multipole(:,:,:)
    CMPLX, allocatable :: zmultipole(:,:,:)
    logical :: cmplxscl

    PUSH_SUB(local_write_multipole)

    cmplxscl = .false.
    if(st%cmplxscl%space) cmplxscl = .true.

    if(mpi_grp_is_root(mpi_world).and. iter == l_start .and. start) then
      do id = 1, nd   
        call local_write_print_header_init(out_multip(id)%handle)
  
        write(aux, '(a15,i2)')      '# nspin        ', st%d%nspin
        call write_iter_string(out_multip(id)%handle, aux)
        call write_iter_nl(out_multip(id)%handle)
  
        write(aux, '(a15,i2)')      '# lmax         ', lmax
        call write_iter_string(out_multip(id)%handle, aux)
        call write_iter_nl(out_multip(id)%handle)

        call kick_write(kick, out = out_multip(id)%handle)

        call write_iter_header_start(out_multip(id)%handle)

        do is = 1, st%d%nspin
          write(aux,'(a18,i1,a1)') 'Electronic charge(', is,')'; call write_iter_header(out_multip(id)%handle, aux)
          if(lmax>0) then
            write(aux, '(a3,a1,i1,a1)') '<x>', '(', is,')'; call write_iter_header(out_multip(id)%handle, aux)
            if(cmplxscl) call write_iter_header(out_multip(id)%handle, ' ')   
            write(aux, '(a3,a1,i1,a1)') '<y>', '(', is,')'; call write_iter_header(out_multip(id)%handle, aux)
            if(cmplxscl) call write_iter_header(out_multip(id)%handle, ' ')   
            write(aux, '(a3,a1,i1,a1)') '<z>', '(', is,')'; call write_iter_header(out_multip(id)%handle, aux)
            if(cmplxscl) call write_iter_header(out_multip(id)%handle, ' ')   
          end if
          do ll = 2, lmax
            do mm = -ll, ll
              write(aux, '(a2,i2,a4,i2,a2,i1,a1)') 'l=', ll, ', m=', mm, ' (', is,')'
              call write_iter_header(out_multip(id)%handle, aux)
            end do
          end do
        end do
        call write_iter_nl(out_multip(id)%handle)

        ! units
        call write_iter_string(out_multip(id)%handle, '#[Iter n.]')
        call write_iter_header(out_multip(id)%handle, '[' // trim(units_abbrev(units_out%time)) // ']')

        do is = 1, st%d%nspin
          do ll = 0, lmax
            do mm = -ll, ll
              select case(ll)
              case(0)
                call write_iter_header(out_multip(id)%handle, 'Electrons')
              case(1)
                call write_iter_header(out_multip(id)%handle, '[' // trim(units_abbrev(units_out%length)) // ']')
                if(cmplxscl) call write_iter_header(out_multip(id)%handle, ' ')   
              case default
                write(aux, '(a,a2,i1)') trim(units_abbrev(units_out%length)), "**", ll
                call write_iter_header(out_multip(id)%handle, '[' // trim(aux) // ']')
                if(cmplxscl) call write_iter_header(out_multip(id)%handle, ' ')   
              end select
            end do
          end do
        end do
        call write_iter_nl(out_multip(id)%handle)

        ! complex quantities
        if(cmplxscl) then
          call write_iter_string(out_multip(id)%handle, '#       _         ')
          call write_iter_header(out_multip(id)%handle, ' ')

          do is = 1, st%d%nspin
            do ll = 0, lmax
              do mm = -ll, ll
                select case(ll)
                case(0)
                  call write_iter_header(out_multip(id)%handle, ' ')
                case(1)
                  call write_iter_header(out_multip(id)%handle, 'Re')
                  call write_iter_header(out_multip(id)%handle, 'Im')   
                case default
                  call write_iter_header(out_multip(id)%handle, 'Re')
                  call write_iter_header(out_multip(id)%handle, 'Im')   
                end select
              end do
            end do
          end do
          call write_iter_nl(out_multip(id)%handle)

        end if
      
        call local_write_print_header_end(out_multip(id)%handle)
        call write_iter_flush(out_multip(id)%handle)
      end do
    end if

    SAFE_ALLOCATE(ionic_dipole(1:gr%mesh%sb%dim, nd))
    SAFE_ALLOCATE(multipole(1:(lmax + 1)**2, 1:st%d%nspin, nd))
    ionic_dipole(:,:) = M_ZERO
    multipole   (:,:,:) = M_ZERO
    if(cmplxscl) then
      SAFE_ALLOCATE(zmultipole(1:(lmax + 1)**2, 1:st%d%nspin, nd))
      zmultipole(:,:,:) = M_z0
    end if

    do is = 1, st%d%nspin
      if(.not. cmplxscl) then
        call dmf_local_multipoles(gr%mesh, nd, st%rho(:,is), lmax, multipole(:,is,:), inside)
      else
        message(1) = 'Local Multipoles is still not implemented for complex densities'
        call messages_fatal(1)
        !FIXME: modify X(mf_local_multipoles) to deal with complex rho
        !call zmf_local_multipoles(gr%mesh, st%zrho%Re(:,is) + M_zI * st%zrho%Im(:,is), lmax,&
        !  zmultipole(:,is), inside, cmplxscl_th = st%cmplxscl%theta, inside)
        !multipole (:,is) = real(zmultipole(:,is)) ! it should be real anyways 
      end if 
    end do
    ! FIXME: with cmplxscl we need to think how to treat 
    ! the ions dipole moment 
    call local_geometry_dipole(nd, domain, geo, ionic_dipole)
    do is = 1, st%d%nspin
      do id = 1, nd
        multipole(2:gr%mesh%sb%dim+1, is, id) = -ionic_dipole(1:gr%mesh%sb%dim, id)/st%d%nspin & 
                                                - multipole(2:gr%mesh%sb%dim+1, is, id)
      end do
    end do

    if(mpi_grp_is_root(mpi_world)) then
      do id = 1, nd
        call write_iter_set(out_multip(id)%handle, iter)
        call write_iter_start(out_multip(id)%handle)
        do is = 1, st%d%nspin
          add_lm = 1
          do ll = 0, lmax
            do mm = -ll, ll
              if(cmplxscl .and. ll > 0 ) then
                message(1) = 'Local Multipoles is still not implemented for complex densities'
                call messages_fatal(1)
                !FIXME: to deal with complex rho
                !call write_iter_double(out_multip(id)%handle, units_from_atomic(units_out%length**ll,&
                !  real(zmultipole(add_lm, is), REAL_PRECISION)), 1)
                !call write_iter_double(out_multip(id)%handle, units_from_atomic(units_out%length**ll,&
                !  aimag(zmultipole(add_lm, is))), 1)
              else
                call write_iter_double(out_multip(id)%handle, units_from_atomic(units_out%length**ll, &
                                        multipole(add_lm, is, id)), 1)
              end if
            add_lm = add_lm + 1
            end do
          end do
        end do
        call write_iter_nl(out_multip(id)%handle)
      end do
    end if

   ! Write multipoles in BILD format
    if(iand(how, C_OUTPUT_HOW_BILD) /= 0 )then
      !FIXME: to include spin larger than 1.
      is = 1
      do id = 1, nd
        call out_bld_multipoles(multipole(2:4, is, id), center(:,id), lab(id), iter)
      end do
    end if

    SAFE_DEALLOCATE_A(ionic_dipole)
    SAFE_DEALLOCATE_A(multipole)
    SAFE_DEALLOCATE_A(zmultipole)
    POP_SUB(local_write_multipole)
  end subroutine local_write_multipole

  ! ---------------------------------------------------------
  subroutine local_geometry_dipole(nd, dom, geo, dipole)
    integer,           intent(in)  :: nd 
    type(box_union_t), intent(in)  :: dom(:)
    type(geometry_t),  intent(in)  :: geo
    FLOAT,             intent(inout) :: dipole(:,:)

    integer :: ia, id

    PUSH_SUB(local_geometry_dipole)

    dipole(:,:) = M_ZERO
    do ia = 1, geo%natoms
      do  id = 1, nd
        if (box_union_inside(dom(id), geo%atom(ia)%x)) then
          dipole(1:geo%space%dim, id) = dipole(1:geo%space%dim, id) + &
          species_zval(geo%atom(ia)%spec)*(geo%atom(ia)%x(1:geo%space%dim))
        end if
      end do
    end do
    dipole = P_PROTON_CHARGE*dipole

    POP_SUB(local_geometry_dipole)
  end subroutine local_geometry_dipole

  ! ---------------------------------------------------------
  subroutine out_bld_multipoles(multipoles, center, label, iter)
    FLOAT,         intent(in) :: multipoles(:)
    FLOAT,         intent(in) :: center(:)
    character(15), intent(in) :: label
    integer,       intent(in) :: iter
   
    integer             :: ll, out_bld
    character(len=80)   :: filename, folder
    FLOAT               :: dipolearrow(3,2)

    PUSH_SUB(out_bld_multipoles)
    
    write(folder,'(a,a)')'local.general/multipoles/',trim(label)
    call io_mkdir(folder)
    write(filename,'(a,a,a,a,i7.7,a)')trim(folder),'/',trim(label),'.',iter,'.bld'
    out_bld = io_open(file=trim(filename), action='write')

    write(out_bld,'(a,a,a,i7)')'.comment ** Arrow for the dipole moment centered at the center of mass for ', &
                        trim(label), ' domain and iteration number: ',iter
    write(out_bld,'(a)')''
    write(out_bld,'(a)')'.color red'
    write(out_bld,'(a,3(f12.6,2x),a)')'.sphere ',(units_from_atomic(units_out%length,center(ll)), ll= 1, 3),' 0.2' 
    do ll = 1, 3
      dipolearrow(ll,1) = units_from_atomic(units_out%length, center(ll) - multipoles(ll))
      dipolearrow(ll,2) = units_from_atomic(units_out%length, center(ll) + multipoles(ll))
    end do
    write(out_bld,'(a,6(f12.6,2x),a)')'.arrow ',(dipolearrow(ll,1), ll= 1, 3), &
                                     (dipolearrow(ll,2), ll= 1, 3), ' 0.1 0.5 0.90'
    call io_close(out_bld)

    POP_SUB(out_bld_multipoles)
  end subroutine out_bld_multipoles

  ! ---------------------------------------------------------
  subroutine local_write_print_header_init(out)
    type(c_ptr), intent(inout) :: out

    PUSH_SUB(local_write_print_header_init)
    call write_iter_clear(out)
    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)
    call write_iter_string(out,'# HEADER')
    call write_iter_nl(out)

    POP_SUB(local_write_print_header_init)
  end subroutine local_write_print_header_init


  ! ---------------------------------------------------------
  subroutine local_write_print_header_end(out)
    type(c_ptr), intent(inout) :: out

    PUSH_SUB(local_write_print_header_end)

    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)

    POP_SUB(local_write_print_header_end)
  end subroutine local_write_print_header_end

end module local_write_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End: