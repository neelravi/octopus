!! Copyright (C) 2009 X. Andrade
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

module memory_m
  use global_m
  use hamiltonian_m
  use mesh_m
  use messages_m
  use states_m
  use system_m
  use unit_system_m
  
  implicit none

  private
  public :: memory_run

contains

  ! ---------------------------------------------------------
  subroutine memory_run(sys)
    type(system_t),      intent(inout) :: sys

    real(8) :: mesh_global, mesh_local, wfns

    PUSH_SUB(memory_run)

    mesh_global = mesh_global_memory(sys%gr%mesh)
    mesh_local  = mesh_local_memory(sys%gr%mesh)

    call messages_write('Mesh')
    call messages_new_line()

    call messages_write('  global  :')
    call messages_write(mesh_global, units = unit_megabytes, fmt = '(f10.1)')
    call messages_new_line()

    call messages_write('  local   :')
    call messages_write(mesh_local, units = unit_megabytes, fmt = '(f10.1)')
    call messages_new_line()

    call messages_write('  total   :')
    call messages_write(mesh_global + mesh_local, units = unit_megabytes, fmt = '(f10.1)')
    call messages_new_line()

    call messages_info()

    wfns = states_wfns_memory(sys%st, sys%gr%mesh)

    call messages_write('States')
    call messages_new_line()

    call messages_write('  real    :')
    call messages_write(wfns, units = unit_megabytes, fmt = '(f10.1)')
    call messages_write(' (par_kpoints + par_states + par_domains)')
    call messages_new_line()    

    call messages_write('  complex :')
    call messages_write(2.0_8*wfns, units = unit_megabytes, fmt = '(f10.1)')
    call messages_write(' (par_kpoints + par_states + par_domains)')
    call messages_new_line()

    call messages_info()

    POP_SUB(memory_run)

  end subroutine memory_run

end module memory_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
