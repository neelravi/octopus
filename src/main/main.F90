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

program octopus
  use calc_mode_m
  use command_line_m
  use datasets_m
  use global_m
  use io_m
  use loct_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use run_m
  use string_m
  use utils_m
  use varinfo_m

  implicit none

  character*256 :: config_str
  integer :: ns, inp_calc_mode, ierr
  type(block_t) :: blk

  call getopt_init(ierr)
  config_str = trim(get_config_opts()) // trim(get_optional_libraries())
  if(ierr  ==  0) call getopt_octopus(trim(config_str))
  call getopt_end()

  call global_init()
  call messages_init()

  !%Variable ReportMemory
  !%Type logical
  !%Default no
  !%Section Execution::Debug
  !%Description
  !% If true, <tt>Octopus</tt> will print as part of the screen output
  !% information about the memory the code is using. The quantity
  !% reported is an approximation to the size of the heap and
  !% generally it is a lower bound to the actual memory <tt>Octopus</tt> is
  !% using. By default this variable is set to false.
  !%End
  call parse_logical('ReportMemory', .false., conf%report_memory)

  ! need to find out calc_mode already here since some of the variables here (e.g.
  ! periodic dimensions) can be different for the subsystems

  !%Variable CalculationMode
  !%Type integer
  !%Default gs
  !%Section Calculation Modes
  !%Description
  !% Decides what kind of calculation is to be performed.
  !%Option gs 01
  !% Calculation of the ground state.
  !%Option unocc 02
  !% Calculation of unoccupied/virtual KS states. Can also be used for a non-self-consistent
  !% calculation of states at arbitrary k-points, if <tt>density.obf</tt> from <tt>gs</tt>
  !% is provided in the <tt>restart/gs</tt> directory.
  !%Option td 03
  !% Time-dependent calculation (experimental for periodic systems).
  !%Option go 05
  !% Optimization of the geometry.
  !%Option opt_control 07
  !% Optimal control.
  !%Option em_resp 08
  !% Calculation of the electromagnetic response: electric
  !% polarizabilities and hyperpolarizabilities and magnetic
  !% susceptibilities (experimental for periodic systems).
  !% Sternheimer variables may be prefixed with "EM".
  !%Option casida 09
  !% Excitations via Casida linear-response TDDFT; for finite systems only.
  !%Option td_transport 10
  !% Time-dependent quantum transport (experimental).
  !%Option vdw 11
  !% Calculate van der Waals coefficients. Sternheimer variables may be prefixed with "VdW".
  !%Option vib_modes 12
  !% Calculation of the vibrational modes. Sternheimer variables may be prefixed with "VM".
  !%Option one_shot 14
  !% Obsolete. Use <tt>gs</tt> with <tt>MaximumIter = 0</tt> instead.
  !%Option kdotp 15
  !% Calculation of effective masses by <i>k.p</i> perturbation theory (experimental).
  !% Sternheimer variables may be prefixed with "KdotP".
  !%Option gcm 16
  !% Generator-Coordinates Method calculation (experimental).
  !% Ref. K. Capelle, <i>J. Chem. Phys.</i> <b>119</b>, 1285 (2003).
  !%Option dummy 17
  !% This calculation mode does nothing. Useful for debugging, testing and benchmarking.  
  !%Option invert_ks 18
  !% Invert the Kohn-Sham equations (experimental).
  !%Option recipe 99
  !% Prints out a tasty recipe.
  !%
  !% May also be used as a block for multi-dataset mode. The first line is a list of calculation modes,
  !% the second is labels (optional; default is <tt>dsXX_</tt> where XX=01, 02, ...), and the third is the
  !% order for the runs (optional; default is the order listed). The labels will be used as prefixes for
  !% the output directories such as <tt>restart</tt>, <tt>static</tt>, <tt>td...</tt> etc., and can also
  !% be used as prefixes to variables to indicate that they apply only to the calculation for that dataset.
  !%
  !% Example:
  !% <pre>%CalculationMode
  !%  gs              | unocc
  !%  "ground_state_" | "excited_states_"
  !%  1               | 2
  !% %
  !% excited_states_RestartDir = "ground_state_restart"
  !% excited_states_ExtraStates = 9
  !% </pre>
  !%End
  if(parse_block('CalculationMode', blk) == 0) then
    call datasets_init(inp_calc_mode, blk)
  else
    call parse_integer('CalculationMode', CM_GS, inp_calc_mode)
    if(.not.varinfo_valid_option('CalculationMode', inp_calc_mode)) call input_error('CalculationMode')
    call datasets_init(inp_calc_mode)
  end if

  ! Now we can initialize the I/O
  call io_init()

  ! loop over all datasets
  datasets: do ns = 1, no_datasets

    ! set system label
    current_dataset = dataset_run_order(ns)
    current_label = trim(dataset_label(current_dataset))
    call calc_mode_init()

    ! datasets have to be available before calling the _init() functions below
    call io_init_datasets()

    ! now we declare octopus as running
    call io_switch_status('running')

    call profiling_init()

    call print_header()

    if(no_datasets > 1) then
      message(1) = 'Info: Multi-Dataset Mode'
      message(2) = 'Info: Running dataset "'//trim(current_label)//'"'
      call messages_info(2, stress = .true.)
    end if

    ! now we really start
    call run_init(dataset_runmode(current_dataset))
    call run()
    call run_end()
    
#if defined(HAVE_MPI)
    ! wait for all processors to finish
    call MPI_Barrier(mpi_world%comm, mpi_err)
#endif

    ! run finished successfully
    call io_switch_status('finished')
    call io_end()

    call profiling_end()

    call calc_mode_end()
    
    call print_date("Calculation ended on ")
    call print_walltime()
  end do datasets

  call datasets_end()
  call messages_end()
  call global_end()

contains

  subroutine print_walltime()
    integer :: days, hours, min, sec, usec

    call loct_gettimeofday(sec, usec)
    call epoch_time_diff(sec, usec)
    
    days  = sec / 86400
    hours = (sec / 3600) - (days * 24)
    min   = (sec / 60) - (days * 1440) - (hours * 60)
    sec   = modulo(sec, 60)

    message(2) = ''
    if(days  > 0) write(message(2), '(i3,a)') days, ' days,'
    if(hours > 0.or.message(2) /= '') &
      write(message(2), '(a,1x,i2.2,a)') trim(message(2)), hours, 'h'
    if(min   > 0.or.message(1) /= '') &
      write(message(2), '(a,1x,i2.2,a)') trim(message(2)), min, 'm'
    write(message(2), '(a,1x,i2.2,a,i3,a)') trim(message(2)), sec, '.', usec/1000, 's'
    message(1) = str_center('Walltime: ' // trim(message(2)), 70)
    call messages_info(1)

  end subroutine print_walltime

end program octopus

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
