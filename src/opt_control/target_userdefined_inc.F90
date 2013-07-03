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
!! $Id: target_userdefined_inc.F90 $


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_userdefined(gr, tg)
    type(grid_t),     intent(in)    :: gr
    type(target_t),   intent(inout) :: tg

    integer             :: no_states, ib, ip, idim, inst, inik, id, ist, ik
    type(block_t)       :: blk
    FLOAT               :: xx(MAX_DIM), rr, psi_re, psi_im
    PUSH_SUB(target_init_userdefined)

    message(1) =  'Info: Target is a user-defined state.'
    call messages_info(1)
      
    !%Variable OCTTargetUserdefined
    !%Type block
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% Define a target state. Syntax follows the one of the <tt>UserDefinedStates</tt> block.
    !% Example:
    !%
    !% <tt>%OCTTargetUserdefined
    !% <br>&nbsp;&nbsp; 1 | 1 | 1 |  "exp(-r^2)*exp(-i*0.2*x)"
    !% <br>%</tt>
    !%  
    !%End
    if(parse_block(datasets_check('OCTTargetUserdefined'), blk) == 0) then
        
      no_states = parse_block_n(blk)
      do ib = 1, no_states
        call parse_block_integer(blk, ib - 1, 0, idim)
        call parse_block_integer(blk, ib - 1, 1, inst)
        call parse_block_integer(blk, ib - 1, 2, inik)

        ! read formula strings and convert to C strings
        do id = 1, tg%st%d%dim
          do ist = 1, tg%st%nst
            do ik = 1, tg%st%d%nik   
                
              ! does the block entry match and is this node responsible?
              if(.not. (id  ==  idim .and. ist  ==  inst .and. ik  ==  inik    &
                .and. tg%st%st_start  <=  ist .and. tg%st%st_end >= ist) ) cycle
              
              ! parse formula string
              call parse_block_string(                            &
                blk, ib - 1, 3, tg%st%user_def_states(id, ist, ik))
              ! convert to C string
              call conv_to_C_string(tg%st%user_def_states(id, ist, ik))
              
              do ip = 1, gr%mesh%np
                xx = gr%mesh%x(ip, :)
                rr = sqrt(sum(xx(:)**2))
                
                ! parse user-defined expressions
                call parse_expression(psi_re, psi_im, &
                  gr%sb%dim, xx, rr, M_ZERO, tg%st%user_def_states(id, ist, ik))
                ! fill state
                tg%st%zpsi(ip, id, ist, ik) = psi_re + M_zI * psi_im
              end do
              ! normalize orbital
              call zstates_normalize_orbital(gr%mesh, tg%st%d%dim, &
                tg%st%zpsi(:,:, ist, ik))
            end do
          end do
        enddo
      end do
      call parse_block_end(blk)
      call density_calc(tg%st, gr, tg%st%rho)
    else
      message(1) = '"OCTTargetUserdefined" has to be specified as block.'
      call messages_fatal(1)
    end if

    POP_SUB(target_init_userdefined)
  end subroutine target_init_userdefined


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_userdefined(tg, gr, psi) result(j1)
    type(target_t),   intent(inout) :: tg
    type(grid_t),     intent(inout) :: gr
    type(states_t),   intent(inout) :: psi

    integer :: ik, ist
    PUSH_SUB(target_j1_userdefined)

    j1 = M_ZERO
    do ik = 1, psi%d%nik
      do ist = psi%st_start, psi%st_end
        j1 = j1 + psi%occ(ist, ik) * &
          abs(zmf_dotp(gr%mesh, psi%d%dim, psi%zpsi(:, :, ist, ik), &
              tg%st%zpsi(:, :, ist, ik)))**2
      end do
    end do

    POP_SUB(target_j1_userdefined)
  end function target_j1_userdefined

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End: