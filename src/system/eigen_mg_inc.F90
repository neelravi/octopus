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

! ---------------------------------------------------------
subroutine X(eigensolver_mg) (der, st, hm, sdiag, niter, ik, diff)
  type(derivatives_t),    intent(in)    :: der
  type(states_t),         intent(inout) :: st
  type(hamiltonian_t),    intent(in)    :: hm
  type(subspace_t),       intent(in)    :: sdiag
  integer,                intent(inout) :: niter
  integer,                intent(in)    :: ik
  FLOAT,                  intent(out)   :: diff(:) !< (1:st%nst)

  integer  :: ist, ist2, iter
  R_TYPE, allocatable :: cc(:, :), aa(:)

  PUSH_SUB(X(eigensolver_mg))

  SAFE_ALLOCATE(cc(1:st%nst, 1:st%nst))
  SAFE_ALLOCATE(aa(1:st%nst))

  cc = R_TOTYPE(M_ZERO)

  do iter = 1, niter

    call X(subspace_diag)(sdiag, der, st, hm, ik, st%eigenval(:, ik), diff)

    if(in_debug_mode) then
      do ist = 1, st%nst
        write(message(1), '(a,i4,a,i4,a,i4,a,es12.6)') &
          'Debug: MG Eigensolver - ik', ik, ' ist ', ist, ' iter ', iter, ' res ', diff(ist)
        call messages_info(1)
      end do
    end if

    do ist = 1, st%nst

      aa(ist) = st%eigenval(ist, ik)

      cc(ist, ist) = M_ONE
      do ist2 = 1, ist - 1
        cc(ist, ist2) = X(mf_dotp)(der%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), st%X(psi)(:, :, ist2, ik))
      end do

    end do

    call X(coordinate_relaxation)(der, hm, st%nst, 10, ik, st%X(psi)(:, :, :, ik), aa, cc)

    ! normalize
    do ist = 1, st%nst      
      call lalg_scal(der%mesh%np, CNST(1.0)/sqrt(cc(ist, ist)), st%X(psi)(:, 1, ist, ik))
    end do

  end do

  call X(subspace_diag)(sdiag, der, st, hm, ik, st%eigenval(:, ik), diff)

  niter = iter*10

  SAFE_DEALLOCATE_A(cc)
  SAFE_DEALLOCATE_A(aa)

  POP_SUB(X(eigensolver_mg))
end subroutine X(eigensolver_mg)

! -------------------------------------------------------------------------

subroutine X(coordinate_relaxation)(der, hm, nst, steps, ik, psi, aa, cc)
  type(derivatives_t),    intent(in)    :: der
  type(hamiltonian_t),    intent(in)    :: hm
  integer,                intent(in)    :: nst
  integer,                intent(in)    :: steps
  integer,                intent(in)    :: ik
  R_TYPE,                 intent(inout) :: psi(:, :, :)
  R_TYPE,                 intent(inout) :: aa(:)
  R_TYPE,                 intent(inout) :: cc(:, :)

  integer :: ip, ist, ist2, iter, is, inb
  FLOAT, parameter :: ww = 5.0
  FLOAT :: pot, dh, vv
  R_TYPE  :: rta, rtb, rtc, alpha, s1, s2, s3
  
  R_TYPE, allocatable :: sigma(:), beta(:), hdiag(:, :)

  PUSH_SUB(X(coordinate_relaxation))

  SAFE_ALLOCATE(sigma(1:nst))
  SAFE_ALLOCATE(beta(1:nst))
  SAFE_ALLOCATE(hdiag(1:der%mesh%np, 1:hm%d%dim))

  call X(hamiltonian_diagonal)(hm, der, hdiag, ik)

  do iter = 1, steps
    
    do ip = 1, der%mesh%np
      
      if(der%mesh%use_curvilinear) then
        vv = sqrt(der%mesh%vol_pp(ip))
      else
        vv = sqrt(der%mesh%volume_element)
      endif

      dh = real(hdiag(ip, 1), REAL_PRECISION)
      pot = hm%vhxc(ip, 1) + hm%ep%vpsl(ip)
      
      do ist = 1, nst
        
        beta(ist) = psi(ip, 1, ist)*vv
        
        ! apply the Hamiltonian in the point
            
        alpha = M_ZERO
        do is = 1, der%lapl%stencil%size
          inb = ip + der%lapl%ri(is, der%lapl%rimap(ip))
          alpha = alpha + der%lapl%w_re(is, 1)*psi(inb, 1, ist)
        end do
        
        alpha = -M_HALF*alpha + pot*psi(ip, 1, ist)
        
        alpha = alpha*vv
        
        s1 = M_ZERO
        s2 = M_ZERO
        s3 = M_ZERO
        do ist2 = 1, ist - 1
          s1 = s1 + ww/cc(ist2, ist2)*beta(ist2)**2
          s2 = s2 + ww/cc(ist2, ist2)*beta(ist2)*cc(ist, ist2)
          s3 = s3 + ww/cc(ist2, ist2)*cc(ist, ist2)**2
        end do
        
        rta = alpha - dh*beta(ist) + s2 - beta(ist)*s1
        rtb = cc(ist, ist)*dh - aa(ist) + cc(ist, ist)*s1 - s3
        rtc = aa(ist)*beta(ist) - cc(ist, ist)*alpha + beta(ist)*s3 - s2
        
        sigma(ist) = -M_TWO*rtc/(rtb + sqrt(rtb**2 - M_FOUR*rta*rtc))
        
        psi(ip, 1, ist) = psi(ip, 1, ist) - sigma(ist)/vv
        
        aa(ist) = aa(ist) - M_TWO*sigma(ist)*alpha + sigma(ist)**2*dh
        cc(ist, ist) = cc(ist, ist) - M_TWO*sigma(ist)*beta(ist)  + sigma(ist)**2
        
        do ist2 = 1, ist - 1
          cc(ist, ist2) = cc(ist, ist2) - sigma(ist)*beta(ist2) - sigma(ist2)*beta(ist) + sigma(ist)*sigma(ist2)
        end do
        
      end do
      
    end do
    
  end do
  
  SAFE_DEALLOCATE_A(sigma)
  SAFE_DEALLOCATE_A(beta)
  SAFE_DEALLOCATE_A(hdiag)

  POP_SUB(X(coordinate_relaxation))
end subroutine X(coordinate_relaxation)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
