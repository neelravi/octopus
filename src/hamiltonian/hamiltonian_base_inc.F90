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

subroutine X(hamiltonian_base_local)(this, mesh, std, ispin, psib, vpsib)
  type(hamiltonian_base_t),    intent(in)    :: this
  type(mesh_t),                intent(in)    :: mesh
  type(states_dim_t),          intent(in)    :: std
  integer,                     intent(in)    :: ispin
  type(batch_t),               intent(in)    :: psib
  type(batch_t),               intent(inout) :: vpsib

  PUSH_SUB(X(hamiltonian_base_local))

  if(batch_status(psib) == BATCH_CL_PACKED) then
#ifdef HAVE_OPENCL
    call X(hamiltonian_base_local_sub)(this%potential, mesh, std, ispin, &
      psib, vpsib, potential_opencl = this%potential_opencl)
#endif
  else
    if(associated(this%Impotential)) then
      call X(hamiltonian_base_local_sub)(this%potential, mesh, std, ispin, &
        psib, vpsib, Impotential = this%Impotential)
    else
      call X(hamiltonian_base_local_sub)(this%potential, mesh, std, ispin, psib, vpsib)
    endif
  endif

  POP_SUB(X(hamiltonian_base_local))
end subroutine X(hamiltonian_base_local)

subroutine X(hamiltonian_base_local_sub)(potential, mesh, std, ispin, psib, vpsib, Impotential, potential_opencl)
  FLOAT,                        intent(in)    :: potential(:,:)
  type(mesh_t),                 intent(in)    :: mesh
  type(states_dim_t),           intent(in)    :: std
  integer,                      intent(in)    :: ispin
  type(batch_t), target,        intent(in)    :: psib
  type(batch_t), target,        intent(inout) :: vpsib
  FLOAT, optional,              intent(in)    :: Impotential(:,:)
  type(opencl_mem_t), optional, intent(in)    :: potential_opencl

  integer :: ist, ip
  R_TYPE, pointer :: psi(:, :), vpsi(:, :)
  R_TYPE  :: psi1, psi2
  FLOAT   :: vv, Imvv
  logical :: cmplxscl
#ifdef HAVE_OPENCL
  integer :: pnp, iprange
#endif

  call profiling_in(prof_vlpsi, "VLPSI")
  PUSH_SUB(X(hamiltonian_base_local_sub))

  cmplxscl = .false.
  if(present(Impotential)) cmplxscl = .true.

  if(batch_is_packed(psib) .or. batch_is_packed(vpsib)) then
    ASSERT(batch_is_packed(psib))
    ASSERT(batch_is_packed(vpsib))
    call batch_pack_was_modified(vpsib)
  end if

  select case(batch_status(psib))
  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
    ASSERT(.not. cmplxscl) ! not implemented

    pnp = opencl_padded_size(mesh%np)

    select case(std%ispin)

    case(UNPOLARIZED, SPIN_POLARIZED)
      call opencl_set_kernel_arg(kernel_vpsi, 0, pnp*(ispin - 1))
      call opencl_set_kernel_arg(kernel_vpsi, 1, mesh%np)
      call opencl_set_kernel_arg(kernel_vpsi, 2, potential_opencl)
      call opencl_set_kernel_arg(kernel_vpsi, 3, psib%pack%buffer)
      call opencl_set_kernel_arg(kernel_vpsi, 4, log2(psib%pack%size_real(1)))
      call opencl_set_kernel_arg(kernel_vpsi, 5, vpsib%pack%buffer)
      call opencl_set_kernel_arg(kernel_vpsi, 6, log2(vpsib%pack%size_real(1)))

      iprange = opencl_max_workgroup_size()/psib%pack%size_real(1)

      call opencl_kernel_run(kernel_vpsi, (/psib%pack%size_real(1), pnp/), (/psib%pack%size_real(1), iprange/))

    case(SPINORS)
      call opencl_set_kernel_arg(kernel_vpsi_spinors, 0, mesh%np)
      call opencl_set_kernel_arg(kernel_vpsi_spinors, 1, potential_opencl)
      call opencl_set_kernel_arg(kernel_vpsi_spinors, 2, pnp)
      call opencl_set_kernel_arg(kernel_vpsi_spinors, 3, psib%pack%buffer)
      call opencl_set_kernel_arg(kernel_vpsi_spinors, 4, psib%pack%size(1))
      call opencl_set_kernel_arg(kernel_vpsi_spinors, 5, vpsib%pack%buffer)
      call opencl_set_kernel_arg(kernel_vpsi_spinors, 6, vpsib%pack%size(1))

      call opencl_kernel_run(kernel_vpsi_spinors, (/psib%pack%size(1)/2, pnp/), &
        (/psib%pack%size(1)/2, 2*opencl_max_workgroup_size()/psib%pack%size(1)/))

    end select

    call opencl_finish()

    call profiling_count_operations((R_MUL*psib%nst)*mesh%np)
    call profiling_count_transfers(mesh%np, M_ONE)
    call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))

#endif
  case(BATCH_PACKED)

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      if(cmplxscl)then
        do ip = 1, mesh%np
          vv = potential(ip, ispin)
          Imvv = Impotential(ip, ispin)
          forall (ist = 1:psib%nst_linear)
            vpsib%pack%X(psi)(ist, ip) = vpsib%pack%X(psi)(ist, ip) + (vv+M_zI*Imvv)*psib%pack%X(psi)(ist, ip)
          end forall
        end do
      else
        !$omp parallel do private(vv, ist)
        do ip = 1, mesh%np
          vv = potential(ip, ispin)
          forall (ist = 1:psib%nst_linear)
            vpsib%pack%X(psi)(ist, ip) = vpsib%pack%X(psi)(ist, ip) + vv*psib%pack%X(psi)(ist, ip)
          end forall
        end do
        !$omp end parallel do
      end if
      call profiling_count_operations((2*R_ADD*psib%nst_linear)*mesh%np)
      call profiling_count_transfers(mesh%np, M_ONE)
      call profiling_count_transfers(mesh%np*psib%nst_linear, R_TOTYPE(M_ONE))

    case(SPINORS)
      ASSERT(mod(psib%nst_linear, 2) == 0)
      !the spinor case is more complicated since it mixes the two components.

      !$omp parallel do private(psi1, psi2, ist)
      do ip = 1, mesh%np
        do ist = 1, psib%nst_linear, 2
          psi1 = psib%pack%zpsi(ist    , ip)
          psi2 = psib%pack%zpsi(ist + 1, ip)
          vpsib%pack%zpsi(ist    , ip) = vpsib%pack%zpsi(ist    , ip) + &
            potential(ip, 1)*psi1 + (potential(ip, 3) + M_zI*potential(ip, 4))*psi2
          vpsib%pack%zpsi(ist + 1, ip) = vpsib%pack%zpsi(ist + 1, ip) + &
            potential(ip, 2)*psi2 + (potential(ip, 3) - M_zI*potential(ip, 4))*psi1
        end do
      end do
      !$omp end parallel do

      call profiling_count_operations((6*R_ADD + 2*R_MUL)*mesh%np*psib%nst)

    end select

  case(BATCH_NOT_PACKED)

    select case(std%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      if(cmplxscl)then
        do ist = 1, psib%nst
          forall (ip = 1:mesh%np)
            vpsib%states(ist)%X(psi)(ip, 1) = vpsib%states(ist)%X(psi)(ip, 1) + &
              (potential(ip, ispin)+ M_zI*Impotential(ip, ispin)) * psib%states(ist)%X(psi)(ip, 1) 
          end forall
        end do
      else
        !$omp parallel do private(ip)
        do ist = 1, psib%nst
          forall (ip = 1:mesh%np)
            vpsib%states(ist)%X(psi)(ip, 1) = vpsib%states(ist)%X(psi)(ip, 1) + &
              potential(ip, ispin) * psib%states(ist)%X(psi)(ip, 1)
          end forall
        end do
        !$omp end parallel do
      end if

      call profiling_count_operations((2*R_ADD*psib%nst)*mesh%np)
      call profiling_count_transfers(mesh%np, M_ONE)
      call profiling_count_transfers(mesh%np*psib%nst, R_TOTYPE(M_ONE))

    case(SPINORS)
      !the spinor case is more complicated since it mixes the two components.
      do ist = 1, psib%nst
        psi  => psib%states(ist)%X(psi)
        vpsi => vpsib%states(ist)%X(psi)

        forall(ip = 1:mesh%np)
          vpsi(ip, 1) = vpsi(ip, 1) + potential(ip, 1)*psi(ip, 1) + &
            (potential(ip, 3) + M_zI*potential(ip, 4))*psi(ip, 2)
          vpsi(ip, 2) = vpsi(ip, 2) + potential(ip, 2)*psi(ip, 2) + &
            (potential(ip, 3) - M_zI*potential(ip, 4))*psi(ip, 1)
        end forall

      end do
      call profiling_count_operations((6*R_ADD + 2*R_MUL)*mesh%np*psib%nst)

    end select

  end select

  call profiling_out(prof_vlpsi)
  POP_SUB(X(hamiltonian_base_local_sub))

end subroutine X(hamiltonian_base_local_sub)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_rashba)(this, der, std, psib, vpsib)
  type(hamiltonian_base_t),    intent(in)    :: this
  type(derivatives_t),         intent(in)    :: der
  type(states_dim_t),          intent(in)    :: std
  type(batch_t), target,       intent(in)    :: psib
  type(batch_t), target,       intent(inout) :: vpsib

  integer :: ist, idim, ip
  R_TYPE, pointer :: psi(:, :), vpsi(:, :)
  R_TYPE, allocatable :: grad(:, :, :)
  PUSH_SUB(X(hamiltonian_base_rashba))

  if(abs(this%rashba_coupling) < M_EPSILON) then
    POP_SUB(X(hamiltonian_base_rashba))
    return
  end if
  ASSERT(std%ispin == SPINORS)

  SAFE_ALLOCATE(grad(1:der%mesh%np, 1:MAX_DIM, 1:std%dim))

  do ist = 1, psib%nst
    psi  => psib%states(ist)%X(psi)
    vpsi => vpsib%states(ist)%X(psi)

    grad = M_ZERO
    do idim = 1, std%dim
      call X(derivatives_grad)(der, psi(:, idim), grad(:, :, idim), ghost_update = .false., set_bc = .false.)
    end do
    grad(:, der%mesh%sb%dim + 1:MAX_DIM, :) = M_ZERO
 
    if(associated(this%vector_potential)) then
      forall(ip = 1:der%mesh%np)
        vpsi(ip, 1) = vpsi(ip, 1) + &
          (this%rashba_coupling) * (this%vector_potential(2, ip) + M_zI * this%vector_potential(1, ip)) * psi(ip, 2)
        vpsi(ip, 2) = vpsi(ip, 2) + &
          (this%rashba_coupling) * (this%vector_potential(2, ip) - M_zI * this%vector_potential(1, ip)) * psi(ip, 1)
      end forall
    end if

    forall(ip = 1:der%mesh%np)
      vpsi(ip, 1) = vpsi(ip, 1) - &
        this%rashba_coupling*( grad(ip, 1, 2) - M_zI*grad(ip, 2, 2) )
      vpsi(ip, 2) = vpsi(ip, 2) + &
        this%rashba_coupling*( grad(ip, 1, 1) + M_zI*grad(ip, 2, 1) )
    end forall

  end do
  
  SAFE_DEALLOCATE_A(grad)
  
  POP_SUB(X(hamiltonian_base_rashba))
end subroutine X(hamiltonian_base_rashba)


subroutine X(hamiltonian_base_magnetic)(this, der, std, ep, ispin, psib, vpsib)
  type(hamiltonian_base_t),    intent(in)    :: this
  type(derivatives_t),         intent(in)    :: der
  type(states_dim_t),          intent(in)    :: std
  type(epot_t),                intent(in)    :: ep
  integer,                     intent(in)    :: ispin
  type(batch_t), target,       intent(in)    :: psib
  type(batch_t), target,       intent(inout) :: vpsib

  integer :: ist, idim, ip
  R_TYPE, pointer :: psi(:, :), vpsi(:, :)
  R_TYPE, allocatable :: grad(:, :, :)
  FLOAT :: cc, b2, bb(1:MAX_DIM)
  CMPLX :: b12

  if(.not. hamiltonian_base_has_magnetic(this)) return

  call profiling_in(prof_magnetic, "MAGNETIC")
  PUSH_SUB(X(hamiltonian_base_magnetic))

  SAFE_ALLOCATE(grad(1:der%mesh%np, 1:MAX_DIM, 1:std%dim))

  do ist = 1, psib%nst
    psi  => psib%states(ist)%X(psi)
    vpsi => vpsib%states(ist)%X(psi)

    do idim = 1, std%dim
      call X(derivatives_grad)(der, psi(:, idim), grad(:, :, idim), ghost_update = .false., set_bc = .false.)
    end do
    grad(:, der%mesh%sb%dim + 1:MAX_DIM, :) = M_ZERO
 
    if(associated(this%vector_potential)) then
      forall (idim = 1:std%dim, ip = 1:der%mesh%np)
        vpsi(ip, idim) = vpsi(ip, idim) + (M_HALF / this%mass) * sum(this%vector_potential(1:MAX_DIM, ip)**2)*psi(ip, idim) &
          + (M_ONE / this%mass) * M_zI*dot_product(this%vector_potential(1:MAX_DIM, ip), grad(ip, 1:MAX_DIM, idim))
      end forall
    end if

    if(associated(this%uniform_magnetic_field).and. std%ispin /= UNPOLARIZED) then
      ! Zeeman term
      cc = M_HALF/P_C*ep%gyromagnetic_ratio*M_HALF
      bb = this%uniform_magnetic_field
      b2 = sqrt(sum(this%uniform_magnetic_field**2))
      b12 = bb(1) - M_ZI*bb(2)

      select case (std%ispin)
      case (SPIN_POLARIZED)
        if(is_spin_down(ispin)) cc = -cc
        
        forall (ip = 1:der%mesh%np)
          vpsi(ip, 1) = vpsi(ip, 1) + cc*b2*psi(ip, 1)
        end forall
        
      case (SPINORS)
        forall (ip = 1:der%mesh%np)
          vpsi(ip, 1) = vpsi(ip, 1) + cc*(bb(3)*psi(ip, 1) + b12*psi(ip, 2))
          vpsi(ip, 2) = vpsi(ip, 2) + cc*(-bb(3)*psi(ip, 2) + conjg(b12)*psi(ip, 1))
        end forall
        
      end select
    end if
  end do
  
  SAFE_DEALLOCATE_A(grad)
  
  POP_SUB(X(hamiltonian_base_magnetic))
  call profiling_out(prof_magnetic)
end subroutine X(hamiltonian_base_magnetic)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_nlocal_start)(this, mesh, std, ik, psib, projection)
  type(hamiltonian_base_t), target, intent(in)    :: this
  type(mesh_t),                     intent(in)    :: mesh
  type(states_dim_t),               intent(in)    :: std
  integer,                          intent(in)    :: ik
  type(batch_t),                    intent(in)    :: psib
  type(projection_t),               intent(out)   :: projection

  integer :: ist, ip, iproj, imat, nreal, iprojection
  integer :: npoints, nprojs, nst
  R_TYPE, allocatable :: psi(:, :)
  R_TYPE :: phase
  type(projector_matrix_t), pointer :: pmat
#ifdef HAVE_OPENCL
  integer :: padnprojs, wgsize, lnprojs, size
  type(profile_t), save :: cl_prof
  type(octcl_kernel_t), save :: ker_proj_bra, ker_proj_bra_phase
  type(cl_kernel) :: kernel
#endif
  if(.not. this%apply_projector_matrices) return

  call profiling_in(prof_vnlpsi_start, "VNLPSI_MAT_BRA")
  PUSH_SUB(X(hamiltonian_base_nlocal_start))

  nst = psib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst
#else
  nreal = nst
#endif

#ifdef HAVE_OPENCL
  if(batch_is_packed(psib) .and. opencl_is_enabled()) then
   
    call opencl_create_buffer(projection%buff_projection, CL_MEM_READ_WRITE, R_TYPE_VAL, &
      this%full_projection_size*psib%pack%size_real(1))
    
    call profiling_in(cl_prof, "CL_PROJ_BRA")

    if(associated(this%projector_phases)) then
      call octcl_kernel_start_call(ker_proj_bra_phase, 'projector.cl', 'projector_bra_phase')
      kernel = octcl_kernel_get_ref(ker_proj_bra_phase)
      size = psib%pack%size(1)
      ASSERT(R_TYPE_VAL == TYPE_CMPLX)
    else
      call octcl_kernel_start_call(ker_proj_bra, 'projector.cl', 'projector_bra')
      kernel = octcl_kernel_get_ref(ker_proj_bra)
      size = psib%pack%size_real(1)
    end if

    call opencl_set_kernel_arg(kernel, 0, this%nprojector_matrices)
    call opencl_set_kernel_arg(kernel, 1, this%buff_offsets)
    call opencl_set_kernel_arg(kernel, 2, this%buff_matrices)
    call opencl_set_kernel_arg(kernel, 3, this%buff_maps)
    call opencl_set_kernel_arg(kernel, 4, this%buff_scals)
    call opencl_set_kernel_arg(kernel, 5, psib%pack%buffer)
    call opencl_set_kernel_arg(kernel, 6, log2(size))
    call opencl_set_kernel_arg(kernel, 7, projection%buff_projection)
    call opencl_set_kernel_arg(kernel, 8, log2(size))

    if(associated(this%projector_phases)) then
      call opencl_set_kernel_arg(kernel, 9, this%buff_projector_phases)
      call opencl_set_kernel_arg(kernel, 10, (ik - std%kpt%start)*this%total_points)
    end if

    padnprojs = pad_pow2(this%max_nprojs)
    lnprojs = min(opencl_kernel_workgroup_size(kernel)/size, padnprojs)

    call opencl_kernel_run(kernel, &
      (/size, padnprojs, this%nprojector_matrices/), (/size, lnprojs, 1/))

    do imat = 1, this%nprojector_matrices
      pmat => this%projector_matrices(imat)
      
      npoints = pmat%npoints
      nprojs = pmat%nprojs
      
      call profiling_count_operations(nreal*nprojs*M_TWO*npoints + nst*nprojs)
    end do

    call opencl_finish()

    if(mesh%parallel_in_domains) then
      SAFE_ALLOCATE(projection%X(projection)(1:psib%pack%size_real(1), 1:this%full_projection_size))
      call opencl_read_buffer(projection%buff_projection, &
        this%full_projection_size*psib%pack%size_real(1), projection%X(projection))
    end if

    call profiling_out(cl_prof)

    POP_SUB(X(hamiltonian_base_nlocal_start))
    call profiling_out(prof_vnlpsi_start)
    return
  end if
#endif

  SAFE_ALLOCATE(projection%X(projection)(1:nst, 1:this%full_projection_size))
  projection%X(projection) = M_ZERO
  iprojection = 0

  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if(npoints /= 0) then

      SAFE_ALLOCATE(psi(1:nst, 1:npoints))

      call profiling_in(prof_gather, "PROJ_MAT_GATHER")

      ! collect all the points we need in a continuous array
      if(batch_is_packed(psib)) then

        !$omp parallel do private(ip, ist)
        do ip = 1, npoints
          forall(ist = 1:nst)
            psi(ist, ip) = psib%pack%X(psi)(ist, pmat%map(ip))
          end forall
        end do
        !$omp end parallel do

      else

        do ist = 1, nst
          !$omp parallel do
          do ip = 1, npoints
            psi(ist, ip) = psib%states_linear(ist)%X(psi)(pmat%map(ip))
          end do
          !$omp end parallel do 
        end do

      end if

      if(associated(this%projector_phases)) then

        !$omp parallel do private(ip, ist, phase)
        do ip = 1, npoints
          phase = this%projector_phases(ip, imat, ik)
          forall(ist = 1:nst)
            psi(ist, ip) = phase*psi(ist, ip)
          end forall
        end do
        !$omp end parallel do

      end if

      call profiling_out(prof_gather)
      
      ! Now matrix-multiply to calculate the projections.
      ! the line below does: projection = matmul(psi, pmat%projectors)
      call blas_gemm('N', 'N', nreal, nprojs, npoints, M_ONE, psi(1, 1), nreal, pmat%projectors(1, 1), npoints, &
        M_ZERO, projection%X(projection)(1, iprojection + 1), nreal)

      ! apply the scale
      !$omp parallel do private(iproj, ist)
      do iproj = 1, nprojs
        do ist = 1, nst
          projection%X(projection)(ist, iprojection + iproj) = projection%X(projection)(ist, iprojection + iproj)*pmat%scal(iproj)
        end do
      end do
      !$omp end parallel do

      call profiling_count_operations(nreal*nprojs*M_TWO*npoints + nst*nprojs)

    end if

    SAFE_DEALLOCATE_A(psi)

    INCR(iprojection, nprojs)
  end do

  POP_SUB(X(hamiltonian_base_nlocal_start))
  call profiling_out(prof_vnlpsi_start)
end subroutine X(hamiltonian_base_nlocal_start)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_nlocal_finish)(this, mesh, std, ik, projection, vpsib)
  type(hamiltonian_base_t), target, intent(in)    :: this
  type(mesh_t),                     intent(in)    :: mesh
  type(states_dim_t),               intent(in)    :: std
  integer,                          intent(in)    :: ik
  type(projection_t),               intent(inout) :: projection
  type(batch_t),                    intent(inout) :: vpsib

  integer :: ist, ip, imat, nreal, iprojection
  integer :: npoints, nprojs, nst
  R_TYPE  :: phase
  R_TYPE, allocatable :: psi(:, :)
  type(projector_matrix_t), pointer :: pmat
#ifdef HAVE_MPI
  integer :: d1
  R_TYPE, allocatable :: projection_red(:, :)
  type(profile_t), save :: reduce_prof
#endif

  if(.not. this%apply_projector_matrices) return

  call profiling_in(prof_vnlpsi_finish, "VNLPSI_MAT_KET")
  PUSH_SUB(X(hamiltonian_base_nlocal_finish))

  nst = vpsib%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst
#else
  nreal = nst
#endif

  ! reduce the projections
#ifdef HAVE_MPI
  if(mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, "VNLPSI_MAT_REDUCE")
    d1 = ubound(projection%X(projection), dim = 1)
    SAFE_ALLOCATE(projection_red(1:d1, 1:this%full_projection_size))
    projection_red = projection%X(projection)
    call MPI_Allreduce(projection_red(1, 1), projection%X(projection)(1, 1), d1*this%full_projection_size, &
      R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
    SAFE_DEALLOCATE_A(projection_red)
    call profiling_out(reduce_prof)
  end if
#endif

#ifdef HAVE_OPENCL
  if(batch_is_packed(vpsib) .and. opencl_is_enabled()) then

    if(mesh%parallel_in_domains) then
      call opencl_write_buffer(projection%buff_projection, &
        this%full_projection_size*vpsib%pack%size_real(1), projection%X(projection))
      SAFE_DEALLOCATE_P(projection%X(projection))
    end if

    call finish_opencl()
    call opencl_release_buffer(projection%buff_projection)

    POP_SUB(X(hamiltonian_base_nlocal_finish))
    call profiling_out(prof_vnlpsi_finish)
    return
  end if
#endif

  ASSERT(associated(projection%X(projection)))

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs

    if(npoints /=  0) then

      SAFE_ALLOCATE(psi(1:nst, 1:npoints))

      ! Matrix-multiply again.
      ! the line below does: psi = matmul(projection, transpose(pmat%projectors))
      call blas_gemm('N', 'T', nreal, npoints, nprojs, &
        M_ONE, projection%X(projection)(1, iprojection + 1), nreal, pmat%projectors(1, 1), npoints, &
        M_ZERO, psi(1, 1), nreal)
      
      call profiling_count_operations(nreal*nprojs*M_TWO*npoints)

      call profiling_in(prof_scatter, "PROJ_MAT_SCATTER")

      if(associated(this%projector_phases)) then
        !$omp parallel do private(ip, ist, phase)
        do ip = 1, npoints
          phase = conjg(this%projector_phases(ip, imat, ik))
          forall(ist = 1:nst)
            psi(ist, ip) = phase*psi(ist, ip)
          end forall
        end do
        !$omp end parallel do
      end if
      
      ! and copy the points from the local buffer to its position
      if(batch_is_packed(vpsib)) then
        !$omp parallel do private(ip, ist)
        do ip = 1, npoints
          forall(ist = 1:nst)
            vpsib%pack%X(psi)(ist, pmat%map(ip)) = vpsib%pack%X(psi)(ist, pmat%map(ip)) + psi(ist, ip)
          end forall
        end do
        !$omp end parallel do

        call batch_pack_was_modified(vpsib)
      else
        
        do ist = 1, nst
          !$omp parallel do
          do ip = 1, npoints
            vpsib%states_linear(ist)%X(psi)(pmat%map(ip)) = vpsib%states_linear(ist)%X(psi)(pmat%map(ip)) + psi(ist, ip)
          end do
          !$omp end parallel do
        end do
        
      end if
      call profiling_count_operations(nst*npoints*R_ADD)
      call profiling_out(prof_scatter)
    end if
    
    SAFE_DEALLOCATE_A(psi)
    
    INCR(iprojection, nprojs)
  end do
  
  SAFE_DEALLOCATE_P(projection%X(projection))
  
  POP_SUB(X(hamiltonian_base_nlocal_finish))
  call profiling_out(prof_vnlpsi_finish)

contains

  subroutine finish_opencl()
#ifdef HAVE_OPENCL
    integer :: wgsize, imat, iregion, size
    type(profile_t), save :: cl_prof
    type(octcl_kernel_t), save :: ker_proj_ket, ker_proj_ket_phase
    type(cl_kernel) :: kernel

    PUSH_SUB(X(hamiltonian_base_nlocal_finish).finish_opencl)

    ! In this case we run one kernel per projector, since all write to
    ! the wave-function. Otherwise we would need to do atomic
    ! operations.

    call profiling_in(cl_prof, "CL_PROJ_KET")

    if(associated(this%projector_phases)) then
      call octcl_kernel_start_call(ker_proj_ket_phase, 'projector.cl', 'projector_ket_phase')
      kernel = octcl_kernel_get_ref(ker_proj_ket_phase)
      size = vpsib%pack%size(1)
      ASSERT(R_TYPE_VAL == TYPE_CMPLX)
    else
      call octcl_kernel_start_call(ker_proj_ket, 'projector.cl', 'projector_ket')
      kernel = octcl_kernel_get_ref(ker_proj_ket)
      size = vpsib%pack%size_real(1)
    end if

    do iregion = 1, this%nregions
      
      call opencl_set_kernel_arg(kernel, 0, this%nprojector_matrices)
      call opencl_set_kernel_arg(kernel, 1, this%regions(iregion) - 1)
      call opencl_set_kernel_arg(kernel, 2, this%buff_offsets)
      call opencl_set_kernel_arg(kernel, 3, this%buff_matrices)
      call opencl_set_kernel_arg(kernel, 4, this%buff_maps)
      call opencl_set_kernel_arg(kernel, 5, projection%buff_projection)
      call opencl_set_kernel_arg(kernel, 6, log2(size))
      call opencl_set_kernel_arg(kernel, 7, vpsib%pack%buffer)
      call opencl_set_kernel_arg(kernel, 8, log2(size))

      if(associated(this%projector_phases)) then
        call opencl_set_kernel_arg(kernel, 9, this%buff_projector_phases)
        call opencl_set_kernel_arg(kernel, 10, (ik - std%kpt%start)*this%total_points)
      end if

      wgsize = opencl_kernel_workgroup_size(kernel)/size    

      call opencl_kernel_run(kernel, &
        (/size, pad(this%max_npoints, wgsize), this%regions(iregion + 1) - this%regions(iregion)/), &
        (/size, wgsize, 1/))
      
      call opencl_finish()
      
    end do
    
    do imat = 1, this%nprojector_matrices
      pmat => this%projector_matrices(imat)
      npoints = pmat%npoints
      nprojs = pmat%nprojs
      call profiling_count_operations(nreal*nprojs*M_TWO*npoints)
      call profiling_count_operations(nst*npoints*R_ADD)
    end do

    call batch_pack_was_modified(vpsib)
    call opencl_finish()

    call profiling_out(cl_prof)

    POP_SUB(X(hamiltonian_base_nlocal_finish).finish_opencl)
#endif
  end subroutine finish_opencl

end subroutine X(hamiltonian_base_nlocal_finish)

! ---------------------------------------------------------------------------------------

subroutine X(hamiltonian_base_nlocal_force)(this, mesh, st, geo, iqn, ndim, psi1b, psi2b, force)
  type(hamiltonian_base_t), target, intent(in)    :: this
  type(mesh_t),                     intent(in)    :: mesh
  type(states_t),                   intent(in)    :: st
  type(geometry_t),                 intent(in)    :: geo
  integer,                          intent(in)    :: iqn
  integer,                          intent(in)    :: ndim
  type(batch_t),                    intent(in)    :: psi1b
  type(batch_t),                    intent(in)    :: psi2b(:)
  FLOAT,                            intent(inout) :: force(:, :)

  integer :: ii, ist, ip, iproj, imat, nreal, iprojection, iatom, idir
  integer :: npoints, nprojs, nst
  R_TYPE, allocatable :: psi(:, :, :), projs(:, :, :), ff(:)
  type(projector_matrix_t), pointer :: pmat

  if(.not. this%apply_projector_matrices) return

  call profiling_in(prof_vnlpsi_start, "VNLPSI_MAT_ELEMENT")
  PUSH_SUB(X(hamiltonian_base_nlocal_force))

  ASSERT(psi1b%nst_linear == psi2b(1)%nst_linear)
  ASSERT(batch_status(psi1b) == batch_status(psi2b(1)))

  nst = psi1b%nst_linear
#ifdef R_TCOMPLEX
  nreal = 2*nst
#else
  nreal = nst
#endif

  iprojection = 0
  do imat = 1, this%nprojector_matrices
    pmat => this%projector_matrices(imat)

    npoints = pmat%npoints
    nprojs = pmat%nprojs

    SAFE_ALLOCATE(projs(0:ndim, 1:nst, 1:nprojs))

    if(npoints /= 0) then

      SAFE_ALLOCATE(psi(0:ndim, 1:nst, 1:npoints))
      
      call profiling_in(prof_gather, "PROJ_MAT_ELEM_GATHER")

      ! collect all the points we need in a continuous array
      if(batch_is_packed(psi1b)) then
        forall(ip = 1:npoints)
          forall(ist = 1:nst)
            psi(0, ist, ip) = psi1b%pack%X(psi)(ist, pmat%map(ip))
            forall(idir = 1:ndim) psi(idir, ist, ip) = psi2b(idir)%pack%X(psi)(ist, pmat%map(ip))
          end forall
        end forall
      else
        forall(ip = 1:npoints)
          forall(ist = 1:nst) 
            psi(0, ist, ip) = psi1b%states_linear(ist)%X(psi)(pmat%map(ip))
            forall(idir = 1:ndim) psi(idir, ist, ip) = psi2b(idir)%states_linear(ist)%X(psi)(pmat%map(ip))
          end forall
        end forall
      end if

      if(associated(this%projector_phases)) then
        forall(ip = 1:npoints)
          forall(ist = 1:nst)
            forall(idir = 0:ndim)
              psi(idir, ist, ip) = this%projector_phases(ip, imat, iqn)*psi(idir, ist, ip)
            end forall
          end forall
        end forall
      end if

      call profiling_out(prof_gather)
      
      ! Now matrix-multiply to calculate the projections. We can do all the matrix multiplications at once
      call blas_gemm('N', 'N', (ndim + 1)*nreal, nprojs, npoints, M_ONE, &
        psi(0, 1, 1), (ndim + 1)*nreal, pmat%projectors(1, 1), npoints, &
        M_ZERO, projs(0, 1, 1), (ndim + 1)*nreal)
      
      call profiling_count_operations(nreal*(ndim + 1)*nprojs*M_TWO*npoints)

    else
      
      projs = CNST(0.0)

    end if

    if(mesh%parallel_in_domains) call comm_allreduce(mesh%mpi_grp%comm, projs)
    
    iatom = this%projector_to_atom(imat)
    
    SAFE_ALLOCATE(ff(1:ndim))
    
    ff(1:ndim) = CNST(0.0)
    
    do ii = 1, psi1b%nst_linear
      ist = batch_linear_to_ist(psi1b, ii)
      do iproj = 1, nprojs
        do idir = 1, ndim
          ff(idir) = ff(idir) - CNST(2.0)*st%d%kweights(iqn)*st%occ(ist, iqn)*pmat%scal(iproj)*mesh%volume_element*&
            R_CONJ(projs(0, ii, iproj))*projs(idir, ii, iproj)
        end do
      end do
    end do
    
    force(1:ndim, iatom) = force(1:ndim, iatom) + ff(1:ndim)
    
    call profiling_count_operations((R_ADD + 2*R_MUL)*nst*ndim*nprojs)

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(projs)
    SAFE_DEALLOCATE_A(ff)

    INCR(iprojection, nprojs)
  end do

  POP_SUB(X(hamiltonian_base_nlocal_force))
  call profiling_out(prof_vnlpsi_start)
end subroutine X(hamiltonian_base_nlocal_force)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
