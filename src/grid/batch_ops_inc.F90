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

!--------------------------------------------------------------

subroutine X(batch_set)(this, np, psi)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: np
  R_TYPE,         intent(in)    :: psi(:, :, :)

  integer :: ist, idim

  PUSH_SUB(X(batch_set))

  do ist = 1, this%nst
    do idim = 1, this%dim
      call lalg_copy(np, psi(:, idim, ist), this%states(ist)%X(psi)(:, idim))
    end do
  end do

  POP_SUB(X(batch_set))
end subroutine X(batch_set)

! --------------------------------------------------------------

subroutine X(batch_axpy_const)(np, aa, xx, yy)
  integer,           intent(in)    :: np
  R_TYPE,            intent(in)    :: aa
  type(batch_t),     intent(in)    :: xx
  type(batch_t),     intent(inout) :: yy

  integer :: ist
#ifdef HAVE_OPENCL
  integer :: localsize
#endif
  CMPLX :: zaa

  PUSH_SUB(X(batch_axpy_const))
  call profiling_in(axpy_const_prof, "BATCH_AXPY_CONST")

  ASSERT(batch_type(yy) == batch_type(xx))
#ifdef R_TCMPLX
  !if aa is complex, the functions must be complex
  ASSERT(batch_type(yy) == TYPE_CMPLX)
#endif
  ASSERT(xx%nst_linear == yy%nst_linear)
  ASSERT(batch_status(xx) == batch_status(yy))

  select case(batch_status(xx))
  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
    if(batch_type(yy) == TYPE_FLOAT) then

      call opencl_set_kernel_arg(kernel_daxpy, 0, aa)
      call opencl_set_kernel_arg(kernel_daxpy, 1, xx%pack%buffer)
      call opencl_set_kernel_arg(kernel_daxpy, 2, log2(xx%pack%size(1)))
      call opencl_set_kernel_arg(kernel_daxpy, 3, yy%pack%buffer)
      call opencl_set_kernel_arg(kernel_daxpy, 4, log2(yy%pack%size(1)))
      
      localsize = opencl_max_workgroup_size()
      call opencl_kernel_run(kernel_daxpy, (/yy%pack%size(1), pad(np, localsize)/), (/yy%pack%size(1), localsize/yy%pack%size(1)/))
      
    else
      zaa = aa
      call opencl_set_kernel_arg(kernel_zaxpy, 0, zaa)
      call opencl_set_kernel_arg(kernel_zaxpy, 1, xx%pack%buffer)
      call opencl_set_kernel_arg(kernel_zaxpy, 2, log2(xx%pack%size(1)))
      call opencl_set_kernel_arg(kernel_zaxpy, 3, yy%pack%buffer)
      call opencl_set_kernel_arg(kernel_zaxpy, 4, log2(yy%pack%size(1)))
      
      localsize = opencl_max_workgroup_size()
      call opencl_kernel_run(kernel_zaxpy, (/yy%pack%size(1), pad(np, localsize)/), (/yy%pack%size(1), localsize/yy%pack%size(1)/))

    end if

    call opencl_finish()
#endif
  case(BATCH_PACKED)
    if(batch_type(yy) == TYPE_CMPLX) then
      call lalg_axpy(xx%pack%size(1), np, aa, xx%pack%zpsi, yy%pack%zpsi)
    else
#ifdef R_TREAL
      call lalg_axpy(xx%pack%size(1), np, aa, xx%pack%dpsi, yy%pack%dpsi)
#endif
    end if

  case(BATCH_NOT_PACKED)
    do ist = 1, yy%nst_linear
      if(batch_type(yy) == TYPE_CMPLX) then
        call lalg_axpy(np, aa, xx%states_linear(ist)%zpsi, yy%states_linear(ist)%zpsi)
      else
#ifdef R_TREAL
        call lalg_axpy(np, aa, xx%states_linear(ist)%dpsi, yy%states_linear(ist)%dpsi)
#endif
      end if
    end do
  end select

  call profiling_count_operations(xx%nst*np*(R_ADD + R_MUL)*types_get_size(batch_type(xx))/types_get_size(TYPE_FLOAT))

  call batch_pack_was_modified(yy)

  call profiling_out(axpy_const_prof)
  POP_SUB(X(batch_axpy_const))
end subroutine X(batch_axpy_const)

! --------------------------------------------------------------

subroutine X(batch_axpy_vec)(np, aa, xx, yy, a_start, a_full)
  integer,           intent(in)    :: np
  R_TYPE,            intent(in)    :: aa(:)
  type(batch_t),     intent(in)    :: xx
  type(batch_t),     intent(inout) :: yy
  integer, optional, intent(in)    :: a_start
  logical, optional, intent(in)    :: a_full

  integer :: ist, ip, effsize, iaa
  R_TYPE, allocatable     :: aa_linear(:)
#ifdef HAVE_OPENCL
  integer :: localsize
  integer :: size_factor
  type(opencl_mem_t)      :: aa_buffer
  FLOAT,  allocatable     :: aa_linear_double(:)
  type(octcl_kernel_t), save :: kernel
  type(cl_kernel)         :: kernel_ref
#endif
  
  PUSH_SUB(X(batch_axpy_vec))
  call profiling_in(axpy_vec_prof, "BATCH_AXPY_VEC")

  ASSERT(batch_type(yy) == batch_type(xx))
#ifdef R_TCMPLX
  !if aa is complex, the functions must be complex
  ASSERT(batch_type(yy) == TYPE_CMPLX)
#endif
  ASSERT(xx%nst_linear == yy%nst_linear)
  ASSERT(batch_status(xx) == batch_status(yy))

  effsize = yy%nst_linear
  if(batch_is_packed(yy)) effsize = yy%pack%size(1)
  SAFE_ALLOCATE(aa_linear(1:effsize))

  aa_linear = M_ZERO
  do ist = 1, yy%nst_linear
    iaa = batch_linear_to_ist(xx, ist) - (optional_default(a_start, 1) - 1)
    if(.not. optional_default(a_full, .true.)) iaa = iaa - (batch_linear_to_ist(xx, 1) - 1)
    aa_linear(ist) = aa(iaa)
  end do

  select case(batch_status(xx))
  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
    call octcl_kernel_start_call(kernel, 'axpy.cl', TOSTRING(X(axpy_vec)), flags = '-D'//R_TYPE_CL)
    kernel_ref = octcl_kernel_get_ref(kernel)

    if(batch_type(yy) == TYPE_CMPLX .and. R_TYPE_VAL == TYPE_FLOAT) then
      size_factor = 2
      SAFE_ALLOCATE(aa_linear_double(1:2*yy%pack%size(1)))
      do ist = 1, yy%pack%size(1)
        aa_linear_double(2*ist - 1) = aa_linear(ist)
        aa_linear_double(2*ist) = aa_linear(ist)
      end do
      call opencl_create_buffer(aa_buffer, CL_MEM_READ_ONLY, TYPE_FLOAT, 2*yy%pack%size(1))
      call opencl_write_buffer(aa_buffer, 2*yy%pack%size(1), aa_linear_double)
      SAFE_DEALLOCATE_A(aa_linear_double)
    else
      size_factor = 1
      call opencl_create_buffer(aa_buffer, CL_MEM_READ_ONLY, R_TYPE_VAL, yy%pack%size(1))
      call opencl_write_buffer(aa_buffer, yy%pack%size(1), aa_linear)
    end if

    call opencl_set_kernel_arg(kernel_ref, 0, aa_buffer)
    call opencl_set_kernel_arg(kernel_ref, 1, xx%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 2, log2(xx%pack%size(1)*size_factor))
    call opencl_set_kernel_arg(kernel_ref, 3, yy%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 4, log2(yy%pack%size(1)*size_factor))

    localsize = opencl_max_workgroup_size()
    call opencl_kernel_run(kernel_ref, (/yy%pack%size(1)*size_factor, pad(np, localsize)/), &
      (/yy%pack%size(1)*size_factor, localsize/(yy%pack%size(1)*size_factor)/))

    call opencl_finish()

    call opencl_release_buffer(aa_buffer)
#endif
  case(BATCH_PACKED)
    if(batch_type(yy) == TYPE_CMPLX) then
      !$omp parallel do private(ip, ist)
      do ip = 1, np
        do ist = 1, yy%pack%size(1)
          yy%pack%zpsi(ist, ip) = aa_linear(ist)*xx%pack%zpsi(ist, ip) + yy%pack%zpsi(ist, ip)
        end do
      end do
    else
#ifdef R_TREAL
      !$omp parallel do private(ip, ist)
      do ip = 1, np
        do ist = 1, yy%pack%size(1)
          yy%pack%dpsi(ist, ip) = aa_linear(ist)*xx%pack%dpsi(ist, ip) + yy%pack%dpsi(ist, ip)
        end do
      end do
#endif
    end if
    
  case(BATCH_NOT_PACKED)
    do ist = 1, yy%nst_linear
      if(batch_type(yy) == TYPE_CMPLX) then
        call lalg_axpy(np, aa_linear(ist), xx%states_linear(ist)%zpsi, yy%states_linear(ist)%zpsi)
      else
#ifdef R_TREAL
        call lalg_axpy(np, aa_linear(ist), xx%states_linear(ist)%dpsi, yy%states_linear(ist)%dpsi)
#endif
      end if
    end do
  end select

  call batch_pack_was_modified(yy)

  SAFE_DEALLOCATE_A(aa_linear)

  call profiling_count_operations(xx%nst*np*(R_ADD + R_MUL)*types_get_size(batch_type(xx))/types_get_size(TYPE_FLOAT))

  call profiling_out(axpy_vec_prof)
  POP_SUB(X(batch_axpy_vec))
end subroutine X(batch_axpy_vec)


! --------------------------------------------------------------

subroutine X(batch_scal_vec)(np, aa, xx, a_start, a_full)
  integer,           intent(in)    :: np
  R_TYPE,            intent(in)    :: aa(:)
  type(batch_t),     intent(inout) :: xx
  integer, optional, intent(in)    :: a_start
  logical, optional, intent(in)    :: a_full

  integer :: ist, ip, effsize, iaa
  R_TYPE, allocatable     :: aa_linear(:)
#ifdef HAVE_OPENCL
  integer :: localsize
  integer :: size_factor
  FLOAT,  allocatable     :: aa_linear_double(:)
  type(opencl_mem_t)      :: aa_buffer
  type(octcl_kernel_t), save :: kernel
  type(cl_kernel)         :: kernel_ref
#endif
  
  PUSH_SUB(X(batch_scal_vec))
  call profiling_in(scal_prof, "BATCH_SCAL")

#ifdef R_TCMPLX
  !if aa is complex, the functions must be complex
  ASSERT(batch_type(xx) == TYPE_CMPLX)
#endif

  effsize = xx%nst_linear
  if(batch_is_packed(xx)) effsize = xx%pack%size(1)
  SAFE_ALLOCATE(aa_linear(1:effsize))

  aa_linear = M_ZERO
  do ist = 1, xx%nst_linear
    iaa = batch_linear_to_ist(xx, ist) - (optional_default(a_start, 1) - 1)
    if(.not. optional_default(a_full, .true.)) iaa = iaa - (batch_linear_to_ist(xx, 1) - 1)
    aa_linear(ist) = aa(iaa)
  end do
  
  select case(batch_status(xx))
  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL

    if(batch_type(xx) == TYPE_CMPLX .and. R_TYPE_VAL == TYPE_FLOAT) then
      size_factor = 2
      SAFE_ALLOCATE(aa_linear_double(1:2*xx%pack%size(1)))
      do ist = 1, xx%pack%size(1)
        aa_linear_double(2*ist - 1) = aa_linear(ist)
        aa_linear_double(2*ist) = aa_linear(ist)
      end do
      call opencl_create_buffer(aa_buffer, CL_MEM_READ_ONLY, TYPE_FLOAT, 2*xx%pack%size(1))
      call opencl_write_buffer(aa_buffer, 2*xx%pack%size(1), aa_linear_double)
      SAFE_DEALLOCATE_A(aa_linear_double)
    else
      size_factor = 1
      call opencl_create_buffer(aa_buffer, CL_MEM_READ_ONLY, R_TYPE_VAL, xx%pack%size(1))
      call opencl_write_buffer(aa_buffer, xx%pack%size(1), aa_linear)
    end if

    call octcl_kernel_start_call(kernel, 'axpy.cl', TOSTRING(X(scal_vec)), flags = '-D'//R_TYPE_CL)
  
    kernel_ref = octcl_kernel_get_ref(kernel)

    call opencl_set_kernel_arg(kernel_ref, 0, aa_buffer)
    call opencl_set_kernel_arg(kernel_ref, 1, xx%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 2, log2(xx%pack%size(1)*size_factor))

    localsize = opencl_max_workgroup_size()
    call opencl_kernel_run(kernel_ref, (/xx%pack%size(1)*size_factor, pad(np, localsize)/), &
      (/xx%pack%size(1)*size_factor, localsize/(xx%pack%size(1)*size_factor)/))

    call opencl_finish()

    call opencl_release_buffer(aa_buffer)
#endif
  case(BATCH_PACKED)
    if(batch_type(xx) == TYPE_CMPLX) then
      do ist = 1, xx%pack%size(1)
        do ip = 1, np
          xx%pack%zpsi(ist, ip) = aa_linear(ist)*xx%pack%zpsi(ist, ip)
        end do
      end do
    else
#ifdef R_TREAL
      do ist = 1, xx%pack%size(1)
        do ip = 1, np
          xx%pack%dpsi(ist, ip) = aa_linear(ist)*xx%pack%dpsi(ist, ip)
        end do
      end do
#endif
    end if
    
  case(BATCH_NOT_PACKED)
    do ist = 1, xx%nst_linear
      if(batch_type(xx) == TYPE_CMPLX) then
        call lalg_scal(np, aa_linear(ist), xx%states_linear(ist)%zpsi)
      else
#ifdef R_TREAL
        call lalg_scal(np, aa_linear(ist), xx%states_linear(ist)%dpsi)
#endif
      end if
    end do
  end select

  call batch_pack_was_modified(xx)

  SAFE_DEALLOCATE_A(aa_linear)

  call profiling_out(scal_prof)
  POP_SUB(X(batch_scal_vec))
end subroutine X(batch_scal_vec)

! --------------------------------------------------------------

subroutine X(batch_xpay_vec)(np, xx, aa, yy, a_start, a_full)
  integer,           intent(in)    :: np
  type(batch_t),     intent(in)    :: xx
  R_TYPE,            intent(in)    :: aa(:)
  type(batch_t),     intent(inout) :: yy
  integer, optional, intent(in)    :: a_start
  logical, optional, intent(in)    :: a_full

  integer :: ist, ip, effsize, iaa
  R_TYPE, allocatable     :: aa_linear(:)
#ifdef HAVE_OPENCL
  integer :: size_factor, localsize
  FLOAT,  allocatable     :: aa_linear_double(:)
  type(opencl_mem_t)      :: aa_buffer
  type(octcl_kernel_t), save :: kernel
  type(cl_kernel)         :: kernel_ref
#endif
  
  PUSH_SUB(X(batch_xpay_vec))
  call profiling_in(xpay_prof, "BATCH_XPAY")

  ASSERT(batch_type(yy) == batch_type(xx))
#ifdef R_TCMPLX
  !if aa is complex, the functions must be complex
  ASSERT(batch_type(yy) == TYPE_CMPLX)
#endif
  ASSERT(xx%nst_linear == yy%nst_linear)
  ASSERT(batch_status(xx) == batch_status(yy))

  effsize = yy%nst_linear
  if(batch_is_packed(yy)) effsize = yy%pack%size(1)
  SAFE_ALLOCATE(aa_linear(1:effsize))

  aa_linear = M_ZERO
  do ist = 1, yy%nst_linear
    iaa = batch_linear_to_ist(xx, ist) - (optional_default(a_start, 1) - 1)
    if(.not. optional_default(a_full, .true.)) iaa = iaa - (batch_linear_to_ist(xx, 1) - 1)
    aa_linear(ist) = aa(iaa)
  end do

  select case(batch_status(xx))
  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL

    if(batch_type(yy) == TYPE_CMPLX .and. R_TYPE_VAL == TYPE_FLOAT) then
      size_factor = 2
      SAFE_ALLOCATE(aa_linear_double(1:2*yy%pack%size(1)))
      do ist = 1, yy%pack%size(1)
        aa_linear_double(2*ist - 1) = aa_linear(ist)
        aa_linear_double(2*ist) = aa_linear(ist)
      end do
      call opencl_create_buffer(aa_buffer, CL_MEM_READ_ONLY, TYPE_FLOAT, 2*yy%pack%size(1))
      call opencl_write_buffer(aa_buffer, 2*yy%pack%size(1), aa_linear_double)
      SAFE_DEALLOCATE_A(aa_linear_double)
    else
      size_factor = 1
      call opencl_create_buffer(aa_buffer, CL_MEM_READ_ONLY, R_TYPE_VAL, yy%pack%size(1))
      call opencl_write_buffer(aa_buffer, yy%pack%size(1), aa_linear)
    end if

    call octcl_kernel_start_call(kernel, 'axpy.cl', TOSTRING(X(xpay_vec)), flags = '-D'//R_TYPE_CL)
  
    kernel_ref = octcl_kernel_get_ref(kernel)

    call opencl_set_kernel_arg(kernel_ref, 0, aa_buffer)
    call opencl_set_kernel_arg(kernel_ref, 1, xx%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 2, log2(xx%pack%size(1)*size_factor))
    call opencl_set_kernel_arg(kernel_ref, 3, yy%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 4, log2(yy%pack%size(1)*size_factor))

    localsize = opencl_max_workgroup_size()
    call opencl_kernel_run(kernel_ref, (/yy%pack%size(1)*size_factor, pad(np, localsize)/), &
      (/yy%pack%size(1)*size_factor, localsize/(yy%pack%size(1)*size_factor)/))

    call opencl_finish()

    call opencl_release_buffer(aa_buffer)
#endif
  case(BATCH_PACKED)
    if(batch_type(yy) == TYPE_CMPLX) then
      do ist = 1, yy%pack%size(1)
        do ip = 1, np
          yy%pack%zpsi(ist, ip) = xx%pack%zpsi(ist, ip) + aa_linear(ist)*yy%pack%zpsi(ist, ip)
        end do
      end do
    else
#ifdef R_TREAL
      do ist = 1, yy%pack%size(1)
        do ip = 1, np
          yy%pack%dpsi(ist, ip) = xx%pack%dpsi(ist, ip) + aa_linear(ist)*yy%pack%dpsi(ist, ip)
        end do
      end do
#endif
    end if
    
  case(BATCH_NOT_PACKED)
    do ist = 1, yy%nst_linear
      if(batch_type(yy) == TYPE_CMPLX) then
        do ip = 1, np
          yy%states_linear(ist)%zpsi(ip) = xx%states_linear(ist)%zpsi(ip) + aa_linear(ist)*yy%states_linear(ist)%zpsi(ip)
        end do
      else
#ifdef R_TREAL
        do ip = 1, np
          yy%states_linear(ist)%dpsi(ip) = xx%states_linear(ist)%dpsi(ip) + aa_linear(ist)*yy%states_linear(ist)%dpsi(ip)
        end do
#endif
      end if
    end do
  end select

  call batch_pack_was_modified(yy)

  SAFE_DEALLOCATE_A(aa_linear)

  call profiling_out(xpay_prof)
  POP_SUB(X(batch_xpay_vec))
end subroutine X(batch_xpay_vec)

! --------------------------------------------------------------

subroutine X(batch_xpay_const)(np, xx, aa, yy)
  integer,           intent(in)    :: np
  type(batch_t),     intent(in)    :: xx
  R_TYPE,            intent(in)    :: aa
  type(batch_t),     intent(inout) :: yy

  integer :: minst, maxst, ii, ist
  R_TYPE, allocatable :: aavec(:)
  
  minst = HUGE(minst)
  maxst = -HUGE(maxst)
  
  do ii = 1, xx%nst_linear
    ist = batch_linear_to_ist(xx, ii)
    minst = min(minst, ist)
    maxst = max(maxst, ist)
  end do


  SAFE_ALLOCATE(aavec(minst:maxst))

  aavec = aa

  call X(batch_xpay_vec)(np, xx, aavec, yy, a_start = minst)

  SAFE_DEALLOCATE_A(aavec)
  
end subroutine X(batch_xpay_const)

! --------------------------------------------------------------

subroutine X(batch_set_state1)(this, ist, np, psi)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: ist
  integer,        intent(in)    :: np
  R_TYPE,         intent(in)    :: psi(:)

  integer :: ip
  type(profile_t), save :: prof
  type(opencl_mem_t) :: tmp

  call profiling_in(prof, "BATCH_SET_STATE")

  PUSH_SUB(X(batch_set_state1))

  ASSERT(ist >= 1 .and. ist <= this%nst_linear)
#ifdef R_TCOMPLEX
  ! cannot set a real batch with complex values
  ASSERT(batch_type(this) /= TYPE_FLOAT)
#endif

  call batch_pack_was_modified(this)

  select case(batch_status(this))
  case(BATCH_NOT_PACKED)
    if(batch_type(this) == TYPE_FLOAT) then
      forall(ip = 1:np) this%states_linear(ist)%dpsi(ip) = psi(ip)
    else
      forall(ip = 1:np) this%states_linear(ist)%zpsi(ip) = psi(ip)
    end if
  case(BATCH_PACKED)
    if(batch_type(this) == TYPE_FLOAT) then
      forall(ip = 1:np) this%pack%dpsi(ist, ip) = psi(ip)
    else
      forall(ip = 1:np) this%pack%zpsi(ist, ip) = psi(ip)
    end if
  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
    call opencl_create_buffer(tmp, CL_MEM_READ_ONLY, batch_type(this), this%pack%size(2))

    call opencl_write_buffer(tmp, np, psi)

    ! now call an opencl kernel to rearrange the data
    call opencl_set_kernel_arg(X(pack), 0, this%pack%size(1))
    call opencl_set_kernel_arg(X(pack), 1, ist - 1)
    call opencl_set_kernel_arg(X(pack), 2, tmp)
    call opencl_set_kernel_arg(X(pack), 3, this%pack%buffer)
    
    call opencl_kernel_run(X(pack), (/this%pack%size(2), 1/), (/opencl_max_workgroup_size(), 1/))
    
    call opencl_finish()

    call opencl_release_buffer(tmp)
#endif
  end select

  call profiling_out(prof)

  POP_SUB(X(batch_set_state1))
end subroutine X(batch_set_state1)

! --------------------------------------------------------------

subroutine X(batch_set_state2)(this, index, np, psi)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: index(:)
  integer,        intent(in)    :: np
  R_TYPE,         intent(in)    :: psi(:)

  PUSH_SUB(X(batch_set_state2))

  ASSERT(this%nst_linear > 0)
  call X(batch_set_state1)(this, batch_inv_index(this, index), np, psi)

  POP_SUB(X(batch_set_state2))
end subroutine X(batch_set_state2)

! --------------------------------------------------------------

subroutine X(batch_set_state3)(this, ii, np, psi)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: ii
  integer,        intent(in)    :: np
  R_TYPE,         intent(in)    :: psi(:, :)

  integer :: i2

  PUSH_SUB(X(batch_set_state3))

  do i2 = 1, this%dim
    call X(batch_set_state1)(this, (ii - 1)*this%dim + i2, np, psi(:, i2))
  end do

  POP_SUB(X(batch_set_state3))
end subroutine X(batch_set_state3)

! --------------------------------------------------------------

subroutine X(batch_get_state1)(this, ist, np, psi)
  type(batch_t),  intent(in)    :: this
  integer,        intent(in)    :: ist
  integer,        intent(in)    :: np
  R_TYPE,         intent(inout) :: psi(:)

  integer :: ip
  type(profile_t), save :: prof 
  type(opencl_mem_t) :: tmp

  PUSH_SUB(X(batch_get_state1))

  call profiling_in(prof, "BATCH_GET_STATE")

  ASSERT(ist >= 1 .and. ist <= this%nst_linear)
#ifdef R_TREAL
  ! cannot get a real value from a complex batch
  ASSERT(batch_type(this) /= TYPE_CMPLX)
#endif

  select case(batch_status(this))
  case(BATCH_NOT_PACKED)
    if(batch_type(this) == TYPE_FLOAT) then
      !$omp parallel do
      do ip = 1, np
        psi(ip) = this%states_linear(ist)%dpsi(ip)
      end do
      !$omp end parallel do
    else
      !$omp parallel do
      do ip = 1, np
        psi(ip) = this%states_linear(ist)%zpsi(ip)
      end do
      !$omp end parallel do
    end if

  case(BATCH_PACKED)
    if(batch_type(this) == TYPE_FLOAT) then
      !$omp parallel do
      do ip = 1, np
        psi(ip) = this%pack%dpsi(ist, ip)
      end do
      !$omp end parallel do
    else
      !$omp parallel do
      do ip = 1, np
        psi(ip) = this%pack%zpsi(ist, ip)
      end do
      !$omp end parallel do
    end if

  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
    call opencl_create_buffer(tmp, CL_MEM_WRITE_ONLY, batch_type(this), this%pack%size(2))

    call opencl_set_kernel_arg(X(unpack), 0, this%pack%size(1))
    call opencl_set_kernel_arg(X(unpack), 1, ist - 1)
    call opencl_set_kernel_arg(X(unpack), 2, this%pack%buffer)
    call opencl_set_kernel_arg(X(unpack), 3, tmp)

    call opencl_kernel_run(X(unpack), (/1, this%pack%size(2)/), (/1, opencl_max_workgroup_size()/))

    call opencl_finish()

    call opencl_read_buffer(tmp, np, psi)

    call opencl_release_buffer(tmp)
#endif
  end select

  call profiling_out(prof)

  POP_SUB(X(batch_get_state1))
end subroutine X(batch_get_state1)

! --------------------------------------------------------------

subroutine X(batch_get_state2)(this, index, np, psi)
  type(batch_t),  intent(in)    :: this
  integer,        intent(in)    :: index(:)
  integer,        intent(in)    :: np
  R_TYPE,         intent(inout) :: psi(:)

  PUSH_SUB(X(batch_get_state2))

  ASSERT(this%nst_linear > 0)
  call X(batch_get_state1)(this, batch_inv_index(this, index), np, psi)

  POP_SUB(X(batch_get_state2))
end subroutine X(batch_get_state2)


! --------------------------------------------------------------

subroutine X(batch_get_state3)(this, ii, np, psi)
  type(batch_t),  intent(in)    :: this
  integer,        intent(in)    :: ii
  integer,        intent(in)    :: np
  R_TYPE,         intent(inout) :: psi(:, :)

  integer :: i2

  PUSH_SUB(X(batch_get_state3))

  do i2 = 1, this%dim
    call X(batch_get_state1)(this, (ii - 1)*this%dim + i2, np, psi(:, i2))
  end do

  POP_SUB(X(batch_get_state3))
end subroutine X(batch_get_state3)

! --------------------------------------------------------------

subroutine X(batch_get_points)(this, sp, ep, psi)
  type(batch_t),  intent(in)    :: this
  integer,        intent(in)    :: sp  
  integer,        intent(in)    :: ep
  R_TYPE,         intent(inout) :: psi(:, :, sp:)

  integer :: idim, ist, ii, ip
  
  PUSH_SUB(X(batch_get_points))
  call profiling_in(get_points_prof, 'GET_POINTS')

#ifdef R_TREAL
  ! cannot get a real value from a complex batch
  ASSERT(batch_type(this) /= TYPE_CMPLX)
#endif

  select case(batch_status(this))
  case(BATCH_NOT_PACKED)

    if(batch_type(this) == TYPE_FLOAT) then
      
      do ii = 1, this%nst_linear
        ist = batch_linear_to_ist(this, ii)
        idim = batch_linear_to_idim(this, ii)
        psi(ist, idim, sp:ep) = this%states_linear(ii)%dpsi(sp:ep)
      end do

    else

      do ii = 1, this%nst_linear
        ist = batch_linear_to_ist(this, ii)
        idim = batch_linear_to_idim(this, ii)
        psi(ist, idim, sp:ep) = this%states_linear(ii)%zpsi(sp:ep)
      end do

    end if

  case(BATCH_PACKED)

    if(batch_type(this) == TYPE_FLOAT) then

      !$omp parallel do private(ip, ii, ist, idim)
      do ip = sp, ep
        do ii = 1, this%nst_linear
          ist = batch_linear_to_ist(this, ii)
          idim = batch_linear_to_idim(this, ii)
          psi(ist, idim, ip) = this%pack%dpsi(ii, ip)
        end do
      end do
      !$omp end parallel do

    else

      !$omp parallel do private(ip, ii, ist, idim)
      do ip = sp, ep
        do ii = 1, this%nst_linear
          ist = batch_linear_to_ist(this, ii)
          idim = batch_linear_to_idim(this, ii)
          psi(ist, idim, ip) = this%pack%zpsi(ii, ip)
        end do
      end do
      !$omp end parallel do

    end if

  case(BATCH_CL_PACKED)
    call messages_not_implemented('batch_get_points for CL packed batches')
  end select

  call profiling_out(get_points_prof)

  POP_SUB(X(batch_get_points))
end subroutine X(batch_get_points)

! --------------------------------------------------------------

subroutine X(batch_set_points)(this, sp, ep, psi)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: sp  
  integer,        intent(in)    :: ep
  R_TYPE,         intent(in)    :: psi(:, :, sp:)

  integer :: idim, ist, ii, ip

  PUSH_SUB(X(batch_set_points))

  call profiling_in(set_points_prof, 'SET_POINTS')

#ifdef R_TCOMPLEX
  ! cannot set a real batch with complex values
  ASSERT(batch_type(this) /= TYPE_FLOAT)
#endif

  call batch_pack_was_modified(this)

  select case(batch_status(this))
  case(BATCH_NOT_PACKED)

    if(batch_type(this) == TYPE_FLOAT) then

      do ii = 1, this%nst_linear
        ist = batch_linear_to_ist(this, ii)
        idim = batch_linear_to_idim(this, ii)
        this%states_linear(ii)%dpsi(sp:ep) = psi(ist, idim, sp:ep)
      end do

    else

      do ii = 1, this%nst_linear
        ist = batch_linear_to_ist(this, ii)
        idim = batch_linear_to_idim(this, ii)
        this%states_linear(ii)%zpsi(sp:ep) = psi(ist, idim, sp:ep)
      end do

    end if

  case(BATCH_PACKED)

    if(batch_type(this) == TYPE_FLOAT) then

      !$omp parallel do private(ip, ii, ist, idim)
      do ip = sp, ep
        do ii = 1, this%nst_linear
          ist = batch_linear_to_ist(this, ii)
          idim = batch_linear_to_idim(this, ii)
          this%pack%dpsi(ii, ip) = psi(ist, idim, ip)
        end do
      end do
      !$omp end parallel do

    else

      !$omp parallel do private(ip, ii, ist, idim)
      do ip = sp, ep
        do ii = 1, this%nst_linear
          ist = batch_linear_to_ist(this, ii)
          idim = batch_linear_to_idim(this, ii)
          this%pack%zpsi(ii, ip) = psi(ist, idim, ip)
        end do
      end do
      !$omp end parallel do

    end if

  case(BATCH_CL_PACKED)
    call messages_not_implemented('batch_set_points for CL packed batches')
  end select

  call profiling_out(set_points_prof)

  POP_SUB(X(batch_set_points))
end subroutine X(batch_set_points)

! --------------------------------------------------------------

subroutine X(batch_mul)(np, ff,  xx, yy)
  integer,           intent(in)    :: np
  R_TYPE,            intent(in)    :: ff(:)
  type(batch_t),     intent(in)    :: xx
  type(batch_t),     intent(inout) :: yy

  integer :: ist, ip
  R_TYPE :: mul

  PUSH_SUB(X(batch_mul))
  call profiling_in(mul_prof, "BATCH_MUL")

  ASSERT(batch_type(yy) == batch_type(xx))
#ifdef R_TCMPLX
  !if aa is complex, the functions must be complex
  ASSERT(batch_type(yy) == TYPE_CMPLX)
#endif
  ASSERT(xx%nst_linear == yy%nst_linear)
  ASSERT(batch_status(xx) == batch_status(yy))

  select case(batch_status(yy))
  case(BATCH_CL_PACKED)
    call messages_not_implemented("OpenCL batch_mul")

  case(BATCH_PACKED)
    if(batch_type(yy) == TYPE_CMPLX) then
      !$omp parallel do private(ip, ist, mul)
      do ip = 1, np
        mul = ff(ip)
        do ist = 1, yy%nst_linear
          yy%pack%zpsi(ist, ip) = mul*xx%pack%zpsi(ist, ip)
        end do
      end do
      !$omp end parallel do
    else
#ifdef R_TREAL
      !$omp parallel do private(ip, ist, mul)
      do ip = 1, np
        mul = ff(ip)
        do ist = 1, yy%nst_linear
          yy%pack%dpsi(ist, ip) = mul*xx%pack%dpsi(ist, ip)
        end do
      end do
      !$omp end parallel do
#endif
    end if

  case(BATCH_NOT_PACKED)
    if(batch_type(yy) == TYPE_CMPLX) then
      do ist = 1, yy%nst_linear
        !$omp parallel do
        do ip = 1, np
          yy%states_linear(ist)%zpsi(ip) = ff(ip)*xx%states_linear(ist)%zpsi(ip)
        end do
        !$omp end parallel do
      end do
    else
#ifdef R_TREAL
      do ist = 1, yy%nst_linear
        !$omp parallel do
        do ip = 1, np
          yy%states_linear(ist)%zpsi(ip) = ff(ip)*xx%states_linear(ist)%zpsi(ip)
        end do
        !$omp end parallel do
      end do

#endif
    end if
  end select

  call batch_pack_was_modified(yy)

  call profiling_out(mul_prof)
  POP_SUB(X(batch_mul))

end subroutine X(batch_mul)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
