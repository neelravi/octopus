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

module metric_m
  use datasets_m
  use geometry_m
  use global_m
  use loct_m
  use math_m
  use messages_m
  use parser_m
  use profiling_m
  use unit_m
  use unit_system_m
  use utils_m

  implicit none

  private
    integer   ::  i, j, info
    real*8,   ::  WORK(3)               ! temporary array for matrix inversion
    real*8,   ::  IPIV(3)               ! temporary array for matrix inversion

  public ::                       &
    metric_t,                     &
    metric_init,                  &
    metric_end,                   &
    metric_copy,                  &
    metric_write_info

  type metric_t
    FLOAT :: tensor(1:3, 1:3)
    FLOAT :: br_vecs(1:3)
    FLOAT :: br_angles(1:3)
  end type metric_t

contains

  subroutine metric_init(this)
    type(metric_t),    intent(out) :: this

    PUSH_SUB(metric_init)

    this%tensor = M_ZERO
!    forall(i=1:3) this%tensor(i:i) = 1.0d0

    this%br_vecs(1:3)   = M_ONE
    this%br_angles(1:3) = M_PI*M_HALF

  !%Variable BravaisLattice
  !%Type integer
  !%Default cubic_primitive
  !%Section Mesh Simulation Box
  !%Description
  !% Keyword for calculations on periodic systems. Sets up Bravais Lattice index, and builds
  !% the lattice vectors accordingly.
  !%Option cubic_primitive 01
  !% Simple Cubic system.                 v1 = a(1,0,0),       v2 = a(0,1,0),      v3 = a(0,0,1)
  !%Option cubic_face_centered 02
  !% Face Centered Cubic system.          v1 = (a/2)(-1,0,1),  v2 = (a/2)(0,1,1),  v3 = (a/2)(-1,1,0)
  !%Option cubic_body_centered 03
  !% Body Centered Cubic system           v1 = (a/2)(1,1,1),   v2 = (a/2)(-1,1,1), v3 = (a/2)(-1,-1,1)
  !%Option Triagonal 04
  !% Triagonal (primitive cell is a simple rhombohedron) 
  !% v1 = a(tx,-ty,tz),   v2 = a(0,2ty,tz),   v3 = a(-tx,-ty,tz)
  !%    where c=cos(alpha) is the cosine of the angle alpha between
  !%    any pair of crystallographic vectors, tx, ty, tz are:
  !%      tx=sqrt((1-c)/2), ty=sqrt((1-c)/6), tz=sqrt((1+2c)/3) 
  !%Option hexagonal 05
  !% Hexagonal Lattice, Lattice constants a = b, c, v1 = a(1,0,0),  v2 = a(-1/2,sqrt(3)/2,0),  v3 = a(0,0,c/a)
  !%Option tetragonal_primitive 06
  !% Simple tetragonal system.            v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,c/a)
  !%Option tetragonal_body_centered 07
  !% Body centered tetragonal system.     v1=(a/2)(1,-1,c/a),  v2=(a/2)(1,1,c/a),  v3=(a/2)(-1,-1,c/a)
  !%Option orthorhombic_primitive  08
  !% Simple (P) orthorhomic system;   Lat. const. a != b != c,   v1 = (a,0,0),  v2 = (0,b,0), v3 = (0,0,c)
  !%Option orthorhombic_base_centered 09
  !% Base centered Orthorhombic;          v1 = (a/2, b/2,0),  v2 = (-a/2,b/2,0),  v3 = (0,0,c)
  !%Option orthorhombic_face_centered 10
  !% Face centered orthorhombic system.     v1 = (a/2,0,c/2),  v2 = (a/2,b/2,0),  v3 = (0,b/2,c/2)
  !%Option orthorhombic_body_centered 11
  !% Body centered orthorhombic system.     v1=(a/2,b/2,c/2),  v2=(-a/2,b/2,c/2),  v3=(-a/2,-b/2,c/2)
  !%Option monoclinic_primitive  12
  !%  monoclinic simple (P)                v1=(a,0,0), v2=(b*cos(gamma),b*sin(gamma),0),  v3 = (0,0,c)
  !%Option monoclinic_base_centered 13
  !% Base centered monoclinic;             v1 = (  a/2,         0,                -c/2),
  !%                                       v2 = (b*cos(gamma),  b*sin(gamma),        0),
  !%                                       v3 = (  a/2,         0,                 c/2),
  !%Option triclinic 14
  !% triclinic; the most general lattice;  v1 = (a, 0, 0),
  !%                                       v2 = (b*cos(gamma), b*sin(gamma), 0)
  !%                                       v3 = (c*cos(beta),  c*(cos(alpha)-cos(beta)cos(gamma))/sin(gamma),
  !%                                            c*sqrt( 1 + 2*cos(alpha)cos(beta)cos(gamma)
  !%                                            - cos(alpha)^2-cos(beta)^2-cos(gamma)^2 )/sin(gamma) )
  !% Example:
  !% <pre>%BravaisLattice
  !%   3.0 | 3.0 | 3.0
  !%    90 |  90 |  90
  !% %
  !% </pre>
  !%End

  ! Following tensor will be used for various transformations of lattice vectors

  this%tensor(1,1) =  this%br_vecs(1)*this%br_vecs(1)
  this%tensor(1,2) =  this%br_vecs(1)*this%br_vecs(2)*cos(this%br_angles(3))
  this%tensor(1,3) =  this%br_vecs(1)*this%br_vecs(3)*cos(this%br_angles(2))
  this%tensor(2,1) =  this%br_vecs(1)*this%br_vecs(2)*cos(this%br_angles(3))
  this%tensor(2,2) =  this%br_vecs(2)*this%br_vecs(2)
  this%tensor(2,3) =  this%br_vecs(2)*this%br_vecs(3)*cos(this%br_angles(1))
  this%tensor(3,1) =  this%br_vecs(1)*this%br_vecs(3)*cos(this%br_angles(2))
  this%tensor(3,2) =  this%br_vecs(2)*this%br_vecs(3)*cos(this%br_angles(1))
  this%tensor(3,3) =  this%br_vecs(3)*this%br_vecs(3)

   ! check inserted by ravindra, print the metric tensor and the lattice vectors

   write(*,*) "printing metric tensor "

   Do i = 1, 3
      write(*,*) (this%tensor(i,j),j=1,3)
   enddo
   

  ! Following part will describe the transformation matrix A that transforms the
  ! coordinates from simple Cartesian to the lattice coordinates. We will enlist
  ! the matrices for all types of bravais lattices. 
 
  ! NOTE: following matrices will be called depending upon the geometry read in
  ! the input file. So I am naming the matrix A the same in every sub part that
  ! follows. (this will be moved to appropriate subroutine at the end if
  ! desired) Matrix B is the the inverse of matrix A; F = B*B(tran)

  ! Case 1: Simple Cubic. There is no tranformation needed. Hence the
  ! transformation matrix is the unit matrix.

   A = 0.d0 ; Do i = 1 , 3 ;  A(i,i)  = 1.d0 ;   enddo ; B = A ; F = A

  ! Case 2: Face Centered Cubic.

   A = 0.d0 ; forall (i=1:3, j=1:3, i .ne. j) A(i,j) = 1.d0/sqrt(2.d0)

   B = -1.d0/sqrt(2.d0) ; forall (i=1:3, j=1:3, i .ne. j) B(i,j) = 1.d0/sqrt(2.d0)

   F = 3.d0/2.d0 ; forall (i=1:3, j=1:3, i .ne. j) F(i,j) = -1.d0/2.d0


  
  ! Case 3: Body Centered Cubic

   A = 1.d0 ; forall (i=1:3, j=1:3, i .eq. j) A(i,j) = -1.d0

   F = 3.d0/2.d0 ; forall (i=1:3, j=1:3, i .ne. j) F(i,j) = -1.d0/2.d0
  

  ! Case 5: Hexagonal Lattice
   A = 0.d0 ; A(1,1) = 1.d0; A(1,2) = 0.5d0 ; A(2,2) = sqrt(3.d0)/0.5d0 ; A(3,3) = 1.d0
   
   F = 0.d0 ; F(1,1) = 4.d0/3.d0 ; F(2,2) =  4.d0/3.d0 ; F(1,2) = -2.d0/3.d0 ; F(2,1) =  -2.d0/3.d0

  ! Case 6: Tetragonal_primitive

   A = 0.d0 ; Do i = 1 , 3 ;  A(i,i)  = 1.d0 ;   enddo ; B = A ; F = A


  ! Case 7: Tetragonal Body Centered

   factor = a*0.5d0/(sqrt(0.5*a*a+0.25*c*c))
   
   A(1,1) = 1.d0;       A(1,2) = 1.d0;    A(1,3) = -1.d0
   A(2,1) = -1.d0;      A(2,2) = 1.d0;    A(2,3) = -1.d0
   A(3,1) = c/a;        A(3,2) = c/a;     A(3,3) = c/a

   A = factor*A

     call DGETRF(3,3,A,3,IPIV,info)
     call DGETRI(3,A,3,IPIV,WORK,3,info)
  
   F = matmul(A,transpose(A)) 

  ! Case 8: Simple Orthorhombic primitiv (P) 
   
    A = 0.d0 ; Do i = 1 , 3 ;  A(i,i)  = 1.d0 ;   enddo ; B = A ; F = A

  ! Case 9: Base Centered Orthorhombic 

   factor = 1.d0/(sqrt(a*a+b*b))
   
   A(1,1) = a;        A(1,2) = -1.d0*a;    A(1,3) = 0.d0
   A(2,1) = b;        A(2,2) = b;          A(2,3) = 0.d0
   A(3,1) = 0.d0;        A(3,2) = 0.d0;    A(3,3) = sqrt(a*a+b*b)

   A = factor*A

     call DGETRF(3,3,A,3,IPIV,info)
     call DGETRI(3,A,3,IPIV,WORK,3,info)

   F = matmul(A,transpose(A))

  ! Case 10: Face centered orthorhombic 
   
   A(1,1) = a/(sqrt(a*a+c*c)) ; A(1,2) = a/(sqrt(a*a+c*c)) ; A(1,3) = 0.d0
   A(2,1) = 0.d0              ; A(2,2) = b/(sqrt(a*a+b*b)) ; A(2,3) = b/(sqrt(b*b+c*c))
   A(3,1) = c/(sqrt(a*a+c*c)) ; A(3,2) = 0.d0              ; A(3,3) = c/(sqrt(b*b+c*c))

     call DGETRF(3,3,A,3,IPIV,info)
     call DGETRI(3,A,3,IPIV,WORK,3,info)

   F = matmul(A,transpose(A))

  ! Case 11: Body Centered orthorhombic

  factor = 1.d0/(sqrt(a*a+b*b+c*c))

   A(1,1) = a ; A(1,2) = -a      ; A(1,3) = -a
   A(2,1) = b ; A(2,2) = b       ; A(2,3) = -b
   A(3,1) = c ; A(3,2) = c       ; A(3,3) = c

     call DGETRF(3,3,A,3,IPIV,info)
     call DGETRI(3,A,3,IPIV,WORK,3,info)

   F = matmul(A,transpose(A))
 
  ! case 12: Monoclinic primitiv
   
   A(1,1) = 1.d0 ; A(1,2) = cos(gamma)  ; A(1,3) = 0.d0
   A(2,1) = 0.d0 ; A(2,2) = sin(gamma)  ; A(2,3) = 0.d0
   A(3,1) = 0.d0 ; A(3,2) = 0.d0        ; A(3,3) = 1.d0

     call DGETRF(3,3,A,3,IPIV,info)
     call DGETRI(3,A,3,IPIV,WORK,3,info)

   F = matmul(A,transpose(A))

 ! Case 13: Base centered monoclinic
   A(1,1) = a                 ; A(1,2) = 0.d0  ; A(1,3) = -c/(sqrt(2.d0))
   A(2,1) = 2.d0*b*cos(gamma) ; A(2,2) = 1.d0  ; A(2,3) = 0.d0
   A(3,1) = a                 ; A(3,2) = 0.d0  ; A(3,3) = c/(sqrt(2.d0))

     call DGETRF(3,3,A,3,IPIV,info)
     call DGETRI(3,A,3,IPIV,WORK,3,info)

   F = matmul(A,transpose(A))
 

 ! Case 14: Triclinic
   A(1,1) = 1.d0 ; A(1,2) = cos(gamma)  ; A(1,3) = cos(beta)
   A(2,1) = 0.d0 ; A(2,2) = sin(gamma)  ; A(2,3) = (cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
   A(3,1) = 0.d0 ; A(3,2) = 0.d0        ; A(3,3) = sqrt(1.d0+2*cos(alpha)*cos(beta)*cos(gamma)  & 
                     -cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma))/sin(gamma)

     call DGETRF(3,3,A,3,IPIV,info)
     call DGETRI(3,A,3,IPIV,WORK,3,info)

   F = matmul(A,transpose(A))


!  Cases of all Bravais Lattices end here

    POP_SUB(metric_init)
  end subroutine metric_init

  ! ---------------------------------------------------------
  subroutine metric_end(this)
    type(metric_t), intent(inout) :: this

    PUSH_SUB(metric_end)


    POP_SUB(metric_end)
  end subroutine metric_end


  ! ---------------------------------------------------------
  subroutine metric_copy(bb, aa)
    type(metric_t), intent(in)  :: bb
    type(metric_t), intent(out) :: aa

    PUSH_SUB(metric_copy)


    POP_SUB(metric_copy)
  end subroutine metric_copy





  ! ---------------------------------------------------------
  subroutine metric_write_info(this, iunit)
    type(metric_t),    intent(in) :: this
    integer,            intent(in) :: iunit

    integer :: ik, idir
    character(len=100) :: str_tmp
    character :: index

    PUSH_SUB(metric_write_info)

    POP_SUB(metric_write_info)
  end subroutine metric_write_info

end module metric_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
