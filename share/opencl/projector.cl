/*
 Copyright (C) 2010 X. Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: projector.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void projector_bra(const int imat,
			    const int npoints,
			    const int nprojs,
			    const __global int * offsets,
			    __global const double * matrix,
			    __global const int * map,
			    __global const double * scal,
			    __global const double * psi, const int ldpsi,
			    __global double * projection, const int ldprojection, const int projection_offset
			    ){
  
  const int ist = get_global_id(0);
  const int ipj = get_global_id(1);

  if(ipj >= nprojs) return;
  
  const int matrix_offset = offsets[4*imat    ];
  const int map_offset    = offsets[4*imat + 1];
  const int scal_offset   = offsets[4*imat + 2];

  double aa = 0.0;
  for(int ip = 0; ip < npoints; ip++){
    aa += matrix[matrix_offset + ip + npoints*ipj]*psi[ldpsi*(map[map_offset + ip] - 1) + ist];
  }
  projection[projection_offset + ist + ldprojection*ipj] = scal[scal_offset + ipj]*aa;

}

__kernel void projector_ket(const int imat, 
			    const int npoints,
			    const int nprojs,
			    const __global int * offsets,
			    __global const double * matrix,
			    __global const int * map,
			    __global const double * projection, const int ldprojection, const int projection_offset,
			    __global double * psi, const int ldpsi
			    ){
  
  int ist = get_global_id(0);
  int ip = get_global_id(1);

  if(ip >= npoints) return;

  const int matrix_offset = offsets[4*imat    ];
  const int map_offset    = offsets[4*imat + 1];

  double aa = 0.0;
  for(int ipj = 0; ipj < nprojs; ipj++){
    aa += matrix[matrix_offset + ip + npoints*ipj]*projection[projection_offset + ist + ldprojection*ipj];
  }
  psi[ldpsi*(map[map_offset + ip] - 1) + ist] += aa;


}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/

