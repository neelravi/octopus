/*
 Copyright (C) 2012 X. Andrade

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

 $Id: pack.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#include <cl_global.h>
#include <cl_complex.h>
#include <cl_rtype.h>

__kernel void X(zmul)(const int np,
		      __global rtype const * restrict op,
		      __global double2 * restrict ff){
  
  const int ip = get_global_id(0);
  if(ip < np) ff[ip] = MUL(op[ip], ff[ip]);
}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
