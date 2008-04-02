/*
 Copyright (C) 2008 Georg Madsen

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
  
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <assert.h>
#include "util.h"

#define XC_GGA_X_PBEA  121 /* Madsen (PBE-like) */

/* PBEA: see PBE for more details */
static inline void 
func(const XC(gga_type) *p, FLOAT x, FLOAT *f, FLOAT *dfdx, FLOAT *ldfdx, FLOAT *d2fdx2)
{
  static const FLOAT kappa = 0.8040;
  static const FLOAT mu = 0.00361218645365094697;
  /* hard-coded alpha*/
  static const FLOAT alpha = 0.5;

  FLOAT dd;

  dd     = 1.0 + mu*x*x/(alpha*kappa);

  *f     = 1.0 + kappa*(1.0 - 1.0/POW(dd, alpha));
  *dfdx  = 2.0*mu*x/POW(dd, alpha + 1.0);
  *ldfdx = mu;
}

#include "work_gga_x.c"

const XC(func_info_type) XC(func_info_gga_x_pbea) = {
  XC_GGA_X_PBEA,
  XC_EXCHANGE,
  "Madsen 07",
  XC_FAMILY_GGA,
  "G Madsen, Phys. Rev. B 75, 195108 (2007)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};
