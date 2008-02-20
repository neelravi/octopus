/*
 Copyright (C) 2006-2007 M.A.L. Marques

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

#define XC_GGA_X_PW91         109 /* Perdew & Wang 91 */
#define XC_GGA_X_mPW91        119 /* Modified form of PW91 by Adamo & Barone */

static inline void 
func(xc_gga_type *p, double x, double *f, double *dfdx, double *ldfdx)
{
  const double x2s     = 0.12827824385304220645; /* 1/(2*(6*pi^2)^(1/3)) */
  const double aa[2]   = {0.19645, 0.21516};
  const double bb      = 7.7956;
  const double cc[2]   = {0.2743, 0.30042};
  const double dd[2]   = {-0.1508, -0.17696};
  const double ff[2]   = {0.004, 0.002285};
  const double alpha   = 100.0;
  const double expo[2] = {4.0, 3.75};

  double ss, ss2, ss4;
  double f1, f2, f3, f4;

  int func;

  switch(p->info->number){
  case XC_GGA_X_mPW91:   func = 1; break;
  default:               func = 0; /* original PW91 */
  }

  ss  = x2s*x;
  ss2 = ss*ss;
  ss4 = pow(ss, expo[func]);

  f1 = dd[func]*exp(-alpha*ss2);
  f2 = aa[func]*asinh(bb*ss);
  f3 = (cc[func] + f1)*ss2 - ff[func]*ss4;
  f4 = 1.0 + ss*f2 + ff[func]*ss4;

  *f     = 1.0 + f3/f4;
  {
    double df3, df4;

    df3 = 2.0*ss*(cc[func] + f1*(1.0 - alpha*ss2) - expo[func]/2.0*ff[func]*pow(ss, expo[func]-2.0));
    df4 = f2 + ss*(aa[func]*bb/sqrt(1.0 + bb*bb*ss2) + expo[func]*ff[func]*pow(ss, expo[func]-2.0));
    *dfdx  = (df3*f4 - f3*df4)/(f4*f4);
    *ldfdx = cc[func] + dd[func];
  }


  *dfdx  *= x2s;
  *ldfdx *= x2s*x2s;
}

#include "work_gga_x.c"

const xc_func_info_type func_info_gga_x_pw91 = {
  XC_GGA_X_PW91,
  XC_EXCHANGE,
  "Perdew & Wang 91",
  XC_FAMILY_GGA,
  "JP Perdew, in Proceedings of the 21st Annual International Symposium on the Electronic Structure of Solids, ed. by P Ziesche and H Eschrig (Akademie Verlag, Berlin, 1991), p. 11.\n"
  "JP Perdew, JA Chevary, SH Vosko, KA Jackson, MR Pederson, DJ Singh, and C Fiolhais, Phys. Rev. B 46, 6671 (1992)\n"
  "JP Perdew, JA Chevary, SH Vosko, KA Jackson, MR Pederson, DJ Singh, and C Fiolhais, Phys. Rev. B 48, 4978(E) (1993)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};


const xc_func_info_type func_info_gga_x_mpw91 = {
  XC_GGA_X_mPW91,
  XC_EXCHANGE,
  "mPW91 of Adamo & Barone",
  XC_FAMILY_GGA,
  "C Adamo and V Barone, J. Chem. Phys. 108, 664 (1998)",
  XC_PROVIDES_EXC | XC_PROVIDES_VXC,
  NULL, NULL, NULL,
  work_gga_x
};
