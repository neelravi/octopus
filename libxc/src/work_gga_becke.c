/*
 Copyright (C) 2006-2008 M.A.L. Marques

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

/************************************************************************
  This file is to be included in meta GGA exchange functionals. As often these
  functionals are written as a function of s = |grad n|/n^(4/3) and tau, this
  routine performs the necessary conversions between a functional of s and tau
  and of rho.
************************************************************************/

static void
work_gga_becke_init(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  p->lda_aux = (XC(lda_type) *) malloc(sizeof(XC(lda_type)));
  XC(lda_init)(p->lda_aux, XC_LDA_C_PW, XC_POLARIZED);
}


static void
work_gga_becke_end(void *p_)
{
  XC(gga_type) *p = (XC(gga_type) *)p_;

  /* XC(lda_end)(p->lda_aux); */
  free(p->lda_aux);
}


static void 
work_gga_becke(const void *p_, const FLOAT *rho, const FLOAT *sigma,
	       FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
	       FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2)
{
  const XC(gga_type) *p = p_;

  FLOAT sfact, sfact2, dens;
  FLOAT ds[2], sigmas[2], x[2], x_avg;
  FLOAT e_LDA_opp, v_LDA_opp[2], f_LDA_opp[3];
  int   is, order;

  order = 0;
  if(vrho   != NULL) order = 1;
  if(v2rho2 != NULL) order = 2;

  sfact = (p->nspin == XC_POLARIZED) ? 1.0 : 2.0;
  sfact2 = sfact*sfact;

  if(p->nspin == XC_POLARIZED){
    ds[0] = rho[0];     ds[1] = rho[1];
    dens  = rho[0] + rho[1];
  }else{
    ds[0] = rho[0]/2.0; ds[1] = ds[0];
    dens  = rho[0];
  }

  /* get spin-polarized LDA */
  switch (order){
  case 0:
    XC(lda_exc)(p->lda_aux, ds, &e_LDA_opp);
    break;
  case 1:
    XC(lda_vxc)(p->lda_aux, ds, &e_LDA_opp, v_LDA_opp);
    break;
  case 2: /* to be implemented */
    XC(lda)(p->lda_aux, ds, &e_LDA_opp, v_LDA_opp, f_LDA_opp, NULL);
    break; 
  }
  e_LDA_opp *= dens;

  if(p->nspin == XC_UNPOLARIZED)
    f_LDA_opp[0] = (f_LDA_opp[0] + 2.0*f_LDA_opp[1] + f_LDA_opp[2])/4.0;

  x_avg = 0.0;
  for(is=0; is<p->nspin; is++){
    FLOAT gdm, rho13, mrho[2];
    FLOAT e_x, e_ss, e_LDA, v_LDA[2], f_LDA[3];
    FLOAT g_x, dg_x, d2g_x, g_ss, dg_ss, d2g_ss;
    int js = (is == 0) ? 0 : 2;

    if(rho[is] < MIN_DENS) continue;

    sigmas[is] = max(MIN_GRAD*MIN_GRAD, sigma[js]/sfact2);
    gdm    = sqrt(sigmas[is]);
  
    rho13 = POW(ds[is], 1.0/3.0);
    x[is] = gdm/(ds[is]*rho13);
    x_avg+= sfact*0.5*x[is]*x[is];

    func_gga_becke_exchange(p, x[is], order, &g_x, &dg_x, &d2g_x);
    e_x  = -sfact*X_FACTOR_C*(ds[is]*rho13);

    func_gga_becke_parallel(p, x[is], order, &g_ss, &dg_ss, &d2g_ss);

    /* get parallel spin LDA energy */
    mrho[0] = ds[is];
    mrho[1] = 0.0;
    switch (order){
    case 0:
      XC(lda_exc)(p->lda_aux, mrho, &e_LDA);
      break;
    case 1:
      XC(lda_vxc)(p->lda_aux, mrho, &e_LDA, v_LDA);
      break;
    case 2:
      XC(lda)(p->lda_aux, mrho, &e_LDA, v_LDA, f_LDA, NULL);
      break;
    }

    e_ss       = sfact*ds[is]*e_LDA;
    e_LDA_opp -= e_ss;

    if(zk != NULL)
      *zk += e_x*g_x + e_ss*g_ss;
 
    if(vrho != NULL){
      vrho[is]      += -4.0/3.0*X_FACTOR_C*rho13*(g_x - dg_x*x[is]);
      vrho[is]      += v_LDA[0]*g_ss - 4.0/3.0*e_LDA*dg_ss*x[is];
      
      vsigma[js]     = (e_x*dg_x + e_ss*dg_ss) * x[is]/(2.0*sfact2*sigmas[is]);

      v_LDA_opp[is] -= v_LDA[0];
    }
    
    if(v2rho2 != NULL){
      int ks2 = (is == 0) ? 0 : 5;

      v2rho2[js] += -4.0/9.0*X_FACTOR_C/(rho13*rho13) *
	(g_x - dg_x*x[is] + 4.0*d2g_x*x[is]*x[is]);

      v2rho2[js] += f_LDA[0]*g_ss - 4.0*x[is]/(3.0*ds[is])*
	(2.0*v_LDA[0]*dg_ss - 7.0*e_LDA*dg_ss/3.0 - 4.0*e_LDA*x[is]*d2g_ss/3.0);

      v2rho2[js] /= sfact;

      v2rhosigma[ks2] += x[is]/(2.0*sfact2*sigmas[is])*
	(4.0/3.0*X_FACTOR_C*rho13 * d2g_x*x[is] +
	 v_LDA[0]*dg_ss - 4.0/3.0*e_LDA*(d2g_ss*x[is] + dg_ss));

      v2sigma2[ks2] = (e_x*(d2g_x*x[is] - dg_x) + e_ss*(d2g_ss*x[is] - dg_ss)) * 
	x[is]/(4.0*sfact2*sfact2*sigmas[is]*sigmas[is]);

      f_LDA_opp[js] -= f_LDA[0]/sfact;
    }
  }

  /* We are now missing the opposite-spin part */
  {
    FLOAT g_ab, dg_ab, d2g_ab, dx_avg[2];
 
    x_avg = sqrt(x_avg);
    func_gga_becke_opposite(p, x_avg, order, &g_ab, &dg_ab, &d2g_ab);

    if(zk != NULL)
      *zk += e_LDA_opp*g_ab;
 
    if(vrho != NULL){
      for(is=0; is<p->nspin; is++){
	FLOAT dd;
	int js = (is == 0) ? 0 : 2;

	vrho[is] += v_LDA_opp[is]*g_ab;

	dd = POW(dens, 4.0/3.0);
	if(x_avg*dd < MIN_GRAD*MIN_GRAD || ds[is] < MIN_DENS) continue;

	dx_avg[is]  = -2.0*x[is]*x[is]/(3.0*x_avg*ds[is]);

	vrho[is]   += e_LDA_opp*dg_ab*dx_avg[is];
	vsigma[js] += e_LDA_opp*dg_ab*POW(ds[is], -8.0/3.0)/(sfact*4.0*x_avg);
      }
    }

    if(v2rho2 != NULL){
      int ks, ncomp;

      ncomp = (p->nspin == XC_POLARIZED) ? 3 : 1;
      for(ks=0; ks<ncomp; ks++){
	static const int sp[][2] = {{0,0}, {0,1}, {1,1}};
	static const int sp2[] = {0, 2, 5};
	FLOAT d2x_avg;
	int is, js, ks2;

	is = sp[ks][0]; js = sp[ks][1];
	
	d2x_avg = -2.0*x[is]*x[js]*sfact/(x_avg*x_avg);
	if(is == js)
	  d2x_avg += 11.0;
	d2x_avg *= 2.0*x[is]*x[js]/(9.0*sfact*x_avg*ds[is]*ds[js]);

	v2rho2[ks] += f_LDA_opp[ks]*g_ab + 
	  dg_ab*(v_LDA_opp[is]*dx_avg[js] + v_LDA_opp[js]*dx_avg[is]) +
	  e_LDA_opp*(d2g_ab*dx_avg[is]*dx_avg[js] + dg_ab*d2x_avg);

	ks2 = sp2[ks];

	v2sigma2[ks2] += e_LDA_opp*POW(ds[is], -8.0/3.0)*POW(ds[js], -8.0/3.0) *
	  (d2g_ab - dg_ab/x_avg)/(sfact2*16.0*x_avg*x_avg);
      }

      ncomp = (p->nspin == XC_POLARIZED) ? 6 : 1;
      for(ks=0; ks<ncomp; ks++){
	static const int sp[][2] = {{0,0}, {-1,-1}, {0,1}, {1,0}, {-1,-1}, {1,1}};
	FLOAT tmp1;
	int is, js;

	is = sp[ks][0]; js = sp[ks][1];
	if(is==-1) continue;

	tmp1 = x[js]/(2.0*sigma[js==0 ? 0 : 2]);
	if(is == js)
	  v2rhosigma[ks] += -4.0*e_LDA_opp*dg_ab*x[is]*tmp1/(3.0*x_avg*ds[is]);

	v2rhosigma[ks] += tmp1*x[js]/(2.0*x_avg)*
	  (v_LDA_opp[is]*dg_ab - e_LDA_opp*2.0*x[is]*x[is]/(3.0*x_avg*ds[is])*(d2g_ab - dg_ab/x_avg));
      }
    }
  }

  if(zk != NULL)
    *zk /= dens; /* we want energy per particle */
}
