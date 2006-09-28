#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

/************************************************************************
 Implements Becke 88 Generalized Gradient Approximation.
************************************************************************/

static func_type func_gga_x_b88 = {
  XC_GGA_X_B88,
  XC_EXCHANGE,
  "Becke 88",
  "GGA",
  "A. D. Becke, Phys. Rev. A 38, 3098-3100 (1988)"
};

/* some parameters */
static const double beta  = 0.0042;

void gga_x_b88_init(gga_type *p)
{
  p->func = &func_gga_x_b88;
  p->lda_aux = (lda_type *) malloc(sizeof(lda_type));
  lda_x_init(p->lda_aux, p->nspin, 3, XC_NON_RELATIVISTIC);
#if defined(LDA_SPEEDUP)
  lda_x_speedup(p->lda_aux, p->nspin, 3, XC_NON_RELATIVISTIC);
#endif
}

void gga_x_b88_end(gga_type *p)
{
  free(p->lda_aux);
}


void gga_x_b88(gga_type *p, double *rho, double *sigma,
	       double *e, double *vrho, double *vsigma)
{
  double sfact, dens, e1;
  int is;

  *e   = 0.0;
  if(p->nspin == XC_POLARIZED){
    sfact     = 1.0;
    vsigma[1] = 0.0; /* there are no cross terms in this functional */
  }else
    sfact     = 2.0;

  lda(p->lda_aux, rho, e, vrho, NULL);

  e1   = 0.0;
  dens = 0.0;
  for(is=0; is<p->nspin; is++){
    double gdm;
    int js = is==0 ? 0 : 2;

    dens += rho[is];
    gdm   = sqrt(sigma[js])/sfact;

    if(gdm>MIN_GRAD){
      double x, f1, f, ds, rho13;
      double dfdx;
      
      ds    = rho[is]/sfact;
      rho13 = pow(ds, 1.0/3.0);
      x     =  gdm/(ds*rho13);
      f1 = (1.0 + 6.0*beta*x*asinh(x));
      f  = x*x/f1;

      e1   -= sfact*beta*(ds*rho13)*f;
    
      dfdx = x*(2.0 + 6.0*beta*(x*asinh(x) - x*x/sqrt(1.0+x*x)))/(f1*f1);
      printf("%f %f\n", x, sfact*beta*rho13);
      
      vrho[is]  += -4.0/3.0*beta*rho13*(f - dfdx*x);
      vsigma[js] = -sfact*beta*(ds*rho13)*dfdx*x/(2.0*sigma[js]);

    }else{
      vsigma[is==0 ? 0 : 2] = 0.0;
    }
  }

  *e += e1/dens; /* we want energy per particle */
}
