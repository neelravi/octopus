/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

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

 $Id$
*/

#include <config.h>

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_multimin.h>


/* ---------------------- Interface to GSL functions ------------------------ */


/* Mathematical Functions */
double FC_FUNC_(oct_asinh, OCT_ASINH)
		 (double *x)
{
  return gsl_asinh(*x);
}

/* Special Functions */
double FC_FUNC_(oct_gamma, OCT_GAMMA)
		 (double *x)
{
  return gsl_sf_gamma(*x);
}

double FC_FUNC_(oct_hypergeometric, OCT_HYPERGEOMETRIC)
                 (double *a, double*b, double *x)
{
  return gsl_sf_hyperg_U(*a, *b, *x);
}

double FC_FUNC_(oct_incomplete_gamma, OCT_INCOMPLETE_GAMMA)
		 (double *a, double *x)
{
  return gsl_sf_gamma_inc_Q(*a, *x);
}

double FC_FUNC_(oct_sph_bessel, OCT_SPH_BESSEL)
     (int *l, double*x)
{
  return gsl_sf_bessel_jl(*l, *x);
}

double FC_FUNC_(oct_bessel, OCT_BESSEL)
     (int *n, double *x)
{
  return gsl_sf_bessel_Jn(*n, *x);
}

double FC_FUNC_(oct_bessel_in, OCT_BESSEL_IN)
     (int *n, double *x)
{
  return gsl_sf_bessel_In(*n, *x);
}

double FC_FUNC_(oct_bessel_j0, OCT_BESSEL_J0)
     (double *x)
{
  return gsl_sf_bessel_J0(*x);
}

double FC_FUNC_(oct_bessel_j1, OCT_BESSEL_J1)
     (double *x)
{
  return gsl_sf_bessel_J1(*x);
}

double FC_FUNC_(oct_bessel_k0, OCT_BESSEL_K0)
     (double *x)
{
  return gsl_sf_bessel_K0(*x);
}

double FC_FUNC_(oct_bessel_k1, OCT_BESSEL_K1)
     (double *x)
{
  return gsl_sf_bessel_K1(*x);
}

double FC_FUNC_(oct_erfc, OCT_ERFC)
		 (double *x)
{
  /* avoid floating invalids in the asymptotic limit */
  if(*x >  20.0) return  0.0;
	if(*x < -20.0) return  2.0;
	/* otherwise call gsl */
	return gsl_sf_erfc(*x);
}

double FC_FUNC_(oct_erf, OCT_ERF)
		 (double *x)
{
  /* avoid floating invalids in the asymptotic limit */
  if(*x >  20.0) return  1.0;
	if(*x < -20.0) return -1.0;
	/* otherwise call gsl */
	return gsl_sf_erf(*x);
}

double FC_FUNC_(oct_legendre_sphplm, OCT_LEGENDRE_SPHPLM)
		 (int *l, int *m, double *x)
{
  return gsl_sf_legendre_sphPlm(*l, *m, *x);
}


/* Vectors and Matrices */


/* Permutations */


/* Linear Algebra */


/* Random Number Generation */
void FC_FUNC_(oct_ran_init, OCT_RAN_INIT)
     (gsl_rng **r)
{
  gsl_rng_env_setup();
  *r = gsl_rng_alloc(gsl_rng_default);
}

void FC_FUNC_(oct_ran_end, OCT_RAN_END)
     (gsl_rng **r)
{
  gsl_rng_free(*r);
}


/* Random Number Distributions */ 
double FC_FUNC_(oct_ran_gaussian, OCT_RAN_GAUSSIAN)
		(gsl_rng **r, double *sigma)
{
  return gsl_ran_gaussian(*r, *sigma);
}


/* Interpolation */
void FC_FUNC_(oct_spline_end, OCT_SPLINE_END)
		 (void **spl, void **acc)
{
	gsl_spline_free((gsl_spline *)(*spl));
	gsl_interp_accel_free((gsl_interp_accel *)(*acc));
}

void FC_FUNC_(oct_spline_fit, OCT_SPLINE_FIT)
		 (int *nrc, double *x, double *y, void **spl, void **acc)
{
	*acc = (void *)gsl_interp_accel_alloc();
	*spl = (void *)gsl_spline_alloc(gsl_interp_cspline, *nrc);	
	gsl_spline_init((gsl_spline *)(*spl), x, y, *nrc);
	fflush(stdout);
}

double FC_FUNC_(oct_spline_eval, OCT_SPLINE_EVAL)
		 (double *x, void **spl, void **acc)
{
	return gsl_spline_eval((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}


/* This function returns the number of points with which a spline
	 was constructed (the size component of the gsl_spline struct). */
int FC_FUNC_(oct_spline_npoints, OCT_SPLINE_NPOINTS)
		 (void **spl)
{
	return (int)((gsl_spline *)(*spl))->size;
}

/* This function places in the x array the x values of a given spline spl*/ 
void FC_FUNC_(oct_spline_x, OCT_SPLINE_X)
     (void **spl, double *x)
{
  int size, i;
	
  size = (int)((gsl_spline *)(*spl))->size;
  for(i=0; i<size; i++)
		x[i] = ((gsl_spline *)(*spl))->x[i];
}

/* This function places in the y array the y values of a given spline spl*/ 
void FC_FUNC_(oct_spline_y, OCT_SPLINE_Y)
     (void **spl, double *y)
{
  int size, i;
	
  size = (int)((gsl_spline *)(*spl))->size;
  for(i=0; i<size; i++)
		y[i] = ((gsl_spline *)(*spl))->y[i];
}

/* Returns the integral of the spline stored in spl, between a and b */
double FC_FUNC_(oct_spline_eval_integ, OCT_SPLINE_EVAL_INTEG)
     (void **spl, double *a, double *b, void **acc)
{
  return gsl_spline_eval_integ((gsl_spline *)(*spl), *a, *b, (gsl_interp_accel *)(* acc));
}

/* Performs the derivative of a spline */
double FC_FUNC_(oct_spline_eval_der, OCT_SPLINE_EVAL_DER)
     (double *x, void **spl, void **acc)
{
  return gsl_spline_eval_deriv((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}

/* Performs the second derivative of a spline */
double FC_FUNC_(oct_spline_eval_der2, OCT_SPLINE_EVAL_DER2)
     (double *x, void **spl, void **acc)
{
  return gsl_spline_eval_deriv2((gsl_spline *)(*spl), *x, (gsl_interp_accel *)(*acc));
}


double my_f (const gsl_vector *v, void *params)
{
  double val;
  double *x, *gradient;
  int i, dim, getgrad;
  void (*para)(int*, double*, double*, int*, double*) = params;

  dim = v->size;
  x = (double *)malloc(dim*sizeof(double));

  for(i=0; i<dim; i++) x[i] = gsl_vector_get(v, i);
  getgrad = 0;
  para(&dim, x, &val, &getgrad, gradient);

  free(x);
  return val;
}

     /* The gradient of f, df = (df/dx, df/dy). */
void my_df (const gsl_vector *v, void *params,
            gsl_vector *df)
{
  double val;
  double *x, *gradient;
  int i, dim, getgrad;
  void (*para)(int*, double*, double*, int*, double*) = params;

  dim = v->size;
  x = (double *)malloc(dim*sizeof(double));
  gradient = (double *)malloc(dim*sizeof(double));

  for(i=0; i<dim; i++) x[i] = gsl_vector_get(v, i);
  getgrad = 1;
  para(&dim, x, &val, &getgrad, gradient);
  for(i=0; i<dim; i++) gsl_vector_set(df, i, gradient[i]);

  free(x); free(gradient);
}

/* Compute both f and df together. */
void my_fdf (const gsl_vector *x, void *params,
             double *f, gsl_vector *df)
{
  *f = my_f(x, params);
  my_df(x, params, df);
}

double FC_FUNC_(oct_minimize, OCT_MINIMIZE)
     (int*method, int *dim, double *point, double *step, double *tol, int *maxiter, void *f)
{

  size_t iter = 0;
  int status;
  double return_value, maxgrad;
  int i;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  gsl_vector *x;
  gsl_vector *grad;
  gsl_vector *absgrad;
  gsl_multimin_function_fdf my_func;

  my_func.f = &my_f;
  my_func.df = &my_df;
  my_func.fdf = &my_fdf;
  my_func.n = *dim;
  my_func.params = f;

  /* Starting point */
  x = gsl_vector_alloc (*dim);
  for(i=0; i<*dim; i++) gsl_vector_set (x, i, point[i]);

  /* Allocate space for the gradient */
  grad = gsl_vector_alloc (*dim);
  absgrad = gsl_vector_alloc (*dim);

  switch(*method){
  case 1: T = gsl_multimin_fdfminimizer_steepest_descent;
  case 2: T = gsl_multimin_fdfminimizer_conjugate_fr;
  case 3: T = gsl_multimin_fdfminimizer_conjugate_pr;
  case 4: T = gsl_multimin_fdfminimizer_vector_bfgs;
  }
  s = gsl_multimin_fdfminimizer_alloc (T, *dim);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, *step, 1e-4);
  do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      for(i=0; i<*dim; i++) point[i] = gsl_vector_get(gsl_multimin_fdfminimizer_gradient(s), i);
      for(i=0; i<*dim; i++) gsl_vector_set(absgrad, i, fabs(point[i]));
      maxgrad = gsl_vector_max(absgrad);

      for(i=0; i<*dim; i++) point[i] = gsl_vector_get(gsl_multimin_fdfminimizer_x(s), i);
      return_value = gsl_multimin_fdfminimizer_minimum(s);

      if (status) break;

      status = (maxgrad > *tol);
      }
  while (status && iter < *maxiter);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x); gsl_vector_free(grad); gsl_vector_free(absgrad);
  return return_value;
}


/* Chebyshev Approximations */
/*
void FC_FUNC_(oct_chebyshev_coeffs, OCT_CHEBYSHEV_COEFFS)
     (gsl_complex *coeffs, int *order)
{
  int i;
  double f (double x, void *p){return cos(x);}
  double g (double x, void *p){return -sin(x);}
  gsl_cheb_series *cs = gsl_cheb_alloc (*order);
  gsl_function F;
  F.function = f;
  F.params = 0;
  gsl_cheb_init (cs, &F, -1.0, 1.0);
  for(i=0; i<=12; i++){GSL_SET_REAL(&coeffs[i], (*cs).c[i]);}
  F.function = g;
  F.params = 0;
  gsl_cheb_init (cs, &F, -1.0, 1.0);
  for(i=0; i<=12; i++){GSL_SET_IMAG(&coeffs[i], (*cs).c[i]);}    
}
*/
