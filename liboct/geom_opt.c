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
*/

#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_min.h>

#include "config.h"

struct geom_oct_type_struct
{
	void (* func) (double *x, double *f, double *df);

	int size;
	double *x, f, *df;
};

typedef struct geom_oct_type_struct geom_oct_type;

static void my_work_fdf(const gsl_vector *x, geom_oct_type *p)
{
	int i;

	for(i=0; i<p->size; i++){
		if(gsl_vector_get(x, i) != p->x[i])
			break;
	}

	if(i<p->size){ /* recalculate */
		for(i=0; i<p->size; i++)
			p->x[i] = gsl_vector_get(x, i);
 
		(*(p->func))(p->x, &(p->f), p->df);
	}
}

static double my_f(const gsl_vector *x, void *params)
{
  geom_oct_type *p = (geom_oct_type *)params;
  
	my_work_fdf(x, p);
	return p->f;
}

static void my_df(const gsl_vector *x, void *params, gsl_vector *df)
{
	int i;
	geom_oct_type *p = (geom_oct_type *)params;
 
	my_work_fdf(x, p);
	for(i=0; i<p->size; i++)
		gsl_vector_set(df, i, p->df[i]);
}

static void my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df) 
{
	int i;
	geom_oct_type *p = (geom_oct_type *)params;
 
	my_work_fdf(x, p);
	*f = p->f;
	for(i=0; i<p->size; i++)
		gsl_vector_set(df, i, p->df[i]);
}

int F90_FUNC_(oct_geom_opt, OCT_GEOM_OPT)
		 (double *x, int *size, int *method, double *tol, int *max_iter,
			void (* cp) (double *, double *, double *))
{
	geom_oct_type c_geom_oct;
	size_t iter = 0;
  int i, status;
	double norm;
  gsl_multimin_fdfminimizer *s;
	gsl_multimin_function_fdf my_func;
  gsl_vector *vx;
	
	c_geom_oct.size = *size;
	c_geom_oct.x  = (double *) malloc((*size)*sizeof(double));
	c_geom_oct.df = (double *) malloc((*size)*sizeof(double));
	c_geom_oct.func = *cp;

	/* setup starting point */
  vx = gsl_vector_calloc(c_geom_oct.size);
	for(i = 0; i<c_geom_oct.size; i++){
		c_geom_oct.x[i] = 0.0;
		gsl_vector_set(vx, i, x[i]);
	}

	/* setup function */
  my_func.f      = &my_f;
  my_func.df     = &my_df;
  my_func.fdf    = &my_fdf;
  my_func.n      = c_geom_oct.size;
  my_func.params = (void *)(&c_geom_oct);
	
	switch(*method){
	case 1:
		s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_steepest_descent, c_geom_oct.size);
		break;
	case 2:
		s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_pr, c_geom_oct.size);
		break;
	case 3:
		s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, c_geom_oct.size);
		break;
	case 4:
		s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, c_geom_oct.size);
		break;
	}
		
  gsl_multimin_fdfminimizer_set (s, &my_func, vx, 0.01, *tol);

	printf("Info: Using %s minimiser\n", gsl_multimin_fdfminimizer_name(s));
	fflush(stdout);

  do{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);

		if (status)
			break;

		//norm = (double)gsl_blas_dnrm2(s->gradient);
		norm = 0.;
		for(i=0; i<s->gradient->size; i++)
			norm += s->gradient->data[i] * s->gradient->data[i];
		norm = sqrt(norm);

		status = (norm < *tol) ? GSL_SUCCESS : GSL_CONTINUE;

		printf("Info: geom_opt norm = %lf, tol = %lf\n", norm, *tol);
		fflush(stdout);

	}while (status == GSL_CONTINUE && iter < *max_iter);

	/* return minimum in x */
	for(i = 0; i<c_geom_oct.size; i++){
		x[i] = gsl_vector_get(s->x, i);
	}

	/* clean up */
	gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(vx);
	free(c_geom_oct.x);
	free(c_geom_oct.df);
	
	return (status == GSL_SUCCESS) ? 1 : 0;
}
