#include <math.h>
#include <stdio.h>

#include "xc.h"

void main()
{
	lda_type l1, l2, l3;
	double rs;

	lda_init(&l1, XC_LDA_C_OB_PZ, XC_POLARIZED);
	lda_init(&l2, XC_LDA_C_GL, XC_POLARIZED);
	lda_init(&l3, XC_LDA_C_VWN, XC_POLARIZED);

	for(rs=0.001; rs<10; rs+=0.1){
		double ec, vc[2], rho[2];
		
		rho[0] = 1.0/(4.0/3.0*M_PI*pow(rs,3));

		rho[1] = 0.9*rho[0];
		rho[0] = 0.1*rho[0];

		lda(&l1, rho, &ec, vc);
		printf("%lf\t%lf\t%lf", rs, ec, vc[1]);

		lda(&l2, rho, &ec, vc);
		printf("\t%lf\t%lf", ec, vc[1]);

		lda(&l3, rho, &ec, vc);
		printf("\t%lf\t%lf\n", ec, vc[1]);


	}
}
