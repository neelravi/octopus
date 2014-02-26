/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include "util.h"

static FLOAT bi0_data[12] = {
  -.07660547252839144951, 1.92733795399380827000,  .22826445869203013390,  .01304891466707290428,  .00043442709008164874,
   .00000942265768600193,  .00000014340062895106,  .00000000161384906966,  .00000000001396650044,  .00000000000009579451,
   .00000000000000053339,  .00000000000000000245
};

static FLOAT ai0_data[21] = {
   .07575994494023796,  .00759138081082334,  .00041531313389237,  .00001070076463439, -.00000790117997921,
  -.00000078261435014,  .00000027838499429,  .00000000825247260, -.00000001204463945,  .00000000155964859,
   .00000000022925563, -.00000000011916228,  .00000000001757854,  .00000000000112822, -.00000000000114684,
   .00000000000027155, -.00000000000002415, -.00000000000000608,  .00000000000000314, -.00000000000000071,
   .00000000000000007
};

static FLOAT ai02_data[22] = {
   .05449041101410882,  .00336911647825569,  .00006889758346918,  .00000289137052082,  .00000020489185893,
   .00000002266668991,  .00000000339623203,  .00000000049406022,  .00000000001188914, -.00000000003149915,
  -.00000000001321580, -.00000000000179419,  .00000000000071801,  .00000000000038529,  .00000000000001539,
  -.00000000000004151, -.00000000000000954,  .00000000000000382,  .00000000000000176, -.00000000000000034,
  -.00000000000000027,  .00000000000000003
};

/* 
   Compute the exponentially scaled modified (hyperbolic)
   Bessel function of the first kind of order zero.

   based on the SLATEC routine by W. Fullerton
*/
FLOAT XC(bessel_I0_scaled)(const FLOAT x)
{
  FLOAT y = ABS(x), r = 0.0;

  if(y < 2.0*SQRT_FLOAT_EPSILON)
    r = 1.0 - y;
  else if(y <= 3.0)
    r = EXP(-y)*(2.75 + XC(cheb_eval)(y*y/4.5-1.0, bi0_data, 12));
  else if(y <= 8.0)
    r = (.375 + XC(cheb_eval)((48.0/y - 11.0)/5.0, ai0_data, 21))/SQRT(y);
  else
    r = (.375 + XC(cheb_eval)(16.0/y - 1.0, ai02_data, 22))/SQRT(y);

  return r;
}


/*
  Compute the hyperbolic Bessel function of the first kind
  of order zero.

  based on the SLATEC routine by W. Fullerton
*/
FLOAT XC(bessel_I0)(const FLOAT x)
{
  FLOAT y = fabs(x), r = 0.0;

  if(y < 2.0*SQRT_FLOAT_EPSILON)
    r = 1.0;
  else if(y <= 3.0)
    r = 2.75 + XC(cheb_eval)(y*y/4.5 - 1.0, bi0_data, 12);
  else if(y < LOG_FLOAT_MAX - 1.0)
    r = EXP(y) * XC(bessel_I0_scaled)(x);
  else
    fprintf(stderr, "Overflow in bessel_I0\n");

  return r;
}

static FLOAT bi1_data[11] = {
  -0.001971713261099859,  0.407348876675464810,  0.034838994299959456,  0.001545394556300123,  0.000041888521098377,
   0.000000764902676483,  0.000000010042493924,  0.000000000099322077,  0.000000000000766380,  0.000000000000004741,
   0.000000000000000024
};

static FLOAT ai1_data[21] = {
  -0.02846744181881479, -0.01922953231443221, -0.00061151858579437, -0.00002069971253350,  0.00000858561914581,
   0.00000104949824671, -0.00000029183389184, -0.00000001559378146,  0.00000001318012367, -0.00000000144842341,
  -0.00000000029085122,  0.00000000012663889, -0.00000000001664947, -0.00000000000166665,  0.00000000000124260,
  -0.00000000000027315,  0.00000000000002023,  0.00000000000000730, -0.00000000000000333,  0.00000000000000071,
  -0.00000000000000006
};

static FLOAT ai12_data[22] = {
   0.02857623501828014, -0.00976109749136147, -0.00011058893876263, -0.00000388256480887, -0.00000025122362377,
  -0.00000002631468847, -0.00000000383538039, -0.00000000055897433, -0.00000000001897495,  0.00000000003252602,
   0.00000000001412580,  0.00000000000203564, -0.00000000000071985, -0.00000000000040836, -0.00000000000002101,
   0.00000000000004273,  0.00000000000001041, -0.00000000000000382, -0.00000000000000186,  0.00000000000000033,
   0.00000000000000028, -0.00000000000000003
};


FLOAT XC(bessel_I1_scaled)(const FLOAT x)
{
  const FLOAT xmin    = 2.0 * FLOAT_MIN;
  const FLOAT x_small = 2.0 * M_SQRT2 * SQRT_FLOAT_EPSILON;
  const FLOAT y = fabs(x);
  FLOAT r = 0.0;

  if(y == 0.0)
    r = 0.0;
  else if(y < xmin)
    fprintf(stderr, "Underflow error in bessel_I1_scaled\n");
  else if(y < x_small)
    r = 0.5*x*EXP(-y);
  else if(y <= 3.0)
    r = x*EXP(-y)*(0.875 + XC(cheb_eval)(y*y/4.5 - 1.0, bi1_data, 11));
  else{
    if(y <= 8.0)
      r = (0.375 + XC(cheb_eval)((48.0/y - 11.0)/5.0, ai1_data, 21))/SQRT(y);
    else
      r = (0.375 + XC(cheb_eval)(16.0/y - 1.0, ai12_data, 22))/SQRT(y);

    r *= (x > 0.0 ? 1.0 : -1.0);
  }

  return r;
}

FLOAT XC(bessel_I1)(const FLOAT x)
{
  const FLOAT xmin    = 2.0 * FLOAT_MIN;
  const FLOAT x_small = 2.0 * M_SQRT2 * SQRT_FLOAT_EPSILON;
  const FLOAT y = fabs(x);
  FLOAT r = 0.0;

  if(y == 0.0)
    r = 0.0;
  else if(y < xmin)
    fprintf(stderr, "Underflow error in bessel_I1\n");
  else if(y < x_small)
    r = 0.5*x;
  else if(y <= 3.0)
    r = x*(0.875 + XC(cheb_eval)(y*y/4.5 - 1.0, bi1_data, 11));
  else
    r = EXP(x)*XC(bessel_I1_scaled(x));

  return r;
}

static FLOAT bk0_data[11] = {
  -0.03532739323390276872,  0.3442898999246284869,   0.03597993651536150163,  0.00126461541144692592,  0.00002286212103119451,
   0.00000025347910790261,  0.00000000190451637722,  0.00000000001034969525,  0.00000000000004259816,  0.00000000000000013744,
   0.00000000000000000035
};

static FLOAT ak0_data[17] = {
  -0.07643947903327941, -0.02235652605699819,  0.00077341811546938, -0.00004281006688886,  0.00000308170017386,
  -0.00000026393672220,  0.00000002563713036, -0.00000000274270554,  0.00000000031694296, -0.00000000003902353,
   0.00000000000506804, -0.00000000000068895,  0.00000000000009744, -0.00000000000001427,  0.00000000000000215,
  -0.00000000000000033,  0.00000000000000005
};

static FLOAT ak02_data[14] = {
  -0.01201869826307592, -0.00917485269102569,  0.00014445509317750, -0.00000401361417543,  0.00000015678318108,
  -0.00000000777011043,  0.00000000046111825, -0.00000000003158592,  0.00000000000243501, -0.00000000000020743,
   0.00000000000001925, -0.00000000000000192,  0.00000000000000020, -0.00000000000000002
};


FLOAT XC(bessel_K0_scaled)(const FLOAT x)
{
  FLOAT r = 0.0;

  if(x <= 0.0)
    fprintf(stderr, "Domain error in bessel_K0_scaled\n");
  else if(x <= 2.0)
    r = EXP(x)*(-LOG(0.5*x)*XC(bessel_I0)(x) - 0.25 + XC(cheb_eval)(0.5*x*x - 1.0, bk0_data, 11));
  else if(x <= 8.0)
    r = (1.25 + XC(cheb_eval)((16.0/x - 5.0)/3.0, ak0_data, 17))/SQRT(x);
  else
    r = (1.25 + XC(cheb_eval)(16.0/x - 1.0, ak02_data, 14))/SQRT(x);

  return r;
}

FLOAT XC(bessel_K0)(const FLOAT x)
{
  FLOAT r = 0.0;

  if(x <= 0.0)
    fprintf(stderr, "Domain error in bessel_K0\n");
  else if(x <= 2.0)
    r = -LOG(0.5*x)*XC(bessel_I0)(x) - 0.25 + XC(cheb_eval)(0.5*x*x - 1.0, bk0_data, 11);
  else
    r = EXP(-x)*XC(bessel_K0_scaled)(x);

  return r;
}


static FLOAT bk1_data[11] = {
   0.0253002273389477705, -0.3531559607765448760, -0.1226111808226571480, -0.0069757238596398643, -0.0001730288957513052,
  -0.0000024334061415659, -0.0000000221338763073, -0.0000000001411488392, -0.0000000000006666901, -0.0000000000000024274,
  -0.0000000000000000070
};

static FLOAT ak1_data[17] = {
   0.27443134069738830,  0.07571989953199368, -0.00144105155647540,  0.00006650116955125, -0.00000436998470952,
   0.00000035402774997, -0.00000003311163779,  0.00000000344597758, -0.00000000038989323,  0.00000000004720819,
  -0.00000000000604783,  0.00000000000081284, -0.00000000000011386,  0.00000000000001654, -0.00000000000000248,
   0.00000000000000038, -0.00000000000000006
};

static FLOAT ak12_data[14] = {
   0.06379308343739001,  0.02832887813049721, -0.00024753706739052,  0.00000577197245160, -0.00000020689392195,
   0.00000000973998344, -0.00000000055853361,  0.00000000003732996, -0.00000000000282505,  0.00000000000023720,
  -0.00000000000002176,  0.00000000000000215, -0.00000000000000022,  0.00000000000000002
};

FLOAT XC(bessel_K1_scaled)(const FLOAT x)
{
  FLOAT r = 0.0;

  if(x <= 0.0)
    fprintf(stderr, "Domain error in bessel_K1_scaled\n");
  else if(x <= 2.0)
    r =  EXP(x)*(LOG(0.5*x)*XC(bessel_I1)(x) +
		 (0.75 + XC(cheb_eval)(.5*x*x - 1.0, bk1_data, 11))/x);
  else if(x <= 8.0)
    r = (1.25 + XC(cheb_eval)((16.0/x - 5.0)/3.0, ak1_data, 17))/SQRT(x);
  else
    r = (1.25 + XC(cheb_eval)(16.0/x - 1.0, ak12_data, 14))/SQRT(x);

  return r;
}


FLOAT XC(bessel_K1)(const FLOAT x)
{
  FLOAT r = 0.0;

  if(x <= 0.0)
    fprintf(stderr, "Domain error in bessel_K1\n");
  else if(x<2.0*FLOAT_MIN)
    fprintf(stderr, "Overflow error in bessel_K1\n");
  else if(x <= 2.0)
    r = LOG(0.5*x)*XC(bessel_I1)(x) + (0.75 + XC(cheb_eval)(0.5*x*x - 1.0, bk1_data, 11))/x;
  else
    r = EXP(-x)*XC(bessel_K1_scaled)(x);

  return r;
}
