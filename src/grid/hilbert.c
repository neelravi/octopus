/*
 Copyright (C) 2013 X. Andrade

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

#include <assert.h>

/* This code is based on the implementation of the Hilbert curves
   presented in:

   J. Skilling, Programming the Hilbert curve, AIP Conf. Proc. 707, 381 (2004); http://dx.doi.org/10.1063/1.1751381

*/

void InttoTranspose(const long long int h, int * x){
  /* the code uses some funny way of storing the bits */

  x[2] =  (((h>>0)&1)<<0)   + (((h>>3)&1)<<1)   + (((h>>6)&1)<<2)   + (((h>>9 )&1)<<3)  + (((h>>12)&1)<<4);
  x[1] =  (((h>>1)&1)<<0)   + (((h>>4)&1)<<1)   + (((h>>7)&1)<<2)   + (((h>>10)&1)<<3)  + (((h>>13)&1)<<4);
  x[0] =  (((h>>2)&1)<<0)   + (((h>>5)&1)<<1)   + (((h>>8)&1)<<2)   + (((h>>11)&1)<<3)  + (((h>>14)&1)<<4);

  x[2] += (((h>>15)&1)<<5)  + (((h>>18)&1)<<6)  + (((h>>21)&1)<<7)  + (((h>>24)&1)<<8)  + (((h>>27)&1)<<9);
  x[1] += (((h>>16)&1)<<5)  + (((h>>19)&1)<<6)  + (((h>>22)&1)<<7)  + (((h>>25)&1)<<8)  + (((h>>28)&1)<<9);
  x[0] += (((h>>17)&1)<<5)  + (((h>>20)&1)<<6)  + (((h>>23)&1)<<7)  + (((h>>26)&1)<<8)  + (((h>>29)&1)<<9);

  x[2] += (((h>>30)&1)<<10) + (((h>>33)&1)<<11) + (((h>>36)&1)<<12) + (((h>>39)&1)<<13) + (((h>>42)&1)<<14);
  x[1] += (((h>>31)&1)<<10) + (((h>>34)&1)<<11) + (((h>>37)&1)<<12) + (((h>>40)&1)<<13) + (((h>>43)&1)<<14);
  x[0] += (((h>>32)&1)<<10) + (((h>>35)&1)<<11) + (((h>>38)&1)<<12) + (((h>>41)&1)<<13) + (((h>>44)&1)<<14);
  
  x[2] += (((h>>45)&1)<<15) + (((h>>48)&1)<<16) + (((h>>51)&1)<<17) + (((h>>54)&1)<<18) + (((h>>57)&1)<<19);
  x[1] += (((h>>46)&1)<<15) + (((h>>49)&1)<<16) + (((h>>52)&1)<<17) + (((h>>55)&1)<<18) + (((h>>58)&1)<<19);
  x[0] += (((h>>47)&1)<<15) + (((h>>50)&1)<<16) + (((h>>53)&1)<<17) + (((h>>56)&1)<<18) + (((h>>59)&1)<<19);
  
  x[2] += (((h>>60)&1)<<20);
  x[1] += (((h>>61)&1)<<20);
  x[0] += (((h>>62)&1)<<20);
  
}


void TransposetoAxes(int* X, int b, int n ){ /* position, #bits, dimension */

  int N = 2 << (b-1), P, Q, t;
  int i;

  /* Gray decode by H ^ (H/2) */
  t = X[n-1] >> 1;
  for(i = n - 1; i > 0; i--) X[i] ^= X[i-1];
  X[0] ^= t;

  /* Undo excess work */
  for( Q = 2; Q != N; Q <<= 1 ) {
    P = Q - 1;
    for( i = n-1; i >= 0 ; i-- ){
      if( X[i] & Q ) {
	X[0] ^= P; /* invert */
      } else{ 
	t = (X[0]^X[i]) & P; X[0] ^= t; X[i] ^= t; /* exchange */
      }
    }
  }

}

void FC_FUNC_(hilbert_index_to_point, HILBERT_INDEX_TO_POINT)(const int * nbits, const long long int * index, int * point){
  
  InttoTranspose(*index, point);
  TransposetoAxes(point, *nbits, 3);
  
}
