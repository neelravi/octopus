/*
 * Copyright (c) 2006 The Trustees of Indiana University and Indiana
 *                    University Research and Technology
 *                    Corporation.  All rights reserved.
 * Copyright (c) 2006 The Technical University of Chemnitz. All 
 *                    rights reserved.
 */
#include "nbc.h"

/****************** THIS FILE is automatically generated *********************
 * changes will be deleted at the next generation of this file - see nbc_op.c.m4 */

int NBC_Operation(void *buf3, void *buf1, void *buf2, MPI_Op op, MPI_Datatype type, int count) {
  int i;
     
  if(type == MPI_INT) { 
    if(op == MPI_MIN) {
      for(i=0; i<count; i++) {
        if(*(((int*)buf1) + i) > *(((int*)buf2) + i)) *(((int*)buf3) + i) = *(((int*)buf2) + i); else *(((int*)buf3) + i) = *(((int*)buf1) + i); 
      }
    } else if(op == MPI_MAX) {
      for(i=0; i<count; i++) {
        if(*(((int*)buf1) + i) < *(((int*)buf2) + i)) *(((int*)buf3) + i) = *(((int*)buf2) + i); else *(((int*)buf3) + i) = *(((int*)buf1) + i); 
      }
    } else if(op == MPI_SUM) {
      for(i=0; i<count; i++) {
        *(((int*)buf3) + i) = *(((int*)buf1) + i) + *(((int*)buf2) + i); 
      }
    } else if(op == MPI_PROD) {
      for(i=0; i<count; i++) {
        *(((int*)buf3) + i) = *(((int*)buf1) + i) * *(((int*)buf2) + i); 
      }
    } else if(op == MPI_LAND) {
      for(i=0; i<count; i++) {
        *(((int*)buf3) + i) = *(((int*)buf1) + i) && *(((int*)buf2) + i); 
      }
    } else if(op == MPI_BAND) {
      for(i=0; i<count; i++) {
        *(((int*)buf3) + i) = *(((int*)buf1) + i) & *(((int*)buf2) + i); 
      }
    } else if(op == MPI_LOR) {
      for(i=0; i<count; i++) {
        *(((int*)buf3) + i) = *(((int*)buf1) + i) || *(((int*)buf2) + i); 
      }
    } else if(op == MPI_BOR) {
      for(i=0; i<count; i++) {
        *(((int*)buf3) + i) = *(((int*)buf1) + i) | *(((int*)buf2) + i); 
      }
    } else if(op == MPI_LXOR) {
      for(i=0; i<count; i++) {
        *(((int*)buf3) + i) = ((*(((int*)buf1) + i) ? 1 : 0) ^ (*(((int*)buf2) + i) ?  1 : 0)); 
      }
    } else if(op == MPI_BXOR) {
      for(i=0; i<count; i++) {
        *(((int*)buf3) + i) = ((*(((int*)buf1) + i)) ^ (*(((int*)buf2) + i))); 
      }
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_LONG) { 
    if(op == MPI_MIN) {
      for(i=0; i<count; i++) {
        if(*(((long*)buf1) + i) > *(((long*)buf2) + i)) *(((long*)buf3) + i) = *(((long*)buf2) + i); else *(((long*)buf3) + i) = *(((long*)buf1) + i); 
      }
    } else if(op == MPI_MAX) {
      for(i=0; i<count; i++) {
        if(*(((long*)buf1) + i) < *(((long*)buf2) + i)) *(((long*)buf3) + i) = *(((long*)buf2) + i); else *(((long*)buf3) + i) = *(((long*)buf1) + i); 
      }
    } else if(op == MPI_SUM) {
      for(i=0; i<count; i++) {
        *(((long*)buf3) + i) = *(((long*)buf1) + i) + *(((long*)buf2) + i); 
      }
    } else if(op == MPI_PROD) {
      for(i=0; i<count; i++) {
        *(((long*)buf3) + i) = *(((long*)buf1) + i) * *(((long*)buf2) + i); 
      }
    } else if(op == MPI_LAND) {
      for(i=0; i<count; i++) {
        *(((long*)buf3) + i) = *(((long*)buf1) + i) && *(((long*)buf2) + i); 
      }
    } else if(op == MPI_BAND) {
      for(i=0; i<count; i++) {
        *(((long*)buf3) + i) = *(((long*)buf1) + i) & *(((long*)buf2) + i); 
      }
    } else if(op == MPI_LOR) {
      for(i=0; i<count; i++) {
        *(((long*)buf3) + i) = *(((long*)buf1) + i) || *(((long*)buf2) + i); 
      }
    } else if(op == MPI_BOR) {
      for(i=0; i<count; i++) {
        *(((long*)buf3) + i) = *(((long*)buf1) + i) | *(((long*)buf2) + i); 
      }
    } else if(op == MPI_LXOR) {
      for(i=0; i<count; i++) {
        *(((long*)buf3) + i) = ((*(((long*)buf1) + i) ? 1 : 0) ^ (*(((long*)buf2) + i) ?  1 : 0)); 
      }
    } else if(op == MPI_BXOR) {
      for(i=0; i<count; i++) {
        *(((long*)buf3) + i) = ((*(((long*)buf1) + i)) ^ (*(((long*)buf2) + i))); 
      }
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_SHORT) { 
    if(op == MPI_MIN) {
      for(i=0; i<count; i++) {
        if(*(((short*)buf1) + i) > *(((short*)buf2) + i)) *(((short*)buf3) + i) = *(((short*)buf2) + i); else *(((short*)buf3) + i) = *(((short*)buf1) + i); 
      }
    } else if(op == MPI_MAX) {
      for(i=0; i<count; i++) {
        if(*(((short*)buf1) + i) < *(((short*)buf2) + i)) *(((short*)buf3) + i) = *(((short*)buf2) + i); else *(((short*)buf3) + i) = *(((short*)buf1) + i); 
      }
    } else if(op == MPI_SUM) {
      for(i=0; i<count; i++) {
        *(((short*)buf3) + i) = *(((short*)buf1) + i) + *(((short*)buf2) + i); 
      }
    } else if(op == MPI_PROD) {
      for(i=0; i<count; i++) {
        *(((short*)buf3) + i) = *(((short*)buf1) + i) * *(((short*)buf2) + i); 
      }
    } else if(op == MPI_LAND) {
      for(i=0; i<count; i++) {
        *(((short*)buf3) + i) = *(((short*)buf1) + i) && *(((short*)buf2) + i); 
      }
    } else if(op == MPI_BAND) {
      for(i=0; i<count; i++) {
        *(((short*)buf3) + i) = *(((short*)buf1) + i) & *(((short*)buf2) + i); 
      }
    } else if(op == MPI_LOR) {
      for(i=0; i<count; i++) {
        *(((short*)buf3) + i) = *(((short*)buf1) + i) || *(((short*)buf2) + i); 
      }
    } else if(op == MPI_BOR) {
      for(i=0; i<count; i++) {
        *(((short*)buf3) + i) = *(((short*)buf1) + i) | *(((short*)buf2) + i); 
      }
    } else if(op == MPI_LXOR) {
      for(i=0; i<count; i++) {
        *(((short*)buf3) + i) = ((*(((short*)buf1) + i) ? 1 : 0) ^ (*(((short*)buf2) + i) ?  1 : 0)); 
      }
    } else if(op == MPI_BXOR) {
      for(i=0; i<count; i++) {
        *(((short*)buf3) + i) = ((*(((short*)buf1) + i)) ^ (*(((short*)buf2) + i))); 
      }
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_UNSIGNED) { 
    if(op == MPI_MIN) {
      for(i=0; i<count; i++) {
        if(*(((unsigned int*)buf1) + i) > *(((unsigned int*)buf2) + i)) *(((unsigned int*)buf3) + i) = *(((unsigned int*)buf2) + i); else *(((unsigned int*)buf3) + i) = *(((unsigned int*)buf1) + i); 
      }
    } else if(op == MPI_MAX) {
      for(i=0; i<count; i++) {
        if(*(((unsigned int*)buf1) + i) < *(((unsigned int*)buf2) + i)) *(((unsigned int*)buf3) + i) = *(((unsigned int*)buf2) + i); else *(((unsigned int*)buf3) + i) = *(((unsigned int*)buf1) + i); 
      }
    } else if(op == MPI_SUM) {
      for(i=0; i<count; i++) {
        *(((unsigned int*)buf3) + i) = *(((unsigned int*)buf1) + i) + *(((unsigned int*)buf2) + i); 
      }
    } else if(op == MPI_PROD) {
      for(i=0; i<count; i++) {
        *(((unsigned int*)buf3) + i) = *(((unsigned int*)buf1) + i) * *(((unsigned int*)buf2) + i); 
      }
    } else if(op == MPI_LAND) {
      for(i=0; i<count; i++) {
        *(((unsigned int*)buf3) + i) = *(((unsigned int*)buf1) + i) && *(((unsigned int*)buf2) + i); 
      }
    } else if(op == MPI_BAND) {
      for(i=0; i<count; i++) {
        *(((unsigned int*)buf3) + i) = *(((unsigned int*)buf1) + i) & *(((unsigned int*)buf2) + i); 
      }
    } else if(op == MPI_LOR) {
      for(i=0; i<count; i++) {
        *(((unsigned int*)buf3) + i) = *(((unsigned int*)buf1) + i) || *(((unsigned int*)buf2) + i); 
      }
    } else if(op == MPI_BOR) {
      for(i=0; i<count; i++) {
        *(((unsigned int*)buf3) + i) = *(((unsigned int*)buf1) + i) | *(((unsigned int*)buf2) + i); 
      }
    } else if(op == MPI_LXOR) {
      for(i=0; i<count; i++) {
        *(((unsigned int*)buf3) + i) = ((*(((unsigned int*)buf1) + i) ? 1 : 0) ^ (*(((unsigned int*)buf2) + i) ?  1 : 0)); 
      }
    } else if(op == MPI_BXOR) {
      for(i=0; i<count; i++) {
        *(((unsigned int*)buf3) + i) = ((*(((unsigned int*)buf1) + i)) ^ (*(((unsigned int*)buf2) + i))); 
      }
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_UNSIGNED_LONG) { 
    if(op == MPI_MIN) {
      for(i=0; i<count; i++) {
        if(*(((unsigned long*)buf1) + i) > *(((unsigned long*)buf2) + i)) *(((unsigned long*)buf3) + i) = *(((unsigned long*)buf2) + i); else *(((unsigned long*)buf3) + i) = *(((unsigned long*)buf1) + i); 
      }
    } else if(op == MPI_MAX) {
      for(i=0; i<count; i++) {
        if(*(((unsigned long*)buf1) + i) < *(((unsigned long*)buf2) + i)) *(((unsigned long*)buf3) + i) = *(((unsigned long*)buf2) + i); else *(((unsigned long*)buf3) + i) = *(((unsigned long*)buf1) + i); 
      }
    } else if(op == MPI_SUM) {
      for(i=0; i<count; i++) {
        *(((unsigned long*)buf3) + i) = *(((unsigned long*)buf1) + i) + *(((unsigned long*)buf2) + i); 
      }
    } else if(op == MPI_PROD) {
      for(i=0; i<count; i++) {
        *(((unsigned long*)buf3) + i) = *(((unsigned long*)buf1) + i) * *(((unsigned long*)buf2) + i); 
      }
    } else if(op == MPI_LAND) {
      for(i=0; i<count; i++) {
        *(((unsigned long*)buf3) + i) = *(((unsigned long*)buf1) + i) && *(((unsigned long*)buf2) + i); 
      }
    } else if(op == MPI_BAND) {
      for(i=0; i<count; i++) {
        *(((unsigned long*)buf3) + i) = *(((unsigned long*)buf1) + i) & *(((unsigned long*)buf2) + i); 
      }
    } else if(op == MPI_LOR) {
      for(i=0; i<count; i++) {
        *(((unsigned long*)buf3) + i) = *(((unsigned long*)buf1) + i) || *(((unsigned long*)buf2) + i); 
      }
    } else if(op == MPI_BOR) {
      for(i=0; i<count; i++) {
        *(((unsigned long*)buf3) + i) = *(((unsigned long*)buf1) + i) | *(((unsigned long*)buf2) + i); 
      }
    } else if(op == MPI_LXOR) {
      for(i=0; i<count; i++) {
        *(((unsigned long*)buf3) + i) = ((*(((unsigned long*)buf1) + i) ? 1 : 0) ^ (*(((unsigned long*)buf2) + i) ?  1 : 0)); 
      }
    } else if(op == MPI_BXOR) {
      for(i=0; i<count; i++) {
        *(((unsigned long*)buf3) + i) = ((*(((unsigned long*)buf1) + i)) ^ (*(((unsigned long*)buf2) + i))); 
      }
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_UNSIGNED_SHORT) { 
    if(op == MPI_MIN) {
      for(i=0; i<count; i++) {
        if(*(((unsigned short*)buf1) + i) > *(((unsigned short*)buf2) + i)) *(((unsigned short*)buf3) + i) = *(((unsigned short*)buf2) + i); else *(((unsigned short*)buf3) + i) = *(((unsigned short*)buf1) + i); 
      }
    } else if(op == MPI_MAX) {
      for(i=0; i<count; i++) {
        if(*(((unsigned short*)buf1) + i) < *(((unsigned short*)buf2) + i)) *(((unsigned short*)buf3) + i) = *(((unsigned short*)buf2) + i); else *(((unsigned short*)buf3) + i) = *(((unsigned short*)buf1) + i); 
      }
    } else if(op == MPI_SUM) {
      for(i=0; i<count; i++) {
        *(((unsigned short*)buf3) + i) = *(((unsigned short*)buf1) + i) + *(((unsigned short*)buf2) + i); 
      }
    } else if(op == MPI_PROD) {
      for(i=0; i<count; i++) {
        *(((unsigned short*)buf3) + i) = *(((unsigned short*)buf1) + i) * *(((unsigned short*)buf2) + i); 
      }
    } else if(op == MPI_LAND) {
      for(i=0; i<count; i++) {
        *(((unsigned short*)buf3) + i) = *(((unsigned short*)buf1) + i) && *(((unsigned short*)buf2) + i); 
      }
    } else if(op == MPI_BAND) {
      for(i=0; i<count; i++) {
        *(((unsigned short*)buf3) + i) = *(((unsigned short*)buf1) + i) & *(((unsigned short*)buf2) + i); 
      }
    } else if(op == MPI_LOR) {
      for(i=0; i<count; i++) {
        *(((unsigned short*)buf3) + i) = *(((unsigned short*)buf1) + i) || *(((unsigned short*)buf2) + i); 
      }
    } else if(op == MPI_BOR) {
      for(i=0; i<count; i++) {
        *(((unsigned short*)buf3) + i) = *(((unsigned short*)buf1) + i) | *(((unsigned short*)buf2) + i); 
      }
    } else if(op == MPI_LXOR) {
      for(i=0; i<count; i++) {
        *(((unsigned short*)buf3) + i) = ((*(((unsigned short*)buf1) + i) ? 1 : 0) ^ (*(((unsigned short*)buf2) + i) ?  1 : 0)); 
      }
    } else if(op == MPI_BXOR) {
      for(i=0; i<count; i++) {
        *(((unsigned short*)buf3) + i) = ((*(((unsigned short*)buf1) + i)) ^ (*(((unsigned short*)buf2) + i))); 
      }
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_FLOAT) { 
    if(op == MPI_MIN) {
      for(i=0; i<count; i++) {
        if(*(((float*)buf1) + i) > *(((float*)buf2) + i)) *(((float*)buf3) + i) = *(((float*)buf2) + i); else *(((float*)buf3) + i) = *(((float*)buf1) + i); 
      }
    } else if(op == MPI_MAX) {
      for(i=0; i<count; i++) {
        if(*(((float*)buf1) + i) < *(((float*)buf2) + i)) *(((float*)buf3) + i) = *(((float*)buf2) + i); else *(((float*)buf3) + i) = *(((float*)buf1) + i); 
      }
    } else if(op == MPI_SUM) {
      for(i=0; i<count; i++) {
        *(((float*)buf3) + i) = *(((float*)buf1) + i) + *(((float*)buf2) + i); 
      }
    } else if(op == MPI_PROD) {
      for(i=0; i<count; i++) {
        *(((float*)buf3) + i) = *(((float*)buf1) + i) * *(((float*)buf2) + i); 
      }
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_DOUBLE) { 
    if(op == MPI_MIN) {
      for(i=0; i<count; i++) {
        if(*(((double*)buf1) + i) > *(((double*)buf2) + i)) *(((double*)buf3) + i) = *(((double*)buf2) + i); else *(((double*)buf3) + i) = *(((double*)buf1) + i); 
      }
    } else if(op == MPI_MAX) {
      for(i=0; i<count; i++) {
        if(*(((double*)buf1) + i) < *(((double*)buf2) + i)) *(((double*)buf3) + i) = *(((double*)buf2) + i); else *(((double*)buf3) + i) = *(((double*)buf1) + i); 
      }
    } else if(op == MPI_SUM) {
      for(i=0; i<count; i++) {
        *(((double*)buf3) + i) = *(((double*)buf1) + i) + *(((double*)buf2) + i); 
      }
    } else if(op == MPI_PROD) {
      for(i=0; i<count; i++) {
        *(((double*)buf3) + i) = *(((double*)buf1) + i) * *(((double*)buf2) + i); 
      }
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_LONG_DOUBLE) { 
    if(op == MPI_MIN) {
      for(i=0; i<count; i++) {
        if(*(((long double*)buf1) + i) > *(((long double*)buf2) + i)) *(((long double*)buf3) + i) = *(((long double*)buf2) + i); else *(((long double*)buf3) + i) = *(((long double*)buf1) + i); 
      }
    } else if(op == MPI_MAX) {
      for(i=0; i<count; i++) {
        if(*(((long double*)buf1) + i) < *(((long double*)buf2) + i)) *(((long double*)buf3) + i) = *(((long double*)buf2) + i); else *(((long double*)buf3) + i) = *(((long double*)buf1) + i); 
      }
    } else if(op == MPI_SUM) {
      for(i=0; i<count; i++) {
        *(((long double*)buf3) + i) = *(((long double*)buf1) + i) + *(((long double*)buf2) + i); 
      }
    } else if(op == MPI_PROD) {
      for(i=0; i<count; i++) {
        *(((long double*)buf3) + i) = *(((long double*)buf1) + i) * *(((long double*)buf2) + i); 
      }
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_BYTE) { 
    if(op == MPI_BAND) {
      for(i=0; i<count; i++) {
        *(((char*)buf3) + i) = *(((char*)buf1) + i) & *(((char*)buf2) + i); 
      }
    } else if(op == MPI_BOR) {
      for(i=0; i<count; i++) {
        *(((char*)buf3) + i) = *(((char*)buf1) + i) | *(((char*)buf2) + i); 
      }
    } else if(op == MPI_BXOR) {
      for(i=0; i<count; i++) {
        *(((char*)buf3) + i) = ((*(((char*)buf1) + i)) ^ (*(((char*)buf2) + i))); 
      }
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_FLOAT_INT) { 
    if(op == MPI_MAXLOC) {
      for(i=0; i<count; i++) {
        typedef struct {
          float val;
          int rank;
        } float_int;
        float_int *ptr1, *ptr2, *ptr3;
                            
        ptr1 = ((float_int*)buf1) + i;
        ptr2 = ((float_int*)buf2) + i;
        ptr3 = ((float_int*)buf3) + i;
      
        if(ptr1->val < ptr2->val) { 
          ptr3->val = ptr2->val; ptr3->rank = ptr2->rank; 
        } else { 
          ptr3->val = ptr1->val; ptr3->rank = ptr1->rank; 
        } 
      }  
    } else if(op == MPI_MINLOC) {
      for(i=0; i<count; i++) {
        typedef struct {
          float val;
          int rank;
        } float_int;
        float_int *ptr1, *ptr2, *ptr3;
                            
        ptr1 = ((float_int*)buf1) + i;
        ptr2 = ((float_int*)buf2) + i;
        ptr3 = ((float_int*)buf3) + i;
      
        if(ptr1->val > ptr2->val) { 
          ptr3->val = ptr2->val; ptr3->rank = ptr2->rank; 
        } else { 
          ptr3->val = ptr1->val; ptr3->rank = ptr1->rank; 
        } 
      }  
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_DOUBLE_INT) { 
    if(op == MPI_MAXLOC) {
      for(i=0; i<count; i++) {
        typedef struct {
          double val;
          int rank;
        } double_int;
        double_int *ptr1, *ptr2, *ptr3;
                            
        ptr1 = ((double_int*)buf1) + i;
        ptr2 = ((double_int*)buf2) + i;
        ptr3 = ((double_int*)buf3) + i;
      
        if(ptr1->val < ptr2->val) { 
          ptr3->val = ptr2->val; ptr3->rank = ptr2->rank; 
        } else { 
          ptr3->val = ptr1->val; ptr3->rank = ptr1->rank; 
        } 
      }  
    } else if(op == MPI_MINLOC) {
      for(i=0; i<count; i++) {
        typedef struct {
          double val;
          int rank;
        } double_int;
        double_int *ptr1, *ptr2, *ptr3;
                            
        ptr1 = ((double_int*)buf1) + i;
        ptr2 = ((double_int*)buf2) + i;
        ptr3 = ((double_int*)buf3) + i;
      
        if(ptr1->val > ptr2->val) { 
          ptr3->val = ptr2->val; ptr3->rank = ptr2->rank; 
        } else { 
          ptr3->val = ptr1->val; ptr3->rank = ptr1->rank; 
        } 
      }  
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_LONG_INT) { 
    if(op == MPI_MAXLOC) {
      for(i=0; i<count; i++) {
        typedef struct {
          long val;
          int rank;
        } long_int;
        long_int *ptr1, *ptr2, *ptr3;
                            
        ptr1 = ((long_int*)buf1) + i;
        ptr2 = ((long_int*)buf2) + i;
        ptr3 = ((long_int*)buf3) + i;
      
        if(ptr1->val < ptr2->val) { 
          ptr3->val = ptr2->val; ptr3->rank = ptr2->rank; 
        } else { 
          ptr3->val = ptr1->val; ptr3->rank = ptr1->rank; 
        } 
      }  
    } else if(op == MPI_MINLOC) {
      for(i=0; i<count; i++) {
        typedef struct {
          long val;
          int rank;
        } long_int;
        long_int *ptr1, *ptr2, *ptr3;
                            
        ptr1 = ((long_int*)buf1) + i;
        ptr2 = ((long_int*)buf2) + i;
        ptr3 = ((long_int*)buf3) + i;
      
        if(ptr1->val > ptr2->val) { 
          ptr3->val = ptr2->val; ptr3->rank = ptr2->rank; 
        } else { 
          ptr3->val = ptr1->val; ptr3->rank = ptr1->rank; 
        } 
      }  
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_2INT) { 
    if(op == MPI_MAXLOC) {
      for(i=0; i<count; i++) {
        typedef struct {
          int val;
          int rank;
        } int_int;
        int_int *ptr1, *ptr2, *ptr3;
                            
        ptr1 = ((int_int*)buf1) + i;
        ptr2 = ((int_int*)buf2) + i;
        ptr3 = ((int_int*)buf3) + i;
      
        if(ptr1->val < ptr2->val) { 
          ptr3->val = ptr2->val; ptr3->rank = ptr2->rank; 
        } else { 
          ptr3->val = ptr1->val; ptr3->rank = ptr1->rank; 
        } 
      }  
    } else if(op == MPI_MINLOC) {
      for(i=0; i<count; i++) {
        typedef struct {
          int val;
          int rank;
        } int_int;
        int_int *ptr1, *ptr2, *ptr3;
                            
        ptr1 = ((int_int*)buf1) + i;
        ptr2 = ((int_int*)buf2) + i;
        ptr3 = ((int_int*)buf3) + i;
      
        if(ptr1->val > ptr2->val) { 
          ptr3->val = ptr2->val; ptr3->rank = ptr2->rank; 
        } else { 
          ptr3->val = ptr1->val; ptr3->rank = ptr1->rank; 
        } 
      }  
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_SHORT_INT) { 
    if(op == MPI_MAXLOC) {
      for(i=0; i<count; i++) {
        typedef struct {
          short val;
          int rank;
        } short_int;
        short_int *ptr1, *ptr2, *ptr3;
                            
        ptr1 = ((short_int*)buf1) + i;
        ptr2 = ((short_int*)buf2) + i;
        ptr3 = ((short_int*)buf3) + i;
      
        if(ptr1->val < ptr2->val) { 
          ptr3->val = ptr2->val; ptr3->rank = ptr2->rank; 
        } else { 
          ptr3->val = ptr1->val; ptr3->rank = ptr1->rank; 
        } 
      }  
    } else if(op == MPI_MINLOC) {
      for(i=0; i<count; i++) {
        typedef struct {
          short val;
          int rank;
        } short_int;
        short_int *ptr1, *ptr2, *ptr3;
                            
        ptr1 = ((short_int*)buf1) + i;
        ptr2 = ((short_int*)buf2) + i;
        ptr3 = ((short_int*)buf3) + i;
      
        if(ptr1->val > ptr2->val) { 
          ptr3->val = ptr2->val; ptr3->rank = ptr2->rank; 
        } else { 
          ptr3->val = ptr1->val; ptr3->rank = ptr1->rank; 
        } 
      }  
    } else return NBC_OP_NOT_SUPPORTED; 
  } else if(type == MPI_LONG_DOUBLE_INT) { 
    if(op == MPI_MAXLOC) {
      for(i=0; i<count; i++) {
        typedef struct {
          long double val;
          int rank;
        } long_double_int;
        long_double_int *ptr1, *ptr2, *ptr3;
                            
        ptr1 = ((long_double_int*)buf1) + i;
        ptr2 = ((long_double_int*)buf2) + i;
        ptr3 = ((long_double_int*)buf3) + i;
      
        if(ptr1->val < ptr2->val) { 
          ptr3->val = ptr2->val; ptr3->rank = ptr2->rank; 
        } else { 
          ptr3->val = ptr1->val; ptr3->rank = ptr1->rank; 
        } 
      }  
    } else if(op == MPI_MINLOC) {
      for(i=0; i<count; i++) {
        typedef struct {
          long double val;
          int rank;
        } long_double_int;
        long_double_int *ptr1, *ptr2, *ptr3;
                            
        ptr1 = ((long_double_int*)buf1) + i;
        ptr2 = ((long_double_int*)buf2) + i;
        ptr3 = ((long_double_int*)buf3) + i;
      
        if(ptr1->val > ptr2->val) { 
          ptr3->val = ptr2->val; ptr3->rank = ptr2->rank; 
        } else { 
          ptr3->val = ptr1->val; ptr3->rank = ptr1->rank; 
        } 
      }  
    } else return NBC_OP_NOT_SUPPORTED; 
  } else return NBC_DATATYPE_NOT_SUPPORTED;
  
  return NBC_OK;
}

