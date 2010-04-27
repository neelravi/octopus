/*
 Copyright (C) 2010 X. Andrade, N. Suberviola

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

 $Id: opencl.h 2146 2006-05-23 17:36:00Z xavier $
*/

#ifndef OCTOPUS_OPENCL_H
#define OCTOPUS_OPENCL_H

#include <config.h>
#include <CL/cl.h>

typedef struct{
  cl_context Context;
  cl_device_id * Devices;
  cl_command_queue CommandQueue;
  char * source_path;
  size_t max_workgroup_size;
} opencl_env_t;

typedef struct{
  cl_kernel kernel_double;
  cl_kernel kernel_complex;
} opencl_kernel_t;

#endif
