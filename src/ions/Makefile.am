## Process this file with automake to produce Makefile.in

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id$

noinst_LIBRARIES = libions.a

libions_a_SOURCES =         \
	read_coords.F90     \
	atom.F90            \
	geometry.F90        \
	xyz_adjust.F90      \
	kpoints.F90         \
	box.F90             \
	box_union.F90       \
	simul_box.F90       \
	periodic_copy.F90   \
	vibrations.F90      \
	ion_dynamics.F90    \
	symmetries.F90      \
	symmetries_finite.c \
        spglib_f.c          \
        metric.F90

noinst_HEADERS =

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libions_a_SOURCES)

read_coords.o : read_coords.F90 
atom.o : atom.F90
geometry.o : geometry.F90 atom.o read_coords.o
xyz_adjust.o : xyz_adjust.F90 geometry.o 
kpoints.o : geometry.o symmetries.o
box.o : box.F90
box_union.o : box_union.F90 box.o
simul_box.o : simul_box.F90 geometry.o kpoints.o symmetries.o atom.o
periodic_copy.o: periodic_copy.F90 simul_box.o
vibrations.o : vibrations.F90 simul_box.o geometry.o
ion_dynamics.o : ion_dynamics.F90 simul_box.o geometry.o read_coords.o
symmetries.o : symmetries.F90 geometry.o
symmetries_finite.o : symmetries_finite.c
spglib_f.o : spglib_f.c
metric.o : geometry.o
