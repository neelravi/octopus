## Process this file with automake to produce Makefile.in

## Copyright (C) 2002-2014 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

noinst_LIBRARIES = libscf.a

libscf_a_SOURCES =       \
	mix.F90          \
	lcao.F90         \
        rdmft.F90        \
        scf.F90          \
	ground_state.F90 \
	ob_lippmann_schwinger.F90

noinst_HEADERS =     \
	mix_inc.F90  \
	lcao_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libscf_a_SOURCES)

mix.o : mix.F90 mix_inc.F90
lcao.o: lcao.F90 lcao_inc.F90
rdmft.o: rdmft.F90 
scf.o: scf.F90 mix.o lcao.o ob_lippmann_schwinger.o
ground_state.o: ground_state.F90 scf.o rdmft.o lcao.o 
ob_lippmann_schwinger.o: ob_lippmann_schwinger.F90

