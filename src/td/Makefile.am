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

noinst_LIBRARIES = libtd.a

libtd_a_SOURCES = \
	pes.F90 spectrum.F90 propagator.F90 td_calc.F90 td_write.F90 td.F90 \
	ob_propagator.F90 ob_src.F90 ob_mem.F90 ob_terms.F90 tdpsf.F90 \
	pes_rc.F90                 \
	pes_mask.F90

noinst_HEADERS = \
	pes_mask_out_inc.F90         \
	td_init_inc.F90              \
	propagator_qoct_inc.F90      \
	propagator_rk_inc.F90        \
	propagator_vksold_inc.F90    \
	propagator_cn_inc.F90        \
	propagator_magnus_inc.F90    \
	propagator_expmid_inc.F90    \
	propagator_etrs_inc.F90

include $(srcdir)/../common-rules.make

depend:
	@makedepf90 $(libtd_a_SOURCES)

pes.o: pes.F90 pes_mask.o pes_rc.o
pes_rc.o : pes_rc.F90
pes_mask.o : pes_mask.F90 pes_mask_out_inc.F90 tdpsf.o
tdpsf.o: tdpsf.F90
spectrum.o: spectrum.F90
propagator.o: propagator.F90 ob_propagator.o propagator_qoct_inc.F90 propagator_rk_inc.F90 propagator_vksold_inc.F90 propagator_cn_inc.F90 propagator_magnus_inc.F90 propagator_expmid_inc.F90 propagator_etrs_inc.F90
ob_propagator.o: ob_propagator.F90 ob_mem.o ob_src.o ob_terms.o
ob_mem.o: ob_mem.F90 ob_terms.o
ob_src.o: ob_src.F90 ob_mem.o ob_terms.o
td_write.o: td_write.F90 spectrum.o td_calc.o
td.o: td.F90 td_init_inc.F90 spectrum.o pes.o td_write.o propagator.o

