## Copyright (C) 2011 J. Alberdi
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
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id: pfft.m4 6722 2010-06-13 12:44:43Z joseba $
##

AC_DEFUN([ACX_PFFT], [
AC_REQUIRE([ACX_FFT])
acx_pfft_ok=no

dnl Check if the library was given in the command line
AC_ARG_WITH(pfft-prefix, [AS_HELP_STRING([--with-pfft-prefix=<lib>], [http://www-user.tu-chemnitz.de/~mpip/software.php])])
case $with_pfft_prefix in
  yes | "") ;;
  no) acx_pfft_ok=disable ;;
  *.a | *.so | *.so.* | *.o) LIBS_PFFT="-L$with_pfft_prefix" ;;
  *) LIBS_PFFT="-L$with_pfft_prefix/lib"; 
     FCFLAGS_PFFT="$ax_cv_f90_modflag$with_pfft_prefix/include" ;;
esac

dnl The include dir must be specified when the library is given with a 
dnl specified file to be compiled static (i.e. *.a etc.)
AC_ARG_WITH(pfft-include, [AS_HELP_STRING([--with-pfft-include=DIR], [PFFT Fortran include files directory.])])
case $with_pfft_include in
  "") if test "x$FCFLAGS_PFFT" == x; then
  FCFLAGS_PFFT="$ax_cv_f90_modflag/usr/include"
  fi;;
  *)  FCFLAGS_PFFT="$ax_cv_f90_modflag$with_pfft_include" ;;
esac


dnl We need to link against the MPI FFTW3 used to compile PFFT 
dnl This test could in principle live on it's own m4 file but since it's just functional
dnl to PFFT for the moment I will leave it here.
AC_ARG_WITH(mpifftw-prefix, [AS_HELP_STRING([--with-mpifftw-prefix=DIR], [MPI FFTW3 libraries directory (the one used to build PFFT).])])
case $with_mpifftw_prefix in
 "") LIBS_MPIFFT="-L/usr/lib -lfftw3_mpi -lfftw3";
     FCFLAGS_MPIFFT="$ax_cv_f90_modflag /usr/include";;
 *)  LIBS_MPIFFT="-L$with_mpifftw_prefix/lib -lfftw3_mpi -lfftw3";
     FCFLAGS_MPIFFT="$ax_cv_f90_modflag$with_mpifftw_prefix/include" ;;
esac
FCFLAGS_PFFT="$FCFLAGS_PFFT $FCFLAGS_MPIFFT"


dnl We cannot use PFFT if MPI is not found
if test "x$acx_mpi_ok" != xyes; then
  acx_pfft_ok=nompi
fi

dnl We cannot use PFFT if FFTW3 is not found
if test "x$acx_fft_ok" != xyes; then
  acx_pfft_ok=nofftw3
fi

dnl Backup LIBS and FCFLAGS
acx_pfft_save_LIBS="$LIBS"
acx_pfft_save_FCFLAGS="$FCFLAGS"


testprogram="AC_LANG_PROGRAM([],[ 
    include 'fftw3.f'
    include 'pfft.f'
    integer :: x = PFFT_REDFT00
    call dpfft_plan_dft_3d()
  ])"


dnl First, check LIBS_PFFT environment variable
if test x"$acx_pfft_ok" = xno; then
  LIBS="$LIBS_PFFT $LIBS_MPIFFT $acx_pfft_save_LIB"
  FCFLAGS="$FCFLAGS_PFFT $acx_pfft_save_FCFLAGS"
  AC_MSG_CHECKING([for pfft library])
  AC_LINK_IFELSE($testprogram, [acx_pfft_ok=yes; LIBS_PFFT="$LIBS_PFFT $LIBS_MPIFFT"], [])
  if test $acx_pfft_ok = no; then
    AC_MSG_RESULT([$acx_pfft_ok])
  else
    AC_MSG_RESULT([$acx_pfft_ok ($LIBS_PFFT)])
  fi
fi


dnl Generic PFFT library 
if test $acx_pfft_ok = no; then
  AC_MSG_CHECKING([for pfft library with -lpfft])
  if test "$LIBS_PFFT" = ""; then
    LIBS="-lpfft $LIBS_MPIFFT $LIBS $acx_pfft_save_LIB"
    FCFLAGS="$FCFLAGS_FFT $acx_pfft_save_FCFLAGS"
    AC_LINK_IFELSE($testprogram, [acx_pfft_ok=yes; LIBS_PFFT="-lpfft $LIBS_MPIFFT"], [])
  else
    LIBS="$LIBS_PFFT -lpfft $LIBS_MPIFFT $acx_pfft_save_LIB"
    FCFLAGS="$FCFLAGS_PFFT $FCFLAGS_FFT $acx_pfft_save_FCFLAGS"    
    AC_LINK_IFELSE($testprogram, [acx_pfft_ok=yes; 
                                  LIBS_PFFT="$LIBS_PFFT -lpfft $LIBS_MPIFFT"], [])  
  fi
  if test $acx_pfft_ok = no; then
    AC_MSG_RESULT([$acx_pfft_ok])
  else
    AC_MSG_RESULT([$acx_pfft_ok ($LIBS_PFFT)])
  fi
fi

dnl Usually fft library are present in the system.
dnl In order to have the correct symbols for PFFT we need to override the default.
if test $acx_pfft_ok = yes; then
  LDFLAGS="$LDFLAGS $LIBS_MPIFFT"
fi

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_pfft_ok" = xyes; then
  AC_DEFINE(HAVE_PFFT,1,[Defined if you have PFFT library.])
  $1
else
  AC_MSG_WARN([Could not find PFFT library. 
               *** Will compile without PFFT support])
  LIBS_PFFT=""
  FCFLAGS_PFFT=""
  $2
fi

AC_SUBST(LIBS_PFFT)
AC_SUBST(FCFLAGS_PFFT)
LIBS="$acx_pfft_save_LIBS"
FCFLAGS="$acx_pfft_save_FCFLAGS"

])dnl ACX_PFFT
