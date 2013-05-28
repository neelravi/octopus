AC_DEFUN([ACX_LIBXC], [
acx_libxc_ok=no

dnl Check if the library was given in the command line
dnl if not, use environment variables or defaults
AC_ARG_WITH(libxc-prefix, [AS_HELP_STRING([--with-libxc-prefix=DIR], [Directory where libxc was installed.])])
case $with_libxc_prefix in
  "") if test x"$LIBS_LIBXC" = x; then
        LIBS_LIBXC="-lxcf90 -lxc"
      fi
      if test x"$FCFLAGS_LIBXC" = x; then
        FCFLAGS_LIBXC="-I/usr/include";
      fi ;;
  *) LIBS_LIBXC="$with_libxc_prefix/lib/libxcf90.a $with_libxc_prefix/lib/libxc.a"; FCFLAGS_LIBXC="$ax_cv_f90_modflag$with_libxc_prefix/include" ;;
esac

AC_ARG_WITH(libxc-include, [AS_HELP_STRING([--with-libxc-include=DIR], [Directory where libxc Fortran headers were installed.])])
case $with_libxc_include in
  "") ;;
  *)  FCFLAGS_LIBXC="$ax_cv_f90_modflag$with_libxc_include" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_libxc_save_LIBS="$LIBS"
acx_libxc_save_FCFLAGS="$FCFLAGS"

dnl The tests
AC_MSG_CHECKING([for libxc])

testprog="AC_LANG_PROGRAM([],[
  use xc_f90_lib_m
  implicit none
  integer :: major
  integer :: minor
  call xc_f90_version(major, minor)])"

# first try linking statically
libxc_fcflags="$FCFLAGS_LIBXC"; libxc_libs="$LIBS_LIBXC"
FCFLAGS="$libxc_fcflags $acx_libxc_save_FCFLAGS"
LIBS="$libxc_libs $acx_libxc_save_LIBS"
AC_LINK_IFELSE($testprog, [acx_libxc_ok=yes; FCFLAGS_LIBXC="$libxc_fcflags"; LIBS_LIBXC="$libxc_libs"], [])
AC_MSG_RESULT([$acx_libxc_ok ($FCFLAGS_LIBXC $LIBS_LIBXC)])

# otherwise try linking dynamically (for when libxc was installed as a package)
if test x"$acx_libxc_ok" = xno; then
  LIBS_LIBXC="-L$with_libxc_prefix/lib -lxcf90 -lxc"
  libxc_libs="$LIBS_LIBXC"
  LIBS="$libxc_libs $acx_libxc_save_LIBS"
  AC_LINK_IFELSE($testprog, [acx_libxc_ok=yes; FCFLAGS_LIBXC="$libxc_fcflags"; LIBS_LIBXC="$libxc_libs"], [])
  AC_MSG_RESULT([$acx_libxc_ok ($FCFLAGS_LIBXC $LIBS_LIBXC)])
fi

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_libxc_ok" = xyes; then
  AC_DEFINE(HAVE_LIBXC,1,[Defined if you have the LIBXC library.])
  $1
else
  AC_MSG_ERROR([Could not find required libxc library ( >= v 2.0.0).])
  FCFLAGS_LIBXC=""
  LIBS_LIBXC=""
  $2
fi

AC_SUBST(FCFLAGS_LIBXC)
AC_SUBST(LIBS_LIBXC)
FCFLAGS="$acx_libxc_save_FCFLAGS"
LIBS="$acx_libxc_save_LIBS"
])dnl ACX_LIBXC
