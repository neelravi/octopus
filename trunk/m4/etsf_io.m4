AC_DEFUN([ACX_ETSF_IO], [
AC_REQUIRE([ACX_NETCDF])
acx_etsf_io_ok=no
FCFLAGS_ETSF_IO=""
LIBS_ETSF_IO=""

dnl Check if the library was given in the command line
AC_ARG_WITH(etsf-io-prefix, [AS_HELP_STRING([--with-etsf-io-prefix=DIR], [Directory where etsf_io was installed.])])
case $with_etsf_io_prefix in
  no ) acx_etsf_io_ok=disabled ;;
  "") LIBS_ETSF_IO="-letsf_io_utils -letsf_io"; FCFLAGS_ETSF_IO="-I/usr/include" ;;
  *) LIBS_ETSF_IO="-L$with_etsf_io_prefix/lib -letsf_io_utils -letsf_io"; 
     FCFLAGS_ETSF_IO="$ax_cv_f90_modflag$with_etsf_io_prefix/include" ;;
esac

dnl We cannot use etsf_io if netcdf is not found
if test "x$acx_netcdf_ok" != xyes; then
  acx_etsf_io_ok=disabled
  FCFLAGS_ETSF_IO=""
  LIBS_ETSF_IO=""
fi

dnl Backup LIBS and FCFLAGS
acx_etsf_io_save_LIBS="$LIBS"
acx_etsf_io_save_FCFLAGS="$FCFLAGS"

dnl The tests
AC_MSG_CHECKING([for etsf_io])
if test "$acx_etsf_io_ok" != disabled; then
  etsf_io_fcflags="$FCFLAGS_ETSF_IO"; etsf_io_libs="$LIBS_ETSF_IO"
  FCFLAGS="$etsf_io_fcflags $acx_etsf_io_save_FCFLAGS"
  LIBS="$etsf_io_libs $acx_etsf_io_save_LIBS $LIBS_NETCDF"
  AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
    use etsf_io
    type(etsf_vars) :: vars
    call etsf_io_vars_free(vars)
  ]), [acx_etsf_io_ok=yes; FCFLAGS_ETSF_IO="$etsf_io_fcflags"; LIBS_ETSF_IO="$etsf_io_libs"], [])
fi
AC_MSG_RESULT([$acx_etsf_io_ok ($LIBS_ETSF_IO)])

AC_SUBST(FCFLAGS_ETSF_IO)
AC_SUBST(LIBS_ETSF_IO)
FCFLAGS="$acx_etsf_io_save_FCFLAGS"
LIBS="$acx_etsf_io_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_etsf_io_ok" = xyes; then
  AC_DEFINE(HAVE_ETSF_IO,1,[Defined if you have the ETSF_IO library.])
  $1
else
  AC_MSG_WARN([Could not find etsf_io library. 
           *** Will compile without etsf_io support])
  $2
fi

])dnl ACX_ETSF_IO
