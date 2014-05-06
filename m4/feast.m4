dnl Copyright (C) 2014, Ask Hjorth Larsen
dnl Shameless cargo-cult bastardization of adjacent .m4 files

AC_DEFUN([ACX_FEAST], [

dnl FEAST works without MPI, but let's not go there.
AC_REQUIRE([ACX_MPI])
AC_REQUIRE([ACX_BLAS])

acx_feast_ok=no

dnl Check if the library was given in the command line
if test $acx_feast_ok = no; then
dnl how the h would acx_feast be anything other than 'no' at this point?
  AC_ARG_WITH(feast, [AS_HELP_STRING([--with-feast=<lib>], [use FEAST eigensolver library <lib>. <lib> should be the path to the parallel library typically libpfeast.a])])
  case $with_feast in
    yes | "") ;;
    no) acx_feast_ok=disable ;;
    -* | */* | *.a | *.so | *.so.* | *.o) LIBS_FEAST="$with_feast" ;;
    *) LIBS_FEAST="-l$with_feast" ;;
  esac
fi

dnl Backup LIBS and (or not) FCFLAGS.  Apparently not
acx_feast_save_LIBS="$LIBS"
dnl acx_feast_save_FCFLAGS="$FCFLAGS"

dnl The tests
AC_MSG_CHECKING([for FEAST])
if test "$acx_feast_ok" != disabled; then
  LIBS="$LIBS_FEAST $LIBS_BLAS $acx_feast_save_LIBS"
  AC_LINK_IFELSE(AC_LANG_PROGRAM([],[
    call zfeast_grcix()
    ]), [acx_feast_ok=yes], [])
fi
AC_MSG_RESULT([$acx_feast_ok ($LIBS_FEAST)])


dnl if test $acx_feast_ok = no; then
dnl   AC_SACRIFICE([GOAT])
dnl fi

LIBS="$acx_feast_save_LIBS"

if test x"$acx_feast_ok" = xyes; then
  AC_DEFINE(HAVE_FEAST,1,[Defined if you have the FEAST library.])
  $1
else
  AC_MSG_WARN([Could not find FEAST library.
      	      *** Will compile without FEAST support])
  $2
fi

])
