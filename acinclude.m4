dnl
dnl Check size of a pointer
dnl
AC_DEFUN(ACX_POINTER_SIZE,
[AC_MSG_CHECKING([for the size of a pointer])
AC_REQUIRE([AC_PROG_CC])
if test -z "$POINTER_SIZE"; then
cat >pointertest.c <<EOF
#include <stdio.h>
void main()
{
  printf("%ld", sizeof(void *));
}
EOF
	ac_try='$CC $CFLAGS -o pointertest.x pointertest.c 1>&AC_FD_CC'
	if AC_TRY_EVAL(ac_try); then
  	ac_try=""
	else
  	echo "configure: failed program was:" >&AC_FD_CC
	  cat pointertest.c >&AC_FD_CC
  	rm -f pointertest*
	  AC_MSG_ERROR(failed to compile c program to find the size of a pointer)
	fi
	ac_pointersize=`./pointertest.x`;
	rm -f pointertest*
	AC_DEFINE_UNQUOTED(POINTER_SIZE, ${ac_pointersize})
	AC_MSG_RESULT(${ac_pointersize} "bytes")
fi
])

dnl
dnl Check f90 compiler
dnl
AC_DEFUN(ACX_PROG_F90,
[AC_BEFORE([$0], [AC_PROG_CPP])dnl
if test -z "$F90"; then
  AC_CHECK_PROGS(F90, f90 xlf90)
    test -z "$F90" && AC_MSG_ERROR([no acceptable Fortran 90 compiler found in \
$PATH])
fi
])

dnl
dnl detect Fortran name-mangling scheme
dnl adapted from fftw aclocal.m4
dnl
AC_DEFUN(ACX_F90_FUNC_MANGLE,
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([ACX_PROG_F90])
AC_MSG_CHECKING(how f90 mangles function names)
cat > mangle-func.f90 <<EOF
subroutine foobar()
  return
end subroutine foobar
subroutine foo_bar()
  return
end subroutine foo_bar
EOF
ac_try='$F90 -c $F90FLAGS mangle-func.f90 1>&AC_FD_CC'
if AC_TRY_EVAL(ac_try); then
  ac_try=""
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat mangle-func.f >&AC_FD_CC
  rm -f mangle-func*
  AC_MSG_ERROR(failed to compile fortran test program)
fi

ac_f90_mangle_type=unknown
AC_LANG_SAVE
AC_LANG_C
ac_save_LIBS="$LIBS"
LIBS="mangle-func.o $F90LIBS $LIBS"
AC_TRY_LINK(,foobar();,
     ac_f90_mangle_type=lowercase,
     AC_TRY_LINK(,foobar_();,
          ac_f90_mangle_type=lowercase-underscore,
          AC_TRY_LINK(,FOOBAR();,
               ac_f90_mangle_type=uppercase,
               AC_TRY_LINK(,FOOBAR_();,
                    ac_f90_mangle_type=uppercase-underscore))))
LIBS="$ac_save_LIBS"
AC_LANG_RESTORE
AC_MSG_RESULT($ac_f90_mangle_type)

mangle_try=unknown
case $ac_f90_mangle_type in
        lowercase)
                AC_DEFINE(FORTRANIZE_LOWERCASE)
                mangle_try=foo_bar_
                ;;
        lowercase-underscore)
                AC_DEFINE(FORTRANIZE_LOWERCASE_UNDERSCORE)
                mangle_try=foo_bar__
                ;;
        uppercase)
                AC_DEFINE(FORTRANIZE_UPPERCASE)
                mangle_try=FOO_BAR_
                ;;
        uppercase-underscore)
                AC_DEFINE(FORTRANIZE_UPPERCASE_UNDERSCORE)
                mangle_try=FOO_BAR__
                ;;
esac

AC_MSG_CHECKING(if f90 functions with an underscore get an extra underscore)

AC_LANG_SAVE
AC_LANG_C
ac_save_LIBS="$LIBS"
LIBS="mangle-func.o $F90LIBS $LIBS"
AC_TRY_LINK(,$mangle_try();,
            [ac_f90_mangle_underscore=yes;
             AC_DEFINE(FORTRANIZE_EXTRA_UNDERSCORE)],
            [ac_f90_mangle_underscore=no])
LIBS="$ac_save_LIBS"
AC_LANG_RESTORE
rm -f mangle-func*
AC_MSG_RESULT($ac_f90_mangle_underscore)
])

dnl AC_LANG_FORTRAN90()
AC_DEFUN(AC_LANG_FORTRAN90,
[define([AC_LANG], [FORTRAN90])dnl
ac_ext=f90
ac_compile='${F90-f90} -c $F90FLAGS conftest.$ac_ext 1>&AC_FD_CC'
ac_link='${F90-f90} -o conftest${ac_exeext} $F90FLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&AC_FD_CC'
cross_compiling=$ac_cv_prog_f90_cross
])

dnl this would be much simpler if autoconf would
dnl support fortran90
AC_DEFUN(ACX_TRY_LINK_F90,
[
cat > conftest.f90 <<EOF
program main
  call [$1]
end program main
EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
  ifelse([$2], , :, [rm -rf conftest*
  $2])
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat conftest.f90 >&AC_FD_CC
ifelse([$3], , , [  rm -rf conftest*
  $3
])dnl
fi
rm -f conftest*])
])

dnl
dnl try to link a function using the f90 compiler
dnl
AC_DEFUN(ACX_CHECK_FUNC_F90,
[AC_MSG_CHECKING([for $1])
AC_CACHE_VAL(ac_cv_func_$1,
[ACX_TRY_LINK_F90([$1()], eval "ac_cv_func_$1=yes", eval "ac_cv_func_$1=no")
])
if eval "test \"`echo '$ac_cv_func_'$1`\" = yes"; then
  AC_MSG_RESULT(yes)
  ifelse([$2], , :, [$2])
else
  AC_MSG_RESULT(no)
ifelse([$3], , , [$3])dnl
fi
])

AC_DEFUN(ACX_CHECK_LIB_F90,
[AC_MSG_CHECKING([for $2 in -l$1])
ac_lib_var=`echo $1['_']$2 | sed 'y%./+-%__p_%'`
AC_CACHE_VAL(ac_cv_lib_$ac_lib_var,
[ac_save_LIBS="$LIBS"
LIBS="-l$1 $5 $LIBS"
ACX_TRY_LINK_F90([$2()], eval "ac_cv_lib_$ac_lib_var=yes", eval "ac_cv_lib_$ac_lib_var=no")
LIBS="$ac_save_LIBS"
])dnl
if eval "test \"`echo '$ac_cv_lib_'$ac_lib_var`\" = yes"; then
  AC_MSG_RESULT(yes)
  ifelse([$3], ,
[changequote(, )dnl
  ac_tr_lib=HAVE_LIB`echo $1 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($ac_tr_lib)
  LIBS="-l$1 $LIBS"
], [$3])
else
  AC_MSG_RESULT(no)
ifelse([$4], , , [$4
])dnl
fi
])
