#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_CONFIG_MACRO_DIR([macros])
AC_INIT([hplll],[1.0.0])
AC_CONFIG_SRCDIR([src/hlll.cc])
AM_INIT_AUTOMAKE

#Recompile pour toute modif dans la chaine autotool
AM_MAINTAINER_MODE(enable)

AC_CONFIG_HEADER(src/config.h)

DEFAULT_CFLAGS="-O3 -Wall"


CXXFLAGS=${CXXFLAGS:-$DEFAULT_CFLAGS}

LT_INIT 

AC_PROG_INSTALL

# Checks for programs.
#AC_LANG([C++])
AC_LANG_CPLUSPLUS
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AX_CXX_COMPILE_STDCXX([11],[noext],[mandatory])


AC_ARG_ENABLE([quadmath],
    [AC_HELP_STRING([--enable-quadmath],
        [enable support for __float128 floating point arithmetics (default test: yes)])],
    [enable_quadmath="${enableval}"],
    [enable_quadmath="yes"]
)

if test "$enable_quadmath" = "yes"; then
  LIBQUAD_LIBADD="-lquadmath"
        AC_CHECK_LIB(quadmath,sinhq)
  AC_SUBST(LIBQUAD_LIBADD)
fi

#AC_ARG_ENABLE(quadmath, 
   #[ --enable-quadmath  Enable quadmath], [
     #LDFLAGS="${LDFLAGS} -lquadmath"
     #LIBS="${BACKUP_LIBS} -lquadmath"
   #])


AC_ARG_WITH(mpfr,
   [  --with-mpfr=DIR          MPFR install directory ], [
     if test "$withval" != "" ; then
      CXXFLAGS="$CXXFLAGS -I$withval/include"
       LDFLAGS="$LDFLAGS -L$withval/lib -Wl,-rpath,$withval/lib"
       mpfr_lib_path="$withval/lib"
      fi
   ])

AC_ARG_WITH(gmp,
   [  --with-gmp=DIR          GMP install directory ], [
      if test "$withval" != "" ; then
      CXXFLAGS="$CXXFLAGS -I$withval/include"
       LDFLAGS="$LDFLAGS -L$withval/lib -Wl,-rpath,$withval/lib"
       gmp_lib_path="$withval/lib"
       fi
   ])


# Checks for libraries.
AC_CHECK_LIB(gmp, __gmpz_init, , [AC_MSG_ERROR(
[GNU MP not found, see http://gmplib.org])])

AC_CHECK_LIB(mpfr, mpfr_add, , [AC_MSG_ERROR(
[MPFR not found, see http://www.mpfr.org])])


AC_ARG_WITH(qd, AS_HELP_STRING([--with-qd=@<:@=DIR@:>@], [libqd install directory]),)

AS_IF([test "x$with_qd" != "xno"], [
       if test "$with_qd" = "" ;  then 
             CPPFLAGS="$CPPFLAGS"
             LDFLAGS="$LDFLAGS"
       else 
             CPPFLAGS="$CPPFLAGS -I$with_qd/include"
             LDFLAGS="$LDFLAGS -L$with_qd/lib"
       fi 
       qd_lib_path="$with_qd/lib"
       AC_CHECK_HEADER(qd/dd_real.h,
           have_libqd="yes",
           have_libqd="no")
      ])

AS_IF([test "x${have_libqd}" = "xyes"], [
    LIBQD_LIBADD="-lqd"
    AC_DEFINE([HPLLL_WITH_QD], [1], [defined when libqd is available])])

AC_SUBST(LIBQD_LIBADD)


AC_ARG_ENABLE(omp, 
   [ --enable-omp Enables OpenMP], [
     CXXFLAGS="$CXXFLAGS  -fopenmp"
   ])


# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h cstdio iostream string limits vector])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_C_VOLATILE

# Checks for library functions.
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow rint sqrt strtol utime])


AC_CHECK_LIB(gmp, __gmp_get_memory_functions, , [AC_MSG_ERROR(
[GMP version too old, need >= 4.2.0, see http://gmplib.org])])

AC_CHECK_LIB(mpfr, mpfr_fms, , [AC_MSG_ERROR(
[MPFR version too old, need >= 2.3.0, see http://www.mpfr.org])])


LB_CHECK_FPLLL 

LB_CHECK_NTL

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile macros/Makefile examples/Makefile drafts/Makefile \
                 fpbench/Makefile])
AC_OUTPUT

