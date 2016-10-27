#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_CONFIG_MACRO_DIR([macros])
AC_INIT([hplll],[0.1.2])
AC_CONFIG_SRCDIR([src/hlll.cc])
AM_INIT_AUTOMAKE

#Recompile pour toute modif dans la chaine autotool
AM_MAINTAINER_MODE(enable)

AC_CONFIG_HEADER(config.h)

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

AC_ARG_ENABLE(omp, 
   [ --enable-omp Enables OpenMP], [
     CXXFLAGS="$CXXFLAGS  -fopenmp"
     #CXXFLAGS="$CXXFLAGS  -fopenmp"
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


AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile macros/Makefile examples/Makefile drafts/Makefile])
AC_OUTPUT

