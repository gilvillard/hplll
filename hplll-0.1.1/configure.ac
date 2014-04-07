#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_CONFIG_MACRO_DIR([m4])
AC_INIT([hplll],[0.1.1])
AC_CONFIG_SRCDIR([src/hlll.cc])
AM_INIT_AUTOMAKE

#Recompile pour toute modif dans la chaine autotool
AM_MAINTAINER_MODE(enable)

AC_CONFIG_HEADER(config.h)

DEFAULT_CFLAGS=" -O2 -Wall"

CFLAGS=${CFLAGS:-$DEFAULT_CFLAGS}
CXXFLAGS=${CXXFLAGS:-$DEFAULT_CFLAGS}

#AC_DISABLE_SHARED
LT_INIT

AC_PROG_INSTALL

# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET

AC_ARG_WITH(fplll,
   [  --with-fplll=DIR          FPLLL install directory ], [
     CXXFLAGS="$CXXFLAGS -I$withval/include"
     LDFLAGS="$LDFLAGS -L$withval/lib"
     fplll_lib_path="$withval/lib"
   ])

AC_ARG_WITH(mpfr,
   [  --with-mpfr=DIR          MPFR install directory ], [
      CXXFLAGS="$CXXFLAGS -I$withval/include"
       LDFLAGS="$LDFLAGS -L$withval/lib"
       mpfr_lib_path="$withval/lib"
   ])

AC_ARG_WITH(gmp,
   [  --with-gmp=DIR          GMP install directory ], [
      CXXFLAGS="$CXXFLAGS -I$withval/include"
       LDFLAGS="$LDFLAGS -L$withval/lib"
       gmp_lib_path="$withval/lib"
   ])


# Checks for libraries.
AC_CHECK_LIB(gmp, __gmpz_init, , [AC_MSG_ERROR(
[GNU MP not found, see http://gmplib.org])])
AC_CHECK_LIB(mpfr, mpfr_add, , [AC_MSG_ERROR(
[MPFR not found, see http://www.mpfr.org])])



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

AC_LANG(C++)
SAVED_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS -lfplll" 
#echo "checking for cputime in -lfplll... "  # -n or \c 
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <fplll.h>]], [[cputime()]])],[AC_MSG_RESULT(checking for cputime in -lfplll... yes)],[AC_MSG_ERROR([FPLLL not found, see http://perso.ens-lyon.fr/damien.stehle/fplll])
   ])
LDFLAGS=$SAVED_LDFLAGS
# Ne va pas car enlève le lfplll


AC_LANG(C++)
SAVED_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS -lntl" 
#echo "checking for to_long in -lntl... "  # -n or \c 
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <NTL/tools.h>]], [[NTL::to_long(2)]])],[AC_MSG_RESULT(checking for to_long in -lntl... yes)],[AC_MSG_ERROR([NTL not found, see http://www.shoup.net/ntl])
   ])
LDFLAGS=$SAVED_LDFLAGS
# Ne va pas car enlève le lntl


AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT

