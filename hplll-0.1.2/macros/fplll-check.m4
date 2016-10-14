# Check for fplll  
# Modified G. Villard from LinBox Fri Feb 13 12:22:54
# Copyright (c) the LinBox group


# Test for Damien Stehlé's FPLLL (Lattice reduction) 

AC_DEFUN([LB_CHECK_FPLLL],
[

AC_ARG_WITH(fplll,
[AC_HELP_STRING([--with-fplll=<path>|yes], [The FPLLL library is mandatory. 
                         If argument is yes or <empty> that means
                         the library is reachable with the standard search path
                         "/usr" or "/usr/local" (set as default). Otherwise you
                         give the <path> to the directory which contain the
                         library.
	     ])],
	     [if test "$withval" = yes ; then
			FPLLL_HOME_PATH="/usr /usr/local"
	      elif test "$withval" != no ; then
			FPLLL_HOME_PATH="$withval /usr /usr/local"
	     fi],
	     [FPLLL_HOME_PATH="/usr /usr/local"])

# TO SEE FOR VERSION testing if needed (see e.g. test NTL) 

min_fplll_version=ifelse([$1], ,4.0,$1)

dnl Check for existence
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LDFLAGS}

if test -n "$FPLLL_HOME_PATH"; then
AC_MSG_CHECKING(for fplll)
fi

for FPLLL_HOME in ${FPLLL_HOME_PATH}
 do
if test -r "$FPLLL_HOME/include/fplll/defs.h"; then

	if test "x$FPLLL_HOME" != "x/usr" -a "x$FPLLL_HOME" != "x/usr/local"; then
		FPLLL_CFLAGS="-I${FPLLL_HOME}/include"
		FPLLL_LIBS="-L${FPLLL_HOME}/lib -lfplll"
	else
		FPLLL_CFLAGS="-I/usr/local/include"
		FPLLL_LIBS="-L/usr/local/lib -lfplll"
	fi
	CXXFLAGS="${BACKUP_CXXFLAGS} ${FPLLL_CFLAGS}"
	LDFLAGS="${BACKUP_LIBS} ${FPLLL_LIBS} -Wl,-rpath,${FPLLL_HOME}/lib"

	AC_TRY_LINK(
	[#include <fplll.h>],
	[FP_NR<double> x;],
	[
	AC_TRY_RUN(
	[#include <fplll.h>
	int main () { FP_NR<double> x; x=1.0; return 0;}
	],[
	fplll_found="yes"
	break
	],[
	fplll_problem="$problem $FPLLL_HOME"
	unset FPLLL_CFLAGS
	unset FPLLL_LIBS
	],[
	fplll_found="yes"
	fplll_cross="yes"
	break
	])
	],
	[
	fplll_found="no"
	fplll_checked="$checked $FPLLL_HOME"
	unset FPLLL_CFLAGS
	unset FPLLL_LIBS
	])
else
	fplll_found="no"
fi
done

if test "x$fplll_found" = "xyes" ; then

	unset FPLLL_CFLAGS
	unset FPLLL_LIBS

        LIBS="${LIBS} -lfplll"
        AC_SUBST(LIBS)

	AC_SUBST(LDFLAGS)

	AC_SUBST(CXXFLAGS)

	AC_DEFINE(HAVE_FPLLL,1,[Define if FPLLL is installed])
	HAVE_FPLLL=yes
	if test "x$fplll_cross" != "xyes"; then
		AC_MSG_RESULT(yes)
	else
		AC_MSG_RESULT(unknown)
		echo "WARNING: You appear to be cross compiling, so there is no way to determine"
		echo "whether your FPLLL version is new enough. I am assuming it is."
	fi
	ifelse([$2], , :, [$2])
elif test -n "$fplll_problem"; then

        CXXFLAGS=${BACKUP_CXXFLAGS}
        LDFLAGS=${BACKUP_LIBS}

	AC_MSG_RESULT(problem)
	#echo "Sorry, your FPLLL version is too old. Disabling."
	AC_MSG_ERROR(Problem with fplll. fplll is mandatory for compiling hplll.)
	ifelse([$3], , :, [$3])
elif test   "x$fplll_found" = "xno";  then

        CXXFLAGS=${BACKUP_CXXFLAGS}
        LDFLAGS=${BACKUP_LIBS}

	AC_MSG_RESULT(not found)

	# GV 
	AC_MSG_ERROR(fplll was not found. fplll is mandatory for compiling hplll.)

	if test "x$FPLLL_HOME" != "x/usr" -a "x$FPLLL_HOME" != "x/usr/local" ; then
	AC_MSG_ERROR(fplll was not found. fplll is mandatory for compiling hplll.)
	fi
	ifelse([$3], , :, [$3])
fi



AM_CONDITIONAL(HPLLL_HAVE_FPLLL, test "x$HAVE_FPLLL" = "xyes")


])

