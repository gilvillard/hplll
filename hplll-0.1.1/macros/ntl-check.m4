# Check for ntl 
# Modified G. Villard from LinBox Fri Feb 13 09:09:23
# Copyright (c) the LinBox group


# Test for Victor Shoup's NTL (Number Theory Library) and define
# NTL_CFLAGS and NTL_LIBS

AC_DEFUN([LB_CHECK_NTL],
[

AC_ARG_WITH(ntl,
[AC_HELP_STRING([--with-ntl=<path>|yes|no], [Use NTL library. If argument is no, you do not have
                            the library installed on your machine, and no check is 
			    done. Nothing specified : check done by default with 
			    no warning. If argument is yes or <empty> that means
			    the library should be reachable with the standard search
			    path (/usr or /usr/local). Otherwise you give the
			    <path> to the directory which contain the library.
	     ])],
	     [if test "$withval" = yes ; then
			# The "/" bad trick for forcing the warning at the end 
			# to be modified 
			NTL_HOME_PATH="/usr /usr/local /"
	      elif test "$withval" != no ; then
			NTL_HOME_PATH="$withval /usr /usr/local /"
	     fi],
	     [NTL_HOME_PATH="/usr /usr/local"])


# ICI A VOIR VERSION 
min_ntl_version=ifelse([$1], ,5.0,$1)

dnl Check for existence
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LDFLAGS}

if test -n "$NTL_HOME_PATH"; then
AC_MSG_CHECKING(for ntl)
fi

for NTL_HOME in ${NTL_HOME_PATH}
 do
if test -r "$NTL_HOME/include/NTL/ZZ.h"; then

	if test "x$NTL_HOME" != "x/usr" -a "x$NTL_HOME" != "x/usr/local"; then
		NTL_CFLAGS="-I${NTL_HOME}/include"
		NTL_LIBS="-L${NTL_HOME}/lib -lntl"
	else
		NTL_CFLAGS=
		NTL_LIBS="-lntl"
	fi
	#CXXFLAGS="${BACKUP_CXXFLAGS} ${NTL_CFLAGS} ${GMP_CFLAGS}"
	# GV
	CXXFLAGS="${BACKUP_CXXFLAGS} ${NTL_CFLAGS}"
	#LDFLAGS="${BACKUP_LIBS} ${NTL_LIBS} ${GMP_LIBS}"
	# GV
	LDFLAGS="${BACKUP_LIBS} ${NTL_LIBS}"

	AC_TRY_LINK(
	[#include <NTL/ZZ.h>],
	[NTL::ZZ a;],
	[
	AC_TRY_RUN(
	[#include <NTL/version.h>
	int main () { if (NTL_MAJOR_VERSION < 5) return -1; else return 0; }
	],[
	ntl_found="yes"
	break
	],[
	ntl_problem="$problem $NTL_HOME"
	unset NTL_CFLAGS
	unset NTL_LIBS
	],[
	ntl_found="yes"
	ntl_cross="yes"
	break
	])
	],
	[
	ntl_found="no"
	ntl_checked="$checked $NTL_HOME"
	unset NTL_CFLAGS
	unset NTL_LIBS
	])
else
	ntl_found="no"
fi
done

if test "x$ntl_found" = "xyes" ; then
	#AC_SUBST(NTL_CFLAGS)
	#AC_SUBST(NTL_LIBS)
	# GV
	unset NTL_CFLAGS
	# GV
	unset NTL_LIBS

	# GV ??
        LIBS="${LIBS} -lntl"
        AC_SUBST(LIBS)

	# GV ??
	AC_SUBST(LDFLAGS)

	# GV ??
	AC_SUBST(CXXFLAGS)

	AC_DEFINE(HAVE_NTL,1,[Define if NTL is installed])
	HAVE_NTL=yes
	if test "x$ntl_cross" != "xyes"; then
		AC_MSG_RESULT(yes)
	else
		AC_MSG_RESULT(unknown)
		echo "WARNING: You appear to be cross compiling, so there is no way to determine"
		echo "whether your NTL version is new enough. I am assuming it is."
	fi
	ifelse([$2], , :, [$2])
elif test -n "$ntl_problem"; then

        CXXFLAGS=${BACKUP_CXXFLAGS}
        LDFLAGS=${BACKUP_LIBS}

	AC_MSG_RESULT(problem)
	AC_MSG_WARN(ntl >= $min_ntl_version was not found. hplll may need the NTL namespace to be enabled.)
	#echo "Sorry, your NTL version is too old. Disabling."
	ifelse([$3], , :, [$3])
elif test   "x$ntl_found" = "xno";  then

        CXXFLAGS=${BACKUP_CXXFLAGS}
        LDFLAGS=${BACKUP_LIBS}

        # Result by default search, hence no warning 
	AC_MSG_RESULT(no)
        # Result for yes, hence warning 
	if test "x$NTL_HOME" != "x/usr" -a "x$NTL_HOME" != "x/usr/local" ; then
	AC_MSG_WARN(ntl >= $min_ntl_version was not found. hplll may need the NTL namespace to be enabled.)
	fi
	ifelse([$3], , :, [$3])
fi



AM_CONDITIONAL(HPLLL_HAVE_NTL, test "x$HAVE_NTL" = "xyes")

])

