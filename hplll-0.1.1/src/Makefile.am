AUTOMAKE_OPTIONS = gnu

include_hpllldir=$(includedir)/hplll

include_hplll_HEADERS=  hplll.h defs.h  mat.h matpe.h  ldpe.h nr-ld.cpp  \
	hlll.h hlll.cc relations.h relations.cc  \
	matgen.h matgen.cc  \
  	mixed-col.h  matmixed.h  tools.h  \
	givtimer.h givtimer.C ratio.h nr_Z_l128.h wrappers.h

#include_hplll_HEADERS=  hplll.h defs.h  mat.h matpe.h  ldpe.h nr-ld.cpp  \
	#hlll.h hlll.cc decomp.h decomp.cc relations.h relations.cc  \
	#matgen.h matgen.cc decompz.h decompz.cc nullspace.h  nullspace.cc \
  	#mixed-col.h  matmixed.h  lehmer.cc  l1.h  block.h  tools.h  \
	#givtimer.h givtimer.C plll.h plll.cc ratio.h nr_Z_l128.h

#EXTRA_DIST=README.html

#lib_LTLIBRARIES=libhplll.la

#libhplll_la_SOURCES=  
#libfplll_la_CPPFLAGS=$(AM_CPPFLAGS)
#libhplll_la_LIBADD= -lgmp -lmpfr -lfplll

#libhplll_la_LDFLAGS=-no-undefined -version-info 1:0:1

bin_PROGRAMS=isreduced hplll

isreduced_SOURCES=isreduced.cc
hplll_SOURCES=hplll.cc

install-data-hook:
	echo "#include <hplll/hplll.h>" > "$(srcdir)/hplll.h.root"
	echo "using namespace hplll;" >> "$(srcdir)/hplll.h.root"
	$(INSTALL) -m 644 "$(srcdir)/hplll.h.root" "${DESTDIR}$(includedir)/hplll.h"
	rm -f "$(srcdir)/hplll.h.root"

uninstall-hook:
	(test -d "${DESTDIR}$(includedir)/hplll" && rmdir "${DESTDIR}$(includedir)/hplll") || true
	rm -f "${DESTDIR}$(includedir)/hplll.h"


