AUTOMAKE_OPTIONS = gnu

include_hpllldir=$(includedir)/hplll

include_hplll_HEADERS=  defs.h  mat.h matpe.h  ldpe.h nr-ld.cpp  \
	hlll.h hlll.cc decomp.h decomp.cc relations.h relations.cc \
	matgen.h matgen.cc decompz.h decompz.cc nullspace.h nullspace.cc \
  	mixed-col.h  matmixed.h  lehmer.cc matpeomp.h

#EXTRA_DIST=README.html

lib_LTLIBRARIES=libhplll.la

libhplll_la_SOURCES= hlll.cc decomp.cc decompz.cc relations.cc nullspace.cc 

libhplll_la_LIBADD= -lgmp -lmpfr -lfplll

libhplll_la_LDFLAGS=-no-undefined -version-info 1:0:1


