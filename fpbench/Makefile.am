AUTOMAKE_OPTIONS = gnu

EXTRA_DIST= fpbench.h  nr_FP_float128.inl matbench.h

AM_CPPFLAGS = -I$(top_srcdir)/src

EXTRA_PROGRAMS = fpbench matbench 

fpbench_SOURCES = fpbench.cc fpbench.h
matbench_SOURCES = matbench.cc matbench.h

fp : fpbench 
	fpbench 

mat : matbench 
	matbench
