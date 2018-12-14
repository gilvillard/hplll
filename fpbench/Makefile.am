AUTOMAKE_OPTIONS = gnu

EXTRA_DIST= fpbench.h  nr_FP_float128.inl

AM_CPPFLAGS = -I$(top_srcdir)/src

EXTRA_PROGRAMS = fpbench 

fpbench_SOURCES = fpbench.cc fpbench.h

fp : fpbench 
	fpbench 

