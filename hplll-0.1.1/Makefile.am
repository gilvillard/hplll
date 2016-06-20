AUTOMAKE_OPTIONS = gnu

ACLOCAL_AMFLAGS= -I macros

SUBDIRS= src  tests macros examples drafts 

nobase_dist_doc_DATA = COPYING.LESSER 

include_hpllldir=$(includedir)/hplll

AM_CXXFLAGS = -I $(srcdir)

perfpublisher : tests/perfpublisher

tests/perfpublisher:
	(cd tests; ${MAKE} perfpublisher)

#EXTRA_DIST=README.html
