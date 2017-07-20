#! /bin/sh

CXXFLAGS += -O2 -Wall -I/usr/local/include/hplll

# Here with standard installation directories

FPLLL_LIBS += -lfplll
MPFR_LIBS +=  -lmpfr
GMP_LIBS   += -lgmp

LOADLIBES +=  ${FPLLL_LIBS} ${MPFR_LIBS} ${GMP_LIBS} 


.SUFFIXES: .o .c .cc  

.cc.o:
	g++  -o $@ -c $<   ${CXXFLAGS} ${LOADLIBES} 
 
.cc:
	g++  -o $@ $<    ${CXXFLAGS}  ${LOADLIBES} 


 
