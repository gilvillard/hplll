#! /bin/sh

CXXFLAGS += -O2 -Wall -L${HOME}/Soft/lib

# Here with standard installation directories

FPLLL_LIBS += -lfplll
MPFR_LIBS +=  -lmpfr
GMP_LIBS   += -lgmp

LOADLIBES +=  ${FPLLL_LIBS} ${MPFR_LIBS} ${GMP_LIBS} 

INC += ${HOME}/Soft/include 


.SUFFIXES: .o .c .cc  

.cc.o:
	g++  -o $@ -c $<   ${CXXFLAGS} ${LOADLIBES} 
 
.cc:
	g++  -o $@ $<   -I${INC} ${CXXFLAGS}  ${LOADLIBES} 


 
