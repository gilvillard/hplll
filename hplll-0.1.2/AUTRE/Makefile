#! /bin/sh

CXXFLAGS += -O2  -std=c++11 -flto    -Wall  

CXXFLAGS += -I/usr/local/include 
LDFLAGS += -L/usr/local/lib


# TODO EX: setenv LD_LIBRARY_PATH ${HOME}/Soft/lib

# If non standard installation directories 
#  for instance, uncomment and fix the following  

#PREFIX = ${HOME}/Soft

#CXXFLAGS += -O2 -Wall -I${PREFIX}/include
#LDFLAGS += -L${PREFIX}/lib

# Here with standard installation directories
#  otherwise modify the compilation variables

FPLLL_LIBS += -lfplll
MPFR_LIBS +=  -lmpfr
GMP_LIBS   += -lgmp

LOADLIBES +=  ${FPLLL_LIBS} ${MPFR_LIBS} ${GMP_LIBS} 


.SUFFIXES: .o .c .cc  

.cc.o:
	g++  -o $@ -c $<   ${CXXFLAGS}
 
.cc:
	g++  -o $@ $<   ${CXXFLAGS}  ${LDFLAGS} ${LOADLIBES} 

 
