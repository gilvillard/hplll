#! /bin/sh

CXXFLAGS += -O3 -std=c++11 -Wall  

CXXFLAGS += -I/opt/local/include 
LDFLAGS += -L/opt/local/lib


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
	#/opt/local/bin/g++  -o $@ -c $<   ${CXXFLAGS}
	clang++ -o $@ -c $<   ${CXXFLAGS}
 
.cc:
	#/opt/local/bin/g++  -o $@ $<   ${CXXFLAGS}  ${LDFLAGS} ${LOADLIBES} 
	clang++ -o $@ $<   ${CXXFLAGS}  ${LDFLAGS} ${LOADLIBES}

 
