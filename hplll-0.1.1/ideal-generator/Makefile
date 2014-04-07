#! /bin/sh


CXXFLAGS += -O2 -Wall 

FPLLL_LIBS += -lfplll
MPFR_LIBS +=  -lmpfr
GMP_LIBS   += -lgmpxx  -lgmp

LOADLIBES +=  ${FPLLL_LIBS} ${MPFR_LIBS} ${GMP_LIBS}



.SUFFIXES: .o .c .cc  

.cc.o:
	g++  -o $@ -c $<   ${CXXFLAGS}
 
.cc:
	g++  -o $@ $<   ${CXXFLAGS}  ${LOADLIBES} 

.c.o:
	gcc  -o $@ -c $<  ${CXXFLAGS}
 
.c:
	gcc  -o $@ $<  ${CXXFLAGS}  ${LOADLIBES} 

ideal: ideal.cc tools.h
	g++ -o ideal ideal.cc ${CXXFLAGS}  -lntl ${LOADLIBES}


clean:
	rm -f *~ *.o
 
