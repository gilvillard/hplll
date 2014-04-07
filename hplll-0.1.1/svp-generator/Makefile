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

chlg: chlg.cc tools.h
	g++ -o chlg chlg.cc ${CXXFLAGS}  -lntl ${LOADLIBES}

generate_random:   generate_random.cc  tools.h
	g++ -o generate_random generate_random.cc   ${CXXFLAGS}  -lntl ${LOADLIBES} 

clean:
	rm -f *~ *.o
 
