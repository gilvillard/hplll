#! /bin/sh


CXXFLAGS += -O2 -Wall 

GMP_LIBS   += -lgmpxx  -lgmp

LOADLIBES +=   ${GMP_LIBS}



.SUFFIXES: .o .c .cc  

.cc.o:
	g++  -o $@ -c $<   ${CXXFLAGS}
 
.cc:
	g++  -o $@ $<   ${CXXFLAGS}  ${LOADLIBES} 

.c.o:
	gcc  -o $@ -c $<  ${CXXFLAGS}
 
.c:
	gcc  -o $@ $<  ${CXXFLAGS}  ${LOADLIBES} 


all: test-recrypt

test-recrypt: test-recrypt.o randomStuff.o fhe-utils.o fhe-keygen.o fhe-enc.o fhe-dec.o fhe-recrypt.o
	g++ $(LDFLAGS) -o test-recrypt test-recrypt.o fhe-keygen.o fhe-utils.o randomStuff.o fhe-enc.o fhe-dec.o fhe-recrypt.o -lntl -lgmp -lm 

test-recrypt.o: test-recrypt.cc fhe-recrypt.cc fhe.h

keygen.o: keygen.cc fhe.h

fhe-utils.o: fhe-utils.cc fhe.h

fhe-keygen.o: fhe-keygen.cc fhe.h randomStuff.h

fhe-enc.o: fhe-enc.cc fhe.h randomStuff.h

fhe-dec.o: fhe-dec.cc fhe.h

fhe-recrypt.o: fhe-recrypt.cc fhe.h
