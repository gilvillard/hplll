AUTOMAKE_OPTIONS = gnu

EXTRA_DIST= 416_in 2_in 3_in C1_out  C2_out C3_out  \
             C4_out C5_out C6_out C7_out 

AM_CPPFLAGS = -I$(top_srcdir)/src


check_PROGRAMS = arithmetic  hllltest relationstest nullspacetest 

noinst_PROGRAMS= main pkc  ajtai block prelimin U  ord 


#INCLUDES = -I$(top_srcdir)/src 

#Mis déjà dans le configure 
LDADD =  -lfplll -lntl #$(top_srcdir)/src/libhplll.la

main_SOURCES=main.cc 
pkc_SOURCES=issac-pkc.cc 
ajtai_SOURCES=issac-ajtai.cc 
block_SOURCES=issac-block.cc 
prelimin_SOURCES=issac-prelimin.cc 
ord_SOURCES=issac-ord.cc 
U_SOURCES=issac-U.cc 


arithmetic_SOURCES=arithmetic.cc  

hllltest_SOURCES=hllltest.cc  

relationstest_SOURCES=relationstest.cc 

nullspacetest_SOURCES=nullspacetest.cc 


check : 
	./arithmetic 
	./hllltest
	./relationstest
	./nullspacetest 


