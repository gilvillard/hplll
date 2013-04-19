AUTOMAKE_OPTIONS = gnu

EXTRA_DIST= 416_in 2_in 3_in C1_out  C2_out C3_out  \
             C4_out C5_out C6_out C7_out 

AM_CPPFLAGS = -I$(top_srcdir)/src


bin_PROGRAMS= main arithmetic hllltest relationstest nullspacetest 
#check_PROGRAMS = arithmetic  hlll-test


#INCLUDES = -I$(top_srcdir)/src 

#Mis déjà dans le configure 
LDADD =  -lfplll $(top_srcdir)/src/libhplll.la

main_SOURCES=main.cc 

arithmetic_SOURCES=arithmetic.cc  

hllltest_SOURCES=hllltest.cc  

relationstest_SOURCES=relationstest.cc 

nullspacetest_SOURCES=nullspacetest.cc 


check : 
	./arithmetic 
	./hllltest
	./relationstest
	./nullspacetest 


