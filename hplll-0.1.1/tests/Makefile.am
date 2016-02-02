
AUTOMAKE_OPTIONS = gnu

EXTRA_DIST= 416_in 2_in 3_in C1_out  C2_out C3_in C3_out  \
             C_huge_in C_huge_out 

AM_CPPFLAGS = -I$(top_srcdir)/src

#nullspacetest : mettre à jour avec hjls.h 


#INCLUDES = -I$(top_srcdir)/src 

#Mis déjà dans le configure 
#LDADD =  -lfplll  


TESTS = arithmetic  hllltest relationstest #relationstesthuge

arithmetic_SOURCES=arithmetic.cc  

hllltest_SOURCES=hllltest.cc

relationstest_SOURCES=relationstest.cc

check_PROGRAMS = $(TESTS)



#hllltest_SOURCES=hllltest.cc  

#relationstest_SOURCES=relationstest.cc 

#relationstesthuge_SOURCES=relationstesthuge.cc 

#nullspacetest_SOURCES=nullspacetest.cc mettre à jour avec hjls.h 

#./relationstest ./nullspacetest 

EXTRA_PROGRAMS = gentable1 gentable2 gentable3 gentable4 gentable5 gentable6

gentable1_SOURCES=gentable1.cc 
gentable2_SOURCES=gentable2.cc
gentable3_SOURCES=gentable3.cc
gentable4_SOURCES=gentable4.cc
gentable5_SOURCES=gentable5.cc
gentable6_SOURCES=gentable6.cc

#check : 
#	./arithmetic 
#	./hllltest
#	./relationstest



