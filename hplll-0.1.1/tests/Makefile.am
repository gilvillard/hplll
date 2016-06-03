
AUTOMAKE_OPTIONS = gnu

EXTRA_DIST= 416_in 2_in 3_in C1_out  C2_out C3_in C3_out  \
             C_huge_in C_huge_out r128int collection 

AM_CPPFLAGS = -I$(top_srcdir)/src

#nullspacetest : mettre à jour avec hjls.h 


#INCLUDES = -I$(top_srcdir)/src 

#Mis déjà dans le configure 
#LDADD =  -lfplll  


TESTS = arithmetic  hllltest relationstest matint128 mpz128inttest

arithmetic_SOURCES=arithmetic.cc  

hllltest_SOURCES=hllltest.cc

relationstest_SOURCES=relationstest.cc

matint128_SOURCES=matint128.cc

mpz128inttest_SOURCES=mpz128inttest.cc

check_PROGRAMS = $(TESTS)



#hllltest_SOURCES=hllltest.cc  

#relationstest_SOURCES=relationstest.cc 

#relationstesthuge_SOURCES=relationstesthuge.cc 

#nullspacetest_SOURCES=nullspacetest.cc mettre à jour avec hjls.h 

#./relationstest ./nullspacetest 

EXTRA_PROGRAMS = relationstesthuge \
                 gen_table_knapsack gen_table_ajtai \
                 fhn_latticec fhn_ideal fhn_ajtai

relationstesthuge_SOURCES=relationstesthuge.cc

gen_table_knapsack_SOURCES = gen_table_knapsack.cc
gen_table_ajtai_SOURCES = gen_table_ajtai.cc

fhn_latticec_SOURCES = fhn_latticec.cc
fhn_ideal_SOURCES = fhn_ideal.cc
fhn_ajtai_SOURCES = fhn_ajtai.cc


benchmarks : gen_table_knapsack  gen_table_ajtai
	./gen_table_knapsack >& table_knapsack.log
	./gen_table_ajtai >& table_ajtai.log

benchntl : fhn_latticec  fhn_ideal fhn_ajtai
	./fhn_latticec >& fhn_latticec.log
	./fhn_ideal >& fhn_ideal.log
	./fhn_ajtai >& fhn_ajtai.log



