
AUTOMAKE_OPTIONS = gnu

EXTRA_DIST= 416_in 2_in 3_in C1_out  C2_out C3_in C3_out  C4_out\
             C_huge_in C_huge_out r128int collection benchmarks_results \
             perfpublisher.sh 

AM_CPPFLAGS = -I$(top_srcdir)/src

#nullspacetest : mettre à jour avec hjls.h 


#INCLUDES = -I$(top_srcdir)/src 

#Mis déjà dans le configure 
#LDADD =  -lfplll  

PERFPUBLISHERFILE=tests-report.xml

TESTS = arithmetic  hllltest relationstest matint128 mpz128inttest

arithmetic_SOURCES=arithmetic.cc  

hllltest_SOURCES=hllltest.cc

relationstest_SOURCES=relationstest.cc

matint128_SOURCES=matint128.cc

mpz128inttest_SOURCES=mpz128inttest.cc

check_PROGRAMS = $(TESTS)


#nullspacetest_SOURCES=nullspacetest.cc mettre à jour avec hjls.h 

#./relationstest ./nullspacetest 

EXTRA_PROGRAMS = relationstesthuge \
                 gen_table_knapsack gen_table_ajtai \
                 fhn_latticec fhn_ideal fhn_ajtai fhn_table_knapsack \
                 fhn_acyclic fhn_al fhn_ntru fhn_ntrul

relationstesthuge_SOURCES=relationstesthuge.cc

gen_table_knapsack_SOURCES = gen_table_knapsack.cc
gen_table_ajtai_SOURCES = gen_table_ajtai.cc

fhn_latticec_SOURCES = fhn_latticec.cc
fhn_ideal_SOURCES = fhn_ideal.cc
fhn_ajtai_SOURCES = fhn_ajtai.cc
fhn_table_knapsack_SOURCES = fhn_table_knapsack.cc
fhn_acyclic_SOURCES = fhn_acyclic.cc
fhn_al_SOURCES = fhn_al.cc
fhn_ntru_SOURCES = fhn_ntru.cc
fhn_ntrul_SOURCES = fhn_ntrul.cc

benchmarks : bench benchntl 

bench : b1 b2

b1 : gen_table_knapsack  
	./gen_table_knapsack >& benchmarks_results/table_knapsack.log

b2 : gen_table_ajtai
	./gen_table_ajtai >& benchmarks_results/table_ajtai.log


benchntl: bfhn1 bfhn2 bfhn3 bfhn4 bfhn5 bfhn6

bfhn1 : fhn_latticec  
	./fhn_latticec >& benchmarks_results/fhn_latticec.log

bfhn2 : fhn_ideal 
	./fhn_ideal >& benchmarks_results/fhn_ideal.log

bfhn3 : fhn_ajtai 
	./fhn_ajtai >& benchmarks_results/fhn_ajtai.log

bfhn4 : fhn_table_knapsack
	./fhn_table_knapsack >& benchmarks_results/fhn_table_knapsack.log

bfhn5 : fhn_acyclic 
	./fhn_acyclic >& benchmarks_results/fhn_acyclic.log

bfhn6 : fhn_al 
	./fhn_al >& benchmarks_results/fhn_al.log

bfhn7 : fhn_ntru 
	./fhn_ntru >& benchmarks_results/fhn_ntru.log

bfhn8 : fhn_ntrul 
	./fhn_ntrul >& benchmarks_results/fhn_ntrul.log


# Perfpublisher script interaction - AC 2016/06/09
perfpublisher:
	+chmod 777 ./perfpublisher.sh
	+./perfpublisher.sh "$(PERFPUBLISHERFILE)" "$(TESTS)" "$(CXX)"


