

#include <omp.h>
#include <iostream>
#include <chrono>

using namespace std;
using namespace chrono; 

int main(int argc, char *argv[])  {


	double tp;

        int ll; 

        tp =  omp_get_wtime();
	auto start = chrono::high_resolution_clock::now();

#pragma omp parallel for 
	for (int h=1; h<2; h++) ll+=1;

	auto finish = chrono::high_resolution_clock::now();
	tp =  omp_get_wtime() - tp;

	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time: " << elapsed.count() << endl; 
	cout << "omp time: " << tp << endl;


	return 0;
}
