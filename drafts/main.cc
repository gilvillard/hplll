

#include <omp.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])  {


	double tp;

        int ll; 

        tp =  omp_get_wtime();

#pragma omp parallel for 
	for (int h=1; h<2; h++) ll+=1;

	tp =  omp_get_wtime() - tp;

	cout << "ptime: " << tp << endl;


	return 0;
}
