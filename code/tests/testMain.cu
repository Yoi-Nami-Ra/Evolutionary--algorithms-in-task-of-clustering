/**
 Created 13/02/2012
 Jaroslaw Wojtasik

CUDA testss

testMain.cu
 **/

// Includes
#include <stdio.h>
#include <cutil_inline.h>
#include <shrQATest.h>
#include "distanceCalculator_test.cuh"

// Host code
int main(int argc, char** argv)
{
    int devID;	
	cudaDeviceProp props;
	
    shrQAStart(argc, argv);

	//Check which GPU is used
	cutilChooseCudaDevice(argc, argv);
	
	//Get GPU information
	cutilSafeCall(cudaGetDevice(&devID));
	cutilSafeCall(cudaGetDeviceProperties(&props, devID));
	printf("Device %d: \"%s\" with Compute %d.%d capability\n", 
			devID, props.name, props.major, props.minor);

	runDistancesTests();

	shrQAFinishExit(argc, (const char **)argv, QA_PASSED);
	return 0;
}