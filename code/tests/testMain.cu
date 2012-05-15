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
#include "errors.cuh"
#include "distanceCalculator_test.cuh"
#include "clustering_test.cuh"

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

	char * logFile = strdup( argv[ 0] );
	logFile[ strlen( logFile ) -3] = 'l';
	logFile[ strlen( logFile ) -2] = 'o';
	logFile[ strlen( logFile ) -1] = 'g';
	SetLogFile( logFile );
	free( logFile );

	logMessage(" ==== Starting Tests ====");
	runDistancesTests();
	runClusteringTests();
	logMessage(" ==== Tests Finished ====\n");

	shrQAFinishExit(argc, (const char **)argv, QA_PASSED);
	return 0;
}