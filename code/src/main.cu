// Includes
#include <stdio.h>
#include <cutil_inline.h>
#include <shrQATest.h>

#include "dataLoader.cuh"
#include "distanceCalculator.cuh"
#include "clustering.cuh"
#include "errors.cuh"

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

	FILE * testFile = fopen("klopik.data", "w");
	fwrite("klopik", 1, 6, testFile);
	fclose(testFile);

	ErrorCode err;
	
	if ( getDistances() == 0 ) {
		err = calculateDistances()
	}	

	if ( err == errOk ) {		
		err = runClustering( 100, 10 );
	}

	releaseDistances();
	return 0;
}