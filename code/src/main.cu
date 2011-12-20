// Includes
#include <stdio.h>
#include <cutil_inline.h>
#include <shrQATest.h>

#include "dataLoader_Iris.cuh"
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

	ErrorCode err = errOk;
	
	if ( getDistances() == 0 ) {
		err = calculateDistances();
	}	

	algResults results;
	results.bdi.max = results.bdi.min = results.bdi.mean = results.bdi.sum = 0.0;
	results.clusters.max = results.clusters.min = results.clusters.mean = results.clusters.sum = 0.0;
	results.di.max = results.di.mean = results.di.min = results.di.sum = 0.0;
	results.rand.max = results.rand.mean = results.rand.min = results.rand.sum = 0.0;
	results.time.max = results.time.mean = results.time.min = results.time.sum = 0.0;

	//= Setup loaders
	//--------------------
	SetupIrisLoader();
	
	//= Load data
	//--------------------
	DataStore dataStore;
	err = GetCalculatedDistances( 1, &dataStore );

	unsigned int popSize = 4;
	unsigned int evoSteps = 2;
	unsigned int repeats = 5;

	if ( err == errOk ) {
		for (int i = 0; i < repeats; i++ ) {
			err = runClustering( popSize, evoSteps, &results );			
		} // for

		results.rand.mean = results.rand.sum / (float)repeats;
		results.bdi.mean = results.bdi.sum / (float)repeats;
		results.di.mean = results.di.sum / (float)repeats;
		results.time.mean = results.time.sum / (float)repeats;
		results.clusters.mean = results.clusters.sum / (float)repeats;
		
		printf( " \n\n====\n pop(%d), steps(%d), repeats(%d)\n Rand(%f) BDI(%f) DI(%f) time(%f) clusters(%f)",
			popSize, evoSteps, repeats, results.rand.mean, results.bdi.mean, results.di.mean, results.time.mean, results.clusters.mean );

	}

	releaseDistances();

	shrQAFinishExit(argc, (const char **)argv, QA_PASSED);
	return 0;
}