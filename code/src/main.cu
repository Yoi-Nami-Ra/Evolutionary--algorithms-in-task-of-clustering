// Includes
#include <stdio.h>
#include <cutil_inline.h>
#include <shrQATest.h>

#include "dataLoader_Iris.cuh"
#include "distanceCalculator.cuh"
#include "clustering.cuh"
#include "errors.cuh"
#include <time.h>

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

	algResults results;

	//= Setup loaders
	//--------------------
	SetupIrisLoader();
	
	//= Load data
	//--------------------
	DataStore dataStore;
	err = GetCalculatedDistances( 0, &dataStore );

	unsigned int popSize = 256;
	unsigned int evoSteps = 502;
	unsigned int repeats = 5;

	CleanAlgResults( results );

	time_t currTime = 0;
	float timeDiff = 0.0;

	if ( err == errOk ) {
		for (int i = 0; i < repeats; i++ ) {
			time( &currTime );
			err = runClustering( popSize, evoSteps, &dataStore, &results );
			timeDiff = difftime( time( NULL ), currTime );

			if ( results.time.min == 0 || results.time.min > timeDiff ) {
				results.time.min = timeDiff;
			}
			if ( results.time.max == 0 || results.time.max < timeDiff ) {
				results.time.max = timeDiff;
			}
			results.time.sum += timeDiff;
		} // for

		results.rand.mean = results.rand.sum / (float)( repeats * popSize );
		results.bdi.mean = results.bdi.sum / (float)( repeats * popSize );
		results.di.mean = results.di.sum / (float)( repeats * popSize );
		results.time.mean = results.time.sum / (float)( repeats );
		results.clusters.mean = results.clusters.sum / (float)( repeats * popSize );
		
		printf( " \n\n====\n pop(%d), steps(%d), repeats(%d)\n Rand(%f) BDI(%f) DI(%f) time(%f) clusters(%f)",
			popSize, evoSteps, repeats, results.rand.mean, results.bdi.mean, results.di.mean, results.time.mean, results.clusters.mean );

		printf( "---------------------------------\n" );
		printf( " medoids: %d clusters: %d neighbours: %d\n", MEDOID_VECTOR_SIZE, MAX_CLUSTER_SIZE, kMaxNeighboursToUSe );
		printf( " popSize: %d steps: %d\n", popSize, evoSteps );
		printf( " Results:\n" );
		printf( "=============================================\n" );
		printf( " BDI:  %f / %f / %f\n", results.bdi.min, results.bdi.mean, results.bdi.max );
		printf( " DI:  %f / %f / %f\n", results.di.min, results.di.mean, results.di.max );
		printf( " Rand:  %f / %f / %f\n", results.rand.min, results.rand.mean, results.rand.max );
		printf( " Time:  %f / %f / %f\n\n", results.time.min, results.time.mean, results.time.max );

	}

	// TODO:
	//releaseDistances();

	shrQAFinishExit(argc, (const char **)argv, QA_PASSED);
	return 0;
}