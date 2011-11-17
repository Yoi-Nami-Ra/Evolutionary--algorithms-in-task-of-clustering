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

	ErrorCode err = errOk;
	
	if ( getDistances() == 0 ) {
		err = calculateDistances();
	}	

	algResults results;
	float randMax, randMin, randMean, randSum = 0;
	float bdiMax, bdiMin, bdiMean, bdiSum = 0;
	float diMax, diMin, diMean, diSum = 0;
	float timeMax, timeMin, timeMean, timeSum = 0;
	float clusterMax, clusterMin, clusterMean, clusterSum = 0;

	unsigned int popSize = 200;
	unsigned int evoSteps = 2000;
	unsigned int repeats = 5;

	if ( err == errOk ) {
		for (int i = 0; i < repeats; i++ ) {
			err = runClustering( popSize, evoSteps, &results );
			if ( i == 0 ) {
				randSum += randMax = randMin = results.rand;
				bdiSum += bdiMax = bdiMin = results.bdi;
				diSum += diMax = diMin = results.di;
				timeSum += timeMax = timeMin = results.time;
				clusterSum += clusterMax = clusterMin = results.k;
			} else {
				if ( randMax < results.rand ) randMax = results.rand;
				if ( randMin > results.rand ) randMin = results.rand;
				randSum += results.rand;
				//-
				if ( bdiMax < results.bdi ) bdiMax = results.bdi;
				if ( bdiMin > results.bdi ) bdiMin = results.bdi;
				bdiSum += results.bdi;
				//-
				if ( diMax < results.di ) diMax = results.di;
				if ( diMin > results.di ) diMin = results.di;
				diSum += results.di;
				//-
				if ( timeMax < results.time ) timeMax = results.time;
				if ( timeMin > results.time ) timeMin = results.time;
				timeSum += results.time;
				//-
				if ( clusterMax < results.k ) clusterMax = results.k;
				if ( clusterMin > results.k ) clusterMin = results.k;
				clusterSum += results.k;
				//-
			}			
		} // for

		randMean = randSum / (float)repeats;
		bdiMean = bdiSum / (float)repeats;
		diMean = diSum / (float)repeats;
		timeMean = timeSum / (float)repeats;
		clusterMean = clusterSum / (float)repeats;

		printf( " \n\n====\n pop(%d), steps(%d), repeats(%d)\n Rand(%f) BDI(%f) DI(%f) time(%f) clusters(%f)",
			popSize, evoSteps, repeats, randMean, bdiMean, diMean, timeMean, clusterMean );

	}

	releaseDistances();

	shrQAFinishExit(argc, (const char **)argv, QA_PASSED);
	return 0;
}