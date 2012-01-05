// Includes
#include <stdio.h>
#include <cutil_inline.h>
#include <shrQATest.h>

#include "dataLoader_Iris.cuh"
#include "dataLoader_Test.cuh"
#include "dataLoader_Wine.cuh"
#include "dataLoader_Cancer.cuh"
#include "distanceCalculator.cuh"
#include "clustering.cuh"
#include "errors.cuh"
#include <time.h>

#define kReportsFileName "_CUDAReportsFile.txt"
#define kReportsFileNameXls "_CUDAReportsFileXls.txt"

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
	SetupTestLoader();
	SetupWineLoader();
	SetupCancerLoader();
	
	//= Load data
	//--------------------
	DataStore dataStore; // 0 - wine
	err = GetCalculatedDistances( 3, &dataStore );

	unsigned int cNeighbours = 0;
    unsigned int cClusters = 0;
    unsigned int cMedoids = 0;
    unsigned int cPopSize = 0;
    unsigned int cSteps = 0;
    unsigned int cRepeat = 0;

	unsigned int stepsMedoids;
	unsigned int stepsClusters;
	unsigned int stepsNeighbours;
		
	/*
	 medoids: 42 clusters: 8 neighbours: 29
 popSize: 256 steps: 1002
	*/
    char stateSaved = 1;
    unsigned int sPopSize = 256;
    unsigned int sSteps = 1002;
	/*
	char stateSaved = 0;
    unsigned int sNeighbours = 0;
    unsigned int sClusters = 0;
    unsigned int sMedoids = 0;
    unsigned int sPopSize = 0;
    unsigned int sSteps = 0;
	*/
	FILE * reportsFile = NULL;
	char * reportsFileName;
	char * xlsReportsFileName;
	unsigned int fileNameLength = strlen( kReportsFileName ) + strlen( dataStore.info.name );
	reportsFileName = (char*)malloc( fileNameLength + 1 );
	sprintf( reportsFileName, "%s%s", dataStore.info.name, kReportsFileName );
	fileNameLength = strlen( kReportsFileNameXls ) + strlen( dataStore.info.name );
	xlsReportsFileName = (char*)malloc( fileNameLength + 1 );
	sprintf( xlsReportsFileName, "%s%s", dataStore.info.name, kReportsFileNameXls );
	
	if ( reportsFile == NULL && !stateSaved ) {
		reportsFile = fopen( reportsFileName, "w" );
	}

	if (reportsFile != NULL ) {
		fclose( reportsFile );
		reportsFile = NULL;
	}

	if ( !stateSaved ) {
		FILE * xlsReportFile = fopen( xlsReportsFileName, "w" );
		if ( xlsReportFile != NULL ) {
			fclose( xlsReportFile );
			xlsReportFile = NULL;
		}
	}


	// calculate how big changes per step
	stepsMedoids = ( ( dataStore.info.numEntries / 4 - 3 ) / 3 );
	if ( stepsMedoids == 0 ) stepsMedoids = 1;
	
	stepsNeighbours = ( kMaxNeighbours - 1 ) / 2;
	if ( stepsNeighbours == 0 ) stepsNeighbours = 1;

	time_t currTime = 0;
	float timeDiff = 0.0;


	if ( err == errOk ) {
		for ( cPopSize = 4; cPopSize <= 256; cPopSize *= 4 ) { // 4 - 16 - 64 - 256
			// now the evolution params                    
            for ( cSteps = 2; cSteps <= 1502; cSteps += 500 ) { // 2 - 502 - 1002
				// medoids <1; numEntries/2>
				if ( stateSaved ) {
                    cPopSize = sPopSize;
                    cSteps = sSteps;
                    stateSaved = 0;
                }
				CleanAlgResults( results );

				if ( reportsFile == NULL ) {
					reportsFile = fopen( reportsFileName, "a" );
				}

				if (reportsFile != NULL ) {
					fprintf( reportsFile, "---------------------------------\n" );
					printf( "---------------------------------\n" );
					fprintf( reportsFile, " medoids: %d clusters: %d neighbours: %d\n", MEDOID_VECTOR_SIZE, MAX_CLUSTER_SIZE, kMaxNeighboursToUSe );
					printf( " medoids: %d clusters: %d neighbours: %d\n", cMedoids, cClusters, cNeighbours );
					fprintf( reportsFile, " popSize: %d steps: %d\n", cPopSize, cSteps );
					printf( " popSize: %d steps: %d\n", cPopSize, cSteps );
					fprintf( reportsFile, " Results:\n" );
					printf( " Results:\n" );
					fclose( reportsFile );
					reportsFile = NULL;
				}

                for ( cRepeat = 0; cRepeat < 5; cRepeat++ ) {
                    time( &currTime );
					err = runClustering( cPopSize, cSteps, &dataStore, &results );
					timeDiff = difftime( time( NULL ), currTime );

					if ( results.time.min == 0 || results.time.min > timeDiff ) {
						results.time.min = timeDiff;
					}
					if ( results.time.max == 0 || results.time.max < timeDiff ) {
						results.time.max = timeDiff;
					}
					results.time.sum += timeDiff;
                                
                    if ( err != errOk ) {
                        break;
                    }
                }
				results.time.mean = results.time.sum / 5.0;
				printf( "=============================================\n" );
                if ( err != errOk ) {
                    break;
                }
                            
				if ( reportsFile == NULL ) {
					reportsFile = fopen( reportsFileName, "a" );
				}

				if (reportsFile != NULL ) {
					fprintf( reportsFile, " BDI:  %f / %f / %f\n", results.bdi.min, results.bdi.mean, results.bdi.max );
					printf( " BDI:  %f / %f / %f\n", results.bdi.min, results.bdi.mean, results.bdi.max );
					fprintf( reportsFile, " DI:  %f / %f / %f\n", results.di.min, results.di.mean, results.di.max );
					printf( " DI:  %f / %f / %f\n", results.di.min, results.di.mean, results.di.max );
					fprintf( reportsFile, " Rand:  %f / %f / %f\n", results.rand.min, results.rand.mean, results.rand.max );
					printf( " Rand:  %f / %f / %f\n", results.rand.min, results.rand.mean, results.rand.max );
					fprintf( reportsFile, " Time:  %f / %f / %f\n\n", results.time.min, results.time.mean, results.time.max );
					printf( " Time:  %f / %f / %f\n\n", results.time.min, results.time.mean, results.time.max );
					fclose( reportsFile );
					reportsFile = NULL;

				} else {
					// Filed to write report
				}

				// xls readable file
				{
					// medoids, clusters, neighbours, popSize, Steps, BDI_min, BDI_mean, BDI_max, DI_min, DI_mean, DI_max, Rand_min, Rand_mean, Rand_max, Time_min, Time_mean, Time_max
					FILE * xlsReportFile = fopen( xlsReportsFileName, "a" );
					if ( xlsReportFile != NULL ) {
						fprintf( xlsReportFile, "%u, %u, %u, %u, %u, ",
							cMedoids, cClusters, cNeighbours, cPopSize, cSteps );
						fprintf( xlsReportFile, "%f, %f, %f, ",
							results.bdi.min, results.bdi.mean, results.bdi.max );
						fprintf( xlsReportFile, "%f, %f, %f, ",
							results.di.min, results.di.mean, results.di.max );
						fprintf( xlsReportFile, "%f, %f, %f, ",
							results.rand.min, results.rand.mean, results.rand.max );
						fprintf( xlsReportFile, "%f, %f, %f\n", results.time.min, results.time.mean, results.time.max );
						fclose( xlsReportFile );
						xlsReportFile = NULL;
					}
				}
            } // for steps
            if ( err != errOk ) {
                break;
            }
        } // pop size
	} // error
	//====
	/*
	if ( err == errOk ) {
		for (int i = 0; i < repeats; i++ ) {
			time( &currTime );
			err = runClustering( popSize, evoSteps, &dataStore, &results );
			timeDiff = timeDiff( time( NULL ), currTime );

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
	*/
	// TODO:
	//releaseDistances();

	shrQAFinishExit(argc, (const char **)argv, QA_PASSED);
	return 0;
}