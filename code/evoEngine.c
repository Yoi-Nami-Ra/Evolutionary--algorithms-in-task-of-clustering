/**
 17/11/2011
 Jaroslaw Wojtasik

 noCuda

 evoEngine.c
 **/

#include "evoEngine.h"
#include "dataLoader.h"
#include "dataLoader_Iris.h"
#include "dataLoader_Test.h"
#include "dataLoader_Wine.h"
#include "dataLoader_Cancer.h"
#include "distanceCalculator.h"
#include "clustering.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define kReportsFileName "_reportsFile.txt"
#define kReportsFileNameXls "_reportsFileXls.txt"
const unsigned char kMaxNeighbours = MAX_NEIGHBOURS;

/**
 * Displays list of available loaders.
 * @return (unsigned int) number of elements on the list
 */
unsigned char DisplayLoadersList( void );
unsigned char DisplayLoadersList( void ) {
	ErrorCode err = errOk;
	unsigned char num = 0;
	char ** list = NULL;
	int i = 0;

	// 1) Get the list
	err = ListAllFunctions( &num, &list );
	if ( err != errOk ) {
		reportError( err, "Error while fetching list of loaders%s.", "" );
	}
	// 2) Print
	for ( i = 0; i < num; i++ ) {
		printf( " %d : %s\n", i+1, list[ i] );
	}

	return num;
}

/**
 * Ask user for selection
 */
unsigned int PromptForSelection( void );
unsigned int PromptForSelection( void ) {
	unsigned int result;

	printf( "\n >> ");
	scanf( "%i", &result );

	return result;
}

/**
 * Print Details
 */
void PrintDetails( unsigned int selected );
void PrintDetails( unsigned int selected ) {
	DataStore dataStore;
	const char * loaderName;

	LoaderDetails( selected - 1, &dataStore );
	loaderName = LoaderName( selected-1, NULL );

	printf( "Loader %i\n Name: %s\n Num Entries: %i\n Dimensions: %i\n",
		selected, loaderName, dataStore.info.numEntries, dataStore.info.dimSize );
}

/**
 * Confirm selection of this loader.
 * 
 * @return True if the user selected it.
 */
char ConfirmSelection( void );
char ConfirmSelection( void ) {
	int result;

	printf( "\n [y,n]>> ");
	do {
		result = getchar();
	} while ( result != 'y' && result != 'n' );

	return ( result == 'y' );
}

void runEvo( void ) {
	ErrorCode err = errOk;
	DataStore dataStore;
	EvolutionProps props;
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
	char stateSaved = 1;
    unsigned int sNeighbours = 15;
    unsigned int sClusters = 15;
    unsigned int sMedoids = 42;
    unsigned int sPopSize = 256;
    unsigned int sSteps = 1002;
	*/
	char stateSaved = 0;
    unsigned int sNeighbours = 0;
    unsigned int sClusters = 0;
    unsigned int sMedoids = 0;
    unsigned int sPopSize = 0;
    unsigned int sSteps = 0;
	
	double diffTime = 0.0;
	time_t currTime = 0;
	double minTime = 0.0;
	double maxTime = 0.0;
	double meanTime = 0.0;
	double sumTime = 0.0;

	FILE * reportsFile = NULL;
    char * reportsFileName;
	char * xlsReportsFileName;
    unsigned int fileNameLength;

	// -- Setting up all loaders
	SetupIrisLoader();
	SetupTestLoader();
	SetupWineLoader();
	SetupCancerLoader();
    
    // hardcoded selection: 2 - Wine
    err = GetCalculatedDistances( 2, &dataStore );

	
	fileNameLength = (unsigned int)strlen( kReportsFileName ) + (unsigned int)strlen( dataStore.info.name );
	reportsFileName = (char*)malloc( fileNameLength + 1 );
	sprintf( reportsFileName, "%s%s", dataStore.info.name, kReportsFileName );
	fileNameLength = (unsigned int)strlen( kReportsFileNameXls ) + (unsigned int)strlen( dataStore.info.name );
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
    
    if ( err != errOk ) {
        // error occured can't continue with the algorithms
        printf( " Error occured while preparing data for algorithms" );
    } else {
		for ( cPopSize = 4; cPopSize <= 256; cPopSize *= 4 ) { // 4 - 16 - 64 - 256
			// now the evolution params                    
            for ( cSteps = 2; cSteps <= 1002; cSteps += 500 ) { // 2 - 502 - 1002
				// medoids <1; numEntries/2>
				for ( cMedoids = 3; cMedoids <= dataStore.info.numEntries / 4; cMedoids += stepsMedoids ) {
					// cluster max size <1; medoidvectorsize>
					stepsClusters = ( cMedoids - 1) / 4;
					if ( stepsClusters == 0 ) stepsClusters = 1;
					for ( cClusters = 1; cClusters <= cMedoids; cClusters += stepsClusters ) {
						// max neighbours <1; hardMax>
						for ( cNeighbours = 1; cNeighbours <= kMaxNeighbours; cNeighbours += stepsNeighbours ) {
							if ( stateSaved ) {
                                cNeighbours = sNeighbours;
                                cClusters = sClusters;
                                cMedoids = sMedoids;
                                cPopSize = sPopSize;
                                cSteps = sSteps;
                                stateSaved = 0;
                            }
							DefaultProps( &props, &dataStore );
							props.evoSteps = cSteps;
                            props.popSize = cPopSize;
                            props.medoidsVectorSize = cMedoids;
                            props.maxNeighbours = cNeighbours;
                            props.maxClusterSize = cClusters;
							ConfigureAlgorithms( &props );
                            
							minTime = 0.0;
							maxTime = 0.0;
							sumTime = 0.0;

							if ( reportsFile == NULL ) {
								reportsFile = fopen( reportsFileName, "a" );
							}

							if (reportsFile != NULL ) {
								fprintf( reportsFile, "---------------------------------\n" );
								fprintf( reportsFile, " medoids: %d clusters: %d neighbours: %d\n", cMedoids, cClusters, cNeighbours );
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
                                err = RunClustering( &props );
								diffTime = difftime( time( NULL ), currTime );
								sumTime += diffTime;
								if ( minTime == 0.0 || minTime > diffTime ) {
									minTime = diffTime;
								}
								if ( maxTime == 0.0 || maxTime < diffTime ) {
									maxTime = diffTime;
								}
                                
                                if ( err != errOk ) {
                                    break;
                                }
                            }
							meanTime = sumTime / 5.0;
							ClearProps( &props );
							printf( "=============================================\n" );
                            if ( err != errOk ) {
                                break;
                            }
                            
							if ( reportsFile == NULL ) {
								reportsFile = fopen( reportsFileName, "a" );
							}

							if (reportsFile != NULL ) {
								fprintf( reportsFile, " BDI:  %f / %f / %f\n", props.resultBDI.min, props.resultBDI.mean, props.resultBDI.max );
								printf( " BDI:  %f / %f / %f\n", props.resultBDI.min, props.resultBDI.mean, props.resultBDI.max );
								fprintf( reportsFile, " DI:  %f / %f / %f\n", props.resultDI.min, props.resultDI.mean, props.resultDI.max );
								printf( " DI:  %f / %f / %f\n", props.resultDI.min, props.resultDI.mean, props.resultDI.max );
								fprintf( reportsFile, " Rand:  %f / %f / %f\n", props.resultRand.min, props.resultRand.mean, props.resultRand.max );
								printf( " Rand:  %f / %f / %f\n", props.resultRand.min, props.resultRand.mean, props.resultRand.max );
								fprintf( reportsFile, " Time:  %f / %f / %f\n\n", minTime, meanTime, maxTime );
								printf( " Time:  %f / %f / %f\n\n", minTime, meanTime, maxTime );
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
										props.resultBDI.min, props.resultBDI.mean, props.resultBDI.max );
									fprintf( xlsReportFile, "%f, %f, %f, ",
										props.resultDI.min, props.resultDI.mean, props.resultDI.max );
									fprintf( xlsReportFile, "%f, %f, %f, ",
										props.resultRand.min, props.resultRand.mean, props.resultRand.max );
									fprintf( xlsReportFile, "%f, %f, %f\n", minTime, meanTime, maxTime );
									fclose( xlsReportFile );
									xlsReportFile = NULL;
								}
							}
                        }
                        if ( err != errOk ) {
                            break;
                        }
                        if ( err != errOk ) {
                            break;
                        }
                    }
                    if ( err != errOk ) {
                        break;
                    }
                }
                if ( err != errOk ) {
                    break;
                }
            }
            if ( err != errOk ) {
                break;
            }
        } // pop size
    }


	// -- 
	return;
}