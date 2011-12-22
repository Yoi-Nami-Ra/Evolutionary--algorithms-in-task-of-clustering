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
#include "distanceCalculator.h"
#include "clustering.h"
#include <stdlib.h>
#include <stdio.h>

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

	LoaderDetails( selected-1, &dataStore );
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
    
    char stateSaved = 0;
    unsigned int sNeighbours = 0;
    unsigned int sClusters = 0;
    unsigned int sMedoids = 0;
    unsigned int sPopSize = 0;
    unsigned int sSteps = 0;
    unsigned int sRepeat = 0;

	// -- Setting up all loaders
	SetupIrisLoader();
	SetupTestLoader();
	SetupWineLoader();
    
    // hardcoded selection: 0 - Iris
    err = GetCalculatedDistances( 0, &dataStore );
    
    if ( err != errOk ) {
        // error occured can't continue with the algorithms
        printf( " Error occured while preparing data for algorithms" );
    } else {
        // medoids <1; numEntries/2>
        for ( cMedoids = 1; cMedoids <= dataStore.info.numEntries / 2; cMedoids++ ) {
            // cluster max size <1; medoidvectorsize>
            for ( cClusters = 1; cClusters <= cMedoids; cClusters++ ) {
                // max neighbours <1; hardMax>
                for ( cNeighbours = 1; cNeighbours <= kMaxNeighbours; cNeighbours++ ) {
                    // now the evolution params
                    for ( cPopSize = 4; cPopSize <= 256; cPopSize += 4 ) {
                        for ( cSteps = 1; cSteps <= 2000; cSteps += 100 ) {
                            if ( stateSaved ) {
                                cNeighbours = sNeighbours;
                                cClusters = sClusters;
                                cMedoids = sMedoids;
                                cPopSize = sPopSize;
                                cSteps = sSteps;
                                cRepeat = sRepeat;
                                stateSaved = 0;
                            }
							DefaultProps( &props, &dataStore );
							props.evoSteps = cSteps;
                            props.popSize = cPopSize;
                            props.medoidsVectorSize = cMedoids;
                            props.maxNeighbours = cNeighbours;
                            props.maxClusterSize = cClusters;
							ConfigureAlgorithms( &props );

                            for ( cRepeat = 0; cRepeat < 5; cRepeat++ ) {
                                                                
                                err = RunClustering( &props );
                                
                                // TODO: count results from repeats here
                            }
                            // TODO: Display/print results
                        }
                    }
                }
            }
        }
    }


	// -- 
	return;
}